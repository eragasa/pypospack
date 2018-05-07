from collections import OrderedDict
from pypospack.pyposmat.engines import PyposmatEngine

class PyposmatFileSampler(PyposmatEngine):

    def __init__(self,
            filename_in = 'pypospack.config.in',
            filename_out = 'pypospack.results.out'):

        # attributes which are marshalled
        self.configuration = None
        self.parameter_names = None
        self.pyposmat_filename_in = filename_in
        self.pyposmat_filename_out = filename_out

        # attributes which are not marshalled
        self.workflow = None
        self.tasks = None
        # attributes which are not marshalled
        self.pyposmat_file_in = None
        self.pyposmat_file_out = None
        self.parameter_df = None

    def add_task(self,task_name,task_directory,structure_filename,task_type):
        if self.tasks is None:
            self.tasks = OrderedDict()

        if task_type == 'gulp_phonon':
            from pypospack.task.gulp import GulpPhononCalculation

            _task_name = task_name
            _task_directory = task_directory
            _structure_filename = structure_filename

            self.tasks[task_name] = GulpPhononCalculation(
                    task_name=_task_name,
                    task_directory=_task_directory,
                    structure_filename=_structure_filename,
                    restart=False)
        elif task_type == 'gulp_gamma_phonons':
            from pypospack.task.gulp import GulpGammaPointPhonons
            _task_name = task_name
            _task_directory = task_directory
            _structure_filename = structure_filename

            self.tasks[task_name] = GulpGammaPointPhonons(
                    task_name=_task_name,
                    task_directory=_task_directory,
                    structure_filename=_structure_filename,
                    restart=False)

    def read_pyposmat_datafile(self,filename=None):
        if filename is not None:
            self.pyposmat_filename_in = filename

        self.pyposmat_file_in = PyposmatDataFile(
                filename=self.pyposmat_filename_in)
        self.pyposmat_file_in.read()

    def sample_from_file(self,filename=None):
        if filename is not None:
            self.pyposmat_filename_in = filename
            self.read_pyposmat_datafile()
        if self.pyposmat_file_in is None:
            self.read_pyposmat_datafile()
        if self.pyposmat_file_out is None:
            self.pyposmat_file_out = PyposmatDataFile
        self.parameter_df = self.pyposmat_file_in.parameter_df

        _parameter_df = self.pyposmat_file_in.parameter_df
        _columns = _parameter_df.columns
        self.parameter_names = list(_columns)
        self.qoi_names = ['MgO_NaCl.ph_{}'.format(i) for i in range(1,6+1)]
        self.file_out = open(self.pyposmat_filename_out,'w')

        str_out = ",".join(
            ['sim_id']
            +self.parameter_names 
            +self.qoi_names
        )
        self.file_out.write(str_out+"\n")

        str_out = ",".join(
            ['sim_id']\
            +len(self.parameter_names)*['param']\
            +len(self.qoi_names)*['qoi'])
        self.file_out.write(str_out+"\n")

        for idx,row in _parameter_df.iterrows():
            parameters = OrderedDict([(col,row[col]) for col in _columns])
            self.evaluate_parameter_set(parameters)
            for task_name,task in self.tasks.items():
                #sim_id self.pyposmat_filename_in.df.loc[idx,'sim_id']]\
                sim_id = idx
                str_out = ",".join(
                        [str(sim_id)]\
                        +[str(v) for k,v in parameters.items()]\
                        +[str(v) for k,v in task.results.items()])+"\n"
                self.file_out.write(str_out)
        self.file_out.close()
    def run(self):
        self.sample_from_file()

    def is_all_tasks_complete(self):
        all_tasks_complete = all(
                [task.status=='FINISHED' for tn,task in self.tasks.items()])
        return all_tasks_complete

    def evaluate_parameter_set(self,parameters):
        for task_name,task in self.tasks.items():
            for f in os.listdir(task.task_directory):
                os.unlink(os.path.join(task.task_directory,f))
            if os.path.exists(task.results_filename):
                os.remove(task.results_filename)
            task.potential = None
            task.update_status()
        self.configuration['parameters'] = copy.deepcopy(parameters)

        while not self.is_all_tasks_complete():
            for task_name,task in self.tasks.items():
                self.tasks[task_name].configuration = self.configuration
                self.tasks[task_name].on_update_status()
