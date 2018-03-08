import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def str__sampling_configuration(sampling_config):
    _n_iterations = sampling_config['n_iterations']
    _mc_seed = sampling_config['n_iterations']
    _sampling_table = [sampling_config[i] for i in range(_n_iterations)]

    s = []
    s.append(80*"-")
    s.append('{:^80}'.format('SAMPLING INFORMATION'))
    s.append(80*"-")
    s.append('{:^3} {:^10} {:^10} {:^50}'.format('i','type','n_samples','filename'))
    s.append('{} {} {} {}'.format(2*'-',10*'-',10*'-',50*'-'))
    for i,v in enumerate(_sampling_table):
        _type=v['type']
        _nsamples = v['n_samples']

        _filename=''
        if 'filename' in v:
            _filename =v['filename']
        s.append('{:^2} {:^10} {:^10} {:^50}'.format(i,_type,_nsamples,_filename))
    s.append('{} {} {} {}'.format(2*'-',10*'-',10*'-',50*'-'))
    s.append("n_iterations: {}".format(sampling_config['n_iterations']))
    s.append("mc_seed:      {}".format(sampling_config['mc_seed']))
    return "\n".join(s)

def str__qoi_information(qois):
    QOI_LINE_FORMAT="{:^14} {:^14} {:^14} {:^14} {:^14}"

    s = []
    s.append(80*'-')
    s.append('{:^80}'.format('QOI INFORMATION'))
    s.append(80*'-')
    s.append(QOI_LINE_FORMAT.format(
        'qoi_name','qoi_type','target','structure','structure_name'))
    s.append(QOI_LINE_FORMAT.format(14*'-',14*'-',14*'-',14*'-',14*'-'))
    for k,v in qois.items():
        qoi_name=k
        qoi_type=v['qoi_type']
        qoi_target=v['target']
        qoi_s = [ks for ks in v['structures']]
        qoi_n = [kv for ks,kv in v['structures'].items()]
        for i in range(len(qoi_s)):
            if i == 0:
                s.append(QOI_LINE_FORMAT.format(
                    qoi_name,qoi_type,qoi_target,qoi_s[i],qoi_n[i]))
            else:
                s.append(QOI_LINE_FORMAT.format(
                    '','','',qoi_s[i],qoi_n[i]))
    s.append(QOI_LINE_FORMAT.format(14*'-',14*'-',14*'-',14*'-',14*'-'))
    return '\n'.join(s)

class PyposmatIterativeDataAnalyzer(object):
    def __init__(self,src_dir):
        self.fn_config = None
        self.configuration = None

        self.iteration = []

        print(80*'-')
        self.load_configuration_file(src_dir)
        self.load_data_files(src_dir)
        print(80*'-')
        self.report()
    @property
    def n_iterations(self):
        return self.configuration.sampling_type['n_iterations']

    @property
    def qoi_targets(self):
        return self.configuration.qoi_targets

    def load_configuration_file(self,src_dir):
        self.fn_config=os.path.join(src_dir,'pyposmat.config.in')
        print('reading {}...'.format(self.fn_config))
        self.configuration=PyposmatConfigurationFile()
        self.configuration.read(filename=self.fn_config)

    def load_data_files(self,src_dir):
        for i in range(self.n_iterations):
            self.iteration.append(OrderedDict())
        self.__load_results_files(src_dir)
        self.__load_kde_files(src_dir)

    def report(self):
        self.print__sampling_configuration()
        self.print__qoi_information()
        self.print__simulation_information()

    def __load_results_files(self,src_dir):
        PYPOSPACK_RESULTS_FN_FMT="pyposmat.results.{}.out"
        for i in range(self.n_iterations):
            fn=os.path.join(src_dir,PYPOSPACK_RESULTS_FN_FMT.format(i))
            print('reading {}...'.format(fn))
            self.iteration[i]['results']=PyposmatDataAnalyzer()
            self.iteration[i]['results'].read_configuration_file(filename=self.fn_config)
            self.iteration[i]['results'].read_data_file(filename=fn)

    def __load_kde_files(self,src_dir):
        PYPOSPACK_KDE_FN_FMT="pyposmat.kde.{}.out"
        for i in range(self.n_iterations):
            fn=os.path.join(src_dir,PYPOSPACK_KDE_FN_FMT.format(i))
            print('reading {}...'.format(fn))
            self.iteration[i]['kde']=PyposmatDataAnalyzer()
            self.iteration[i]['kde'].read_configuration_file(filename=self.fn_config)
            try:
                self.iteration[i]['kde'].read_data_file(filename=fn)
            except FileNotFoundError as e:
                self.iteration[i]['kde']=None

    def print__sampling_configuration(self):
        print(str__sampling_configuration(
                self.configuration.sampling_type
            )
        )

    def print__qoi_information(self):
        print(str__qoi_information(
                self.configuration.qois
            )
        )

    def print__simulation_information(self):

        n_simulation_info = OrderedDict()
        for i in range(n_iterations):
            _total=configuration.sampling_type[i]['n_samples']
            try:
                _results,ncols=data_results[i].df.shape
                _kde,ncols=data_kde[i].df.shape
            except AttributeError as e:
                _kde=0

            n_simulation_info[i] = OrderedDict()
            n_simulation_info[i]['total']=_total
            n_simulation_info[i]['results']=_results
            n_simulation_info[i]['kde']=_kde

        s = []
        s.append('{:^5} {:^10} {:^10} {:^10}'.format(
            'i','total','results','kde'))
        s.append('{:^5} {:^10} {:^10} {:^10}'.format(
            5*'-',10*'-',10*'-',10*'-'))
        for i in range(n_iterations):
            s.append('{:^5} {:>10} {:>10} {:>10}'.format(
                i,
                n_simulation_info[i]['total'],
                n_simulation_info[i]['results'],
                n_simulation_info[i]['kde']))
        print('\n'.join(s))


if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument("--src_dir",
    #        action='store',
    #        dest='data_dir',
    #        type=str,
    #        help="location of the data directory")
    #parse_args = parser.parse_args()

    src_dir = "data__Ni__eam__born_exp_bjs_00"
    a = PyposmatIterativeDataAnalyzer(src_dir)
    exit()
