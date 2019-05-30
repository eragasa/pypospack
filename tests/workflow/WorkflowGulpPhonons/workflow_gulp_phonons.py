from collections import OrderedDict
import numpy as np
import pandas as pd

class PyposmatData(object):

    def __init__(self,data_directory):
        self.RESULTS_FILENAME_FORMAT = "results_{:03d}.out"
        self.PARETO_FILENAME_FORMAT = "pareto_{:03d}.out"
        self.CULLED_FILENAME_FORMAT = "culled_{:03d}.out"
        self._data_directory = None
        self.data_directory = data_directory
        self.n_iterations = None
        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None

        self.qoi_references = OrderedDict()
        self.parameter_references = OrderedDict()

    def read(self,data_directory=None):
        if data_directory is not None:
            self.data_directory = data_directory

        data_filenames = os.listdir(self.data_directory)
        n_culled = len([s for s in data_filenames if s.startswith('culled')])
        n_pareto = len([s for s in data_filenames if s.startswith('pareto')])
        n_results = len([s for s in data_filenames if s.startswith('results')])

        if n_culled == n_pareto and n_culled == n_results:
            self.n_iterations = n_culled
        else:
            raise ValueError("missing some files")

        self.culled = []
        self.results = []
        self.pareto = []
        for i in range(self.n_iterations):
            self.culled.append(PyposmatDataFile(
                filename=os.path.join(
                    self.data_directory,
                    self.CULLED_FILENAME_FORMAT.format(i))))
            self.results.append(PyposmatDataFile(
                filename=os.path.join(
                    self.data_directory,
                    self.RESULTS_FILENAME_FORMAT.format(i))))
            self.pareto.append(PyposmatDataFile(
                filename=os.path.join(
                    self.data_directory,
                    self.PARETO_FILENAME_FORMAT.format(i))))

        for filename in self.results: 
            filename.read()
            filename.qoi_references = self.qoi_references
        for filename in self.pareto:
            filename.read()
            filename.qoi_references = self.qoi_references
        for filename in self.culled:
            filename.read()
            filename.qoi_references = self.qoi_references
   
    @property
    def data_directory(self):
        return self._data_directory

    @data_directory.setter
    def data_directory(self,directory):
        if os.path.isdir(directory):
            self._data_directory = data_directory
        else:
            raise

class PyposmatDataFile(object):

    def __init__(self,filename):
        self.filename = filename
        self.parameter_names = None
        self.qoi_names = None
        self.err_names = None
   
        self.df = None
        self.parameter_df = None
        self.error_df = None
        self.qoi_df = None
        self.rescaled_error_df = None
        
        self.optimal_indices = None
        self.optimal_df = None
        self.optimal_parameter_df = None
        self.optimal_qoi_df = None
        self.optimal_error_df = None

    def read(self,filename=None):
        if filename is not None:
            self.filename = filename

        with open(self.filename,'r') as f:
            lines = f.readlines()

        self.names = [s.strip() for s in lines[0].strip().split(',')]
        self.types = [s.strip() for s in lines[1].strip().split(',')]
        
        self.values = []
        for i in range(2,len(lines)):
            line = lines[i].strip()
            values = [float(s.strip()) for s in line.split(',')]
            values[0] = int(values[0])
            self.values.append(list(values))
        self.values = np.array(self.values)

        self.parameter_names = [
                n for i,n in enumerate(self.names) \
                    if self.types[i] == 'param']
        self.qoi_names = [
                n for i,n in enumerate(self.names) \
                        if self.types[i] == 'qoi']
        self.err_names = [
                n for i,n in enumerate(self.names) \
                        if self.types[i] == 'err']
        
        self.df = pd.DataFrame(data=self.values,
                columns=self.names,copy=True)
        self.df.set_index('sim_id')
        self.parameter_df = self.df[self.parameter_names] 
        self.error_df = self.df[self.err_names]
        self.qoi_df = self.df[self.qoi_names]

    def create_optimal_population(
            self,
            n=1,
            scaling_factors='DFT',
            err_type='abs'):
        """
            Args:
                scaling_factors(dict): the key is the error name, the value 
                    is a scalar vaue from which the errors will be divided for 
                    the purposes of scaling.
                n(int): the number of points to return
                result(str): should be either results, pareto or culled.  
                    Default is culled.
        """
       
        # here we create the scaling factors
        if type(scaling_factors) == str:
            str_sf = scaling_factors
            self.scaling_factors = OrderedDict()
            for col in self.error_df:
                if col != 'sim_id':
                    qn = '{}.{}'.format(col.split('.')[0],col.split('.')[1])
                    self.scaling_factors[col] = self.qoi_references[str_sf][qn]
        elif isinstance(scaling_factors,dict):
            self.scaling_factors = copy.deepcopy(scaling_factors)
        else:
            raise ValueError()
        
        self.rescaled_error_df = self.error_df.copy(deep=True)
        for col in self.rescaled_error_df:
            if col != 'sim_id':
                sf = self.scaling_factors[col]
                self.rescaled_error_df[col] = self.rescaled_error_df[col].abs() / sf

        # our metric is the sum of the rescaled errors
        self.rescaled_error_df['d_metric'] = self.rescaled_error_df[
                self.rescaled_error_df.columns].sum(axis=1)

        self.optimal_indices = self.rescaled_error_df.nsmallest(n,'d_metric').index

        self.optimal_df = self.df.loc[self.optimal_indices]
        #self.optimal_error_df = self.error_df.loc[self.optimal_indices]
        #self.optimal_parameter_df = self.parameter_df.loc[self.optimal_indices]
        #self.optimal_qoi_df = self.qoi_df.loc[self.optimal_indices]

    def write_optimal_population(self,filename):
        print(','.join([n for n in self.names]))
        print(','.join([t for t in self.types]))
        for row in data.culled[9].optimal_df.iterrows():
            # rows is a tuple, unpack the tuple
            _row = [a for i,a in enumerate(row[1])]
            # cast sim_id back into an int
            _row[0] = int(_row[0])
            print(','.join([str(s) for s in _row]))
def get_minimum_population(d_vector, n = 1):
        """
        Args:
            d_vector(numpy.ndarray): a numpy array of metrics upon
                which to score upon.
            n(int): number of points closest to that vector
        return an index of the dvector closest to the utopia point
        """
        idx_shortest = np.where(d_vector = d_vector.min())
        idx_shortest_pop = np.argpartition(d_vector,n)[:n_potentials]
        return idx_shortest, idx_shortest_pop

def get_absolute_scaling_factors(sr):
    """Determines scaling factor

    If q is a quantity of interest and eta is the scaling factor, then the 
    scaled error should be q/eta ~= 1
    Returns:
        sf - (list of float) a vector of scaling factors
    """
    sf = {} 
    for k,v in sr.qoi_ref.items():
        err_name = k + ".err"
        sf[err_name] = v

    return sf

if __name__ == "__main__":
    import os
    data_directory = os.path.join('test_WorkflowGulpPhonons','data','output')
    names = [
            'chrg_Mg','chrg_O',
            'MgMg_A','MgMg_rho','MgMg_C',
            'MgO_A','MgO_rho','MgO_C',
            'OO_A','OO_rho','OO_C',
            'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
            'MgO_NaCl.B','MgO_NaCl.G',
            'MgO_NaCl.fr_a','MgO_NaCl.fr_c','MgO_NaCl.sch',
            'MgO_NaCl.001s',
            'MgO_NaCl.a0.err','MgO_NaCl.c11.err','MgO_NaCl.c12.err',
            'MgO_NaCl.c44.err','MgO_NaCl.B.err','MgO_NaCl.G.err',
            'MgO_NaCl.fr_a.err','MgO_NaCl.fr_c.err','MgO_NaCl.sch.err',
            'MgO_NaCl.001s.err']
    latex_labels = {
            'MgO_NaCl.a0':r'$a_0$',
            'MgO_NaCl.c11':r'$c_{11}$',
            'MgO_NaCl.c12':r'$c_{12}$',
            'MgO_NaCl.c44':r'$c_{44}$',
            'MgO_NaCl.B':r'$B$',
            'MgO_NaCl.G':r'$G$',
            'MgO_NaCl.fr_a':r'$E_{fr,a}$',
            'MgO_NaCl.fr_c':r'$E_{fr,c}$',
            'MgO_NaCl.sch':r'$E_{sch}$',
            'MgO_NaCl.001s':r'$\gamma_{001}$'}
    parameter_names = ['chrg_Mg','chrg_O',
            'MgMg_A','MgMg_rho','MgMg_C',
            'MgO_A','MgO_rho','MgO_C',
            'OO_A','OO_rho','OO_C']
    qoi_names = [
            'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
            'MgO_NaCl.B','MgO_NaCl.G', 'MgO_NaCl.fr_a','MgO_NaCl.fr_c',
            'MgO_NaCl.sch','MgO_NaCl.001s']
    error_names = ['{}.err'.format(s) for s in qoi_names]
    qoi_reference_dft = {
            'MgO_NaCl.a0': 4.246,
            'MgO_NaCl.c11': 277.00031,
            'MgO_NaCl.c12': 91.67016,
            'MgO_NaCl.c44': 144.00722,
            'MgO_NaCl.B': 153.4468767,
            'MgO_NaCl.G': 92.665075,
            'MgO_NaCl.fr_a': 10.9781666,
            'MgO_NaCl.fr_c': 8.98642095,
            'MgO_NaCl.sch':5.067179685,
            'MgO_NaCl.001s': 0.055950069}
    reference_LC = [2.0,-2.0,0.0,0.5,0.0,821.6,0.3242,0.0,22764.0,0.149,27.88,4.21078276561128,307.5718095128,171.135602331774,168.168424521017,216.61433805878266,68.21810359051298,9.68024989411606,9.810715180656189,5.797051474639375,0.06783817649966246,-0.035217234388720264,30.571809512799973,79.46560233177401,24.158424521016997,63.164338058782675,-24.441896409487015,-1.2977501058839405,0.8247151806561881,0.7300514746393745,0.011888176499662464]
    reference_BG1 = [2.0,-2.0,0.0,0.5,0.0,1279.69,0.29969,0.0,9547.96,0.21916,32.0,4.20923604431415,383.274119165401,169.434215310753,179.601185701851,240.71418326230233,106.91995192732399,12.419259511088967,11.869114175328832,7.198887069605007,0.08070791160146304,-0.036763955685850114,106.27411916540098,77.76421531075299,35.591185701851,87.26418326230234,14.259951927323996,1.4412595110889672,2.8831141753288314,2.131887069605007,0.02475791160146304]
    reference_BG2 = [1.7,-1.7,0.0,0.5,0.0,929.69,0.29909,0.0,4870,0.2679,77.0,4.222448,301.315822490901,150.827961179874,142.471471673523,200.990581616883,75.2439306555135,10.434727086962994,8.526633932683126,5.509135247188169,0.0692527749868838,-0.02355200000000046,24.315822490900985,59.15796117987399,-1.5385283264769782,47.540581616883,-17.416069344486502,-0.543272913037006,-0.4593660673168749,0.442135247188169,0.0133027749868838]

    data = PyposmatData(data_directory=data_directory)
    data.parameter_names = list(parameter_names)
    data.qoi_names = list(qoi_names)
    data.error_names = list(error_names)
   
    # qoi reference values
    # DFT - from density functional theory, GGA
    # LC - Lewis and Catlow
    # BG1 - Ball and Grimes, 1
    # BG2 - Ball and Grimes, 2
    data.qoi_references['DFT'] = OrderedDict()
    data.qoi_references['LC'] = OrderedDict()
    data.qoi_references['BG1'] = OrderedDict()
    data.qoi_references['BG2'] = OrderedDict()
    for qn in qoi_names:
        qn_idx = names.index(qn)
        data.qoi_references['DFT'][qn] = qoi_reference_dft[qn]
        data.qoi_references['LC'][qn] = reference_LC[qn_idx]
        data.qoi_references['BG1'][qn] = reference_BG1[qn_idx]
        data.qoi_references['BG2'][qn] = reference_BG2[qn_idx]
    
    # parameter reference values
    data.parameter_references['LC'] = OrderedDict()
    data.parameter_references['BG1'] = OrderedDict()
    data.parameter_references['BG2'] = OrderedDict()
    for pn in parameter_names:
        pn_idx = names.index(pn)
        data.parameter_references['LC'][pn] = reference_LC[pn_idx]
        data.parameter_references['BG1'][pn] = reference_BG1[pn_idx]
        data.parameter_references['BG2'][pn] = reference_BG2[pn_idx]

    data.read()
    data.culled[9].create_optimal_population(n=10)
    #print(data.culled[9].optimal_indices)
    #print(data.culled[9].optimal_error_df)
    #print(data.culled[9].optimal_parameter_df)
    #print(data.culled[9].optimal_qoi_df)
    data.culled[9].write_optimal_population(filename='optimal_10.out')

    #print('type(optimal_indices):',type(data.culled[9].optimal_indices))
    #print('type(optimal_df):',type(data.culled[9].optimal_df))
    #print(data.culled[9].optimal_df)


    #print('type(df):',type(data.culled[9].df))
    #print(data.culled[9].rescaled_error_df)
