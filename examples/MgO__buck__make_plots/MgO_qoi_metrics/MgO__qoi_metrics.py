import os
import numpy as np
import pandas as pd
import scipy.stats
from pypospack.pyposmat.visualization.plots_1d import Pyposmat1DHistogramWithDensityPlots


#from pypospack.pyposmat.visualization import PyposmatDataFileVisualization

from post_processor import PyposmatPostProcessor

class PyposmatQoiPredictionAnalyzer(PyposmatPostProcessor):

    def get_mean(self,name):
        _mu = self.df[name].mean()
        return _mu

    def get_standard_deviation(self,name):
        _sigma = self.df[name].std()
        return _sigma
   
    def get_confidence_interval(self,mu=None,sigma=None,name=None,zscore=1.96):
        if name is not None:
            _n_samples,_s = self.df.shape
            # Define local variables
            _name = name
            _mu = self.get_mean(_name)
            _sigma = self.get_standard_deviation(_name)
            # Raise errors if method is improperly used
            if mu is not None: raise ValueError('cannot have both name and mu arguments')
            if sigma is not None: raise ValueError('cannot have both name and sigma arguments')

        elif name is None:

            # Raise errors if method is improperly used
            if mu is None: raise ValueError('must either specify name or mu/sigma.')
            if sigma is None: raise ValueError('must either specify name or mu/sigma.')
           
            # Define local variables
            _mu = mu
            _sigma = sigma

        ci_lo = _mu - zscore*_sigma/np.sqrt(_n_samples)
        ci_hi = _mu + zscore*_sigma/np.sqrt(_n_samples)
        return ci_lo,ci_hi

    def get_qoi_ttest(self,name):
        _qoi_target = self.qoi_targets[name]
        [t_stat,p_value] = scipy.stats.ttest_1samp(
                self.df[name],
                _qoi_target)
        return p_value

    def do_analysis(self):
        pass

def test__test_get_confidence_interval():
    o = PyposmatQoiPredictionAnalyzer()
    o.read_configuration(filename=config_fn)
    o.read_datafile(filename=data_fn)

    for q in o.qoi_fitting_names:
        ci_lo, ci_hi = o.get_confidence_interval()
if __name__ == "__main__":

    pyposmat_root = [v for v in os.environ['PYTHONPATH'].strip().split(':') if v.endswith('pypospack')][0]
    data_directory = os.path.join(pyposmat_root,'data','MgO_pareto_data')

    config_fn = os.path.join(pyposmat_root,
            'examples/MgO__buck__add_additional_qoi/data/pyposmat.config.in')
    data_fn = os.path.join(data_directory,"qoiplus_005.out")

    o = PyposmatQoiPredictionAnalyzer()
    o.read_configuration(filename=config_fn)
    o.read_datafile(filename=data_fn)


    s = '{:^20} {:^11} {:^11} {:^11} {:^11} {:^11} {:^11}'.format(
            'qoi_name','mu','sigma','qoi_target','95_CI_low','95_CI_high','p_value')
    print(s) 
    if o.qoi_fitting_names is not None:
        for q in o.qoi_fitting_names:
            _mu  = o.get_mean(q)
            _sigma = o.get_standard_deviation(q)
            _q_target = o.qoi_fitting_targets[q]
            _ci_lo, _ci_hi = o.get_confidence_interval(name=q,zscore=1.96)
            _p_value = o.get_qoi_ttest(q)
            _p_value_truth = _p_value < 0.05
            s = '{:20} {:+10.4e} {:+10.4e} {:+10.4e} {:+10.4e} {:+10.4e} {:+10.4f} {}'.format(
                    q,_mu,_sigma,_q_target,_ci_lo,_ci_hi,_p_value, _p_value_truth)
            print(s)
    else:
        print("No fitting qois")
    if o.qoi_testing_names is not None:
        for q in o.qoi_testing_names:
            _mu  = o.get_mean(q)
            _sigma = o.get_standard_deviation(q)
            _q_target = o.qoi_testing_targets[q]
            _ci_lo, _ci_hi = o.get_confidence_interval(name=q,zscore=1.96)
            _p_value = o.get_qoi_ttest(q)
            _p_value_truth = _p_value < 0.05
            s = '{:20} {:+10.4e} {:+10.4e} {:+10.4e} {:+10.4e} {:+10.4e} {:+10.4f} {}'.format(
                    q,_mu,_sigma,_q_target,_ci_lo,_ci_hi,_p_value, _p_value_truth)
            print(s)
    else:
        print("No testing qois")


