import os
from collections import OrderedDict

from pypospack.pyposmat.engines import PyposmatEngine

testing_set = OrderedDict()
testing_set['base_directory'] = 'rank_test'
testing_set['config_fn'] = os.path.join('data','pyposmat.config.in')
testing_set['parameters'] = OrderedDict([
    ('SiSiSi_epsilon', 2.412326433321953), 
    ('SiSiSi_sigma', 2.047323691368865), 
    ('SiSiSi_a', 1.927303534824485), 
    ('SiSiSi_lambda', 12.109230909182934), 
    ('SiSiSi_gamma', 2.242334213331301), 
    ('SiSiSi_costheta0', -0.3333333333333333), 
    ('SiSiSi_A', 19.19755654803364), 
    ('SiSiSi_B', 0.16452955359488666), 
    ('SiSiSi_p', 4.313046208558489), 
    ('SiSiSi_q', 0.7433313641477534), 
    ('SiSiSi_tol', 0.0)
])

def cleanup(testing_set):
    import shutil

    if os.path.exists(testing_set['base_directory']):
        if os.path.isdir(testing_set['base_directory']):
            shutil.rmtree(testing_set['base_directory'])
        elif os.path.isfile(testing_set['base_directory']):
            os.remove(testing_set['base_directory'])
        else:
            raise ValueError('path:{}'.format(testing_set['base_directory']))

if __name__ == "__main__":
    cleanup(testing_set)
    o = PyposmatEngine(
            filename_in = testing_set['config_fn'],
            base_directory = testing_set['base_directory'])
    o.configure()
    o.evaluate_parameter_set(parameters=testing_set['parameters'])
