import pypospack.potfit as potfit
import os
import numpy as np

if __name__ == "__main__":
    data_dir = 'data'
    filename = 'culled_009.out'

    # configure the fitting engine
    eip_config = potfit.EipFittingEngine(\
            fname_config_potential = 'pypospack.buckingham.yaml',
            fname_config_qoi = 'pypospack.qoi.yaml',
            fname_config_structures = 'pypospack.structure.yaml')
   
    fname = os.path.join(data_dir,filename)

    # persistent variables
    names = []
    name_types = []
    data = []
    # read in the line
    lines = None
    try:
        with open(fname,'r') as f:
            lines = f.readlines()
    except:
        raise
    
    for i,line in enumerate(lines):
        line = line.strip()
        if i == 0:
            names = [v.strip() for v in line.split(',')]
        elif i == 1:
            name_types = [v.strip() for v in line.split(',')]
        else:
            data_line = [float(v.strip()) for v in line.split(',')]
            data.append(data)

    assert len(names) == len(name_types)

    # recast list of lists into numpy array
    data = np.array(data)

    param_names = []
    param_key_index = []
    for i,v in enumerate(names):
        if name_types == 'param':
            param_names.append(v)
            param_key_index.append(i)

    param_data = data[:,param_key_index]
    # this is broken for seaton to fix, and talk about later
    param_df = pd.DataFrame(param_data,param_names)

    print('names:',names)
    print('name_types:',name_types)
    print('data:\n\tnumber of entries:',len(data))
