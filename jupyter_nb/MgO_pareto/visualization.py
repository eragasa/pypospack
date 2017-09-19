import pypospack.potfit as potfit
import os

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

    print('names:',names)
    print('name_types:',name_types)
    print('data:\n\tnumber of entries:',len(data))
