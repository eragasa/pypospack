import copy, yaml
from collections import OrderedDict
import pypospack.pyposmat as pyposmat

class QoiDatabase(object):
    """ Qoi Database 
    
        Attributes:
            filename(str): file to read/write to yaml file
    """
        
    def __init__(self):
        self.filename = 'pypospack.qoi.yaml'
        self.qois = {}

    def add_qoi(self,
            name,
            qoi_type,
            structures,
            target):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """
        assert isinstance(name,str)
        assert isinstance(qoi_type,str)
        assert any([
            isinstance(structures,str),
            isinstance(structures,list),
            isinstance(structures,dict),
            ])
        assert isinstance(target,float)
        _structures = None

        self.qois[name] = {\
                'qoi_type':qoi_type,
                'structures':copy.deepcopy(structures),
                'target':target}

    def read(self,fname=None): 
        """ read qoi configuration from yaml file

        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname
        try:
            self.qois = yaml.load(open(self.filename))
        except:
            raise

    def write(self,fname=None):
        """ write qoi configuration to yaml file

        Args:
            fname(str): file to write yaml from from.  If no argument is passed
               then use the filename attribute.  If the filename is set, then 
               the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dictionary
        self.qoi_db = copy.deepcopy(self.qois)

        # dump to yaml file
        with open(self.filename,'w') as f:
            yaml.dump(self.qoi_db,f, default_flow_style=False)

if __name__ == '__main__':
    qoi_db = QoiDatabase()
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.a0', 
            qoi_type = 'a0_min_all',
            structures = ['MgO_NaCl'],
            target = 4.246)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c11', 
            qoi_type = 'c11',
            structures = ['MgO_NaCl'],
            target = 277.00)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c12', 
            qoi_type = 'c12',
            structures = ['MgO_NaCl'],
            target = 91.67)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c44', 
            qoi_type = 'c44',
            structures = ['MgO_NaCl'],
            target = 144.01)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.B', 
            qoi_type = 'bulk_modulus',
            structures = ['MgO_NaCl'],
            target = 153.45)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.G', 
            qoi_type = 'shear_modulus',
            structures = ['MgO_NaCl'],
            target = 92.66)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.fr_a', 
            qoi_type = 'defect_energy',
            structures = OrderedDict((
                ['defect','MgO_NaCl_fr_a'],
                ['reservoir','MgO_NaCl']
            )),
            target = 10.978)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.fr_c', 
            qoi_type = 'defect_energy',
            structures = OrderedDict((
                ['defect','MgO_NaCl_fr_c'],
                ['reservoir','MgO_NaCl']
            )),
            target = 8.986)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.sch', 
            qoi_type = 'defect_energy',
            structures = OrderedDict((
                ['defect','MgO_NaCl_sch'],
                ['reservoir','MgO_NaCl']
            )),
            target = 5.067)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.001s', 
            qoi_type = 'surface_energy',
            structures = ['MgO_NaCl_001_s','MgO_NaCl'],
            target = 0.05595)

    qoi_db.write(fname='pypospack.qoi.yaml')

    copy_structure_db = QoiDatabase()
    copy_structure_db.read('pypospack.qoi.yaml')
    for k,v in copy_structure_db.qois.items():
        print(k,v)
