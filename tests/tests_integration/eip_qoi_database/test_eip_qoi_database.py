import copy, yaml
import pypospack.potfit as potfit

class QoiDatabase(object):
    """ Qoi Database 
    
        Attributes:
            filename(str): file to read/write to yaml file

    """
        
    def __init__(self):
        self.filename = 'pypospack.qoi.yaml'
        self.qois = {}

    def add_qoi(self,name,qoi,structures,target):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """

        # make a copy of structures
        _structures = None
        if isinstance(structures,str):
            _structures = [structures]
        else:
            _structures = list(structures)

        self.qois[name] = {\
                'qoi':qoi,
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
            name = 'MgO_NaCl.a0', qoi = 'a0',
            structures = ['MgO_NaCl'],
            target = 4.246)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c11', qoi = 'c11',
            structures = ['MgO_NaCl'],
            target = 277.00)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c12', qoi = 'c12',
            structures = ['MgO_NaCl'],
            target = 91.67)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.c44', qoi = 'c44',
            structures = ['MgO_NaCl'],
            target = 144.01)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.B', qoi = 'bulk_modulus',
            structures = ['MgO_NaCl'],
            target = 153.45)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.G', qoi = 'shear_modulus',
            structures = ['MgO_NaCl'],
            target = 92.66)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.fr_a', qoi = 'defect_energy',
            structures = ['MgO_NaCl_fr_a','MgO_NaCl'],
            target = 10.978)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.fr_c', qoi = 'defect_energy',
            structures = ['MgO_NaCl_fr_c', 'MgO_NaCl'],
            target = 8.986)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.sch', qoi = 'defect_energy',
            structures = ['MgO_NaCl_sch','MgO_NaCl'],
            target = 5.067)
    qoi_db.add_qoi(\
            name = 'MgO_NaCl.001s', qoi = 'surface_energy',
            structures = ['MgO_NaCl_001_s','MgO_NaCl'],
            target = 0.05595)

    qoi_db.write(fname='pypospack.qoi.yaml')

    copy_structure_db = QoiDatabase()
    copy_structure_db.read('pypospack.qoi.yaml')
    for k,v in copy_structure_db.qois.items():
        print(k,v)
