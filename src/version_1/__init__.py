
# notes here for

# input file

class Atat():
    def __init__(self): pass
       self.lattice_file = None
       self.structure_file = None

    def create_input_files(self):
        pass

# SQS
# References:
# S.-H. Wei et al., Phys. Rev. B 42, 9622 (1990)
# A. Zunger et al., Phys. Rev. Lett. 65, 353 (1990)

class AtatLatticeFile():
    def __init__(self):
        self.path = None

    def get_coordinate_system(self,cell):
        if isinstance(simulation_cell):
            a1 = simulation_cell.a1
            a2 = simulation_cell.a2
            a3 = simulation_cell.a3

            self.a = np.sqrt(np.dot(a1,a1))
            self.b = np.sqrt(np.dot(a2,a2))
            self.c = np.sqrt(np.dot(a3,a3))

            self.alpha = 'todo'
            self.beta = 'todo'
            self.gamma = 'todo'
# https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/tutorial/

class AtatStructureFile():
    def __init__(self):
        self.path = None

    def read(self,path):
        with open(path,'r') as f:
            lines =  f.readlines()

    def write(self,path):
        raise NotImplementedError

class AtatClustersFile():
    def __init__(self):
        self.path = 'clusters.out'

    def read(self,path):
        with open(path,'r') as f:
            lines =  f.readlines()

    def write(self,path):
        raise NotImplementedError

class AtatEciFile():
    def __init_(self):
        self.path = 'eci.out'

    def read(self,path):
        with open(path,'r') as f:
            lines =  f.readlines()

    def write(self,path):
        raise NotImplementedError

class AtatFvibFile():
    def __init__(self):
        self.path = "fvib.eci"

    def read(self,path):
        with open(path,'r') as f:
            lines =  f.readlines()

    def write(self,path):
        raise NotImplementedError

class AtatTeciFile():
    def __init__(self):
        self.path = "teci.out"

    def read(self,path):
        with open(path,'r') as f:
            lines =  f.readlines()

    def write(self,path):
        raise NotImplementedError
