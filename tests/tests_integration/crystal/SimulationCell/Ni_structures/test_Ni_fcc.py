import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import ase.build
import ase.lattice.cubic
import numpy as np
from numpy.linalg import inv

"""

References:
    for the use of the ASE package:

Todo:
    (1)  Need a more general approach for creating surfaces, the ASE toolkit
         is somewhat lacking in generality unless you already know the
         directions of cell you want.
         https://wiki.fysik.dtu.dk/ase/ase/lattice.html

"""
def get_atomic_radius(symbol):
    if symbol == 'Ni':
        return 1.35

class FaceCenteredCubic(crystal.SimulationCell):

    fcc_111_ortho_vectors = [
            [1,-1,0],
            [1,1,-2],
            [1,1,1]
        ]
    def __init__(self,symbols,structure_type = 'cubic'):
        self.supported_structure_types = ['cubic','orthorhombic','primitive','111','110','100']

        atomic_radius = get_atomic_radius(symbols)
        lattice_parameter = 4*atomic_radius/(2**0.5)
        if structure_type == 'cubic':
            atoms =nsl ase.build.bulk(symbols,'fcc', cubic=True)
            crystal.SimulationCell.__init__(self,atoms)
            self.normalize_h_matrix()
        elif structure_type == 'primitive':
            atoms = ase.build.bulk(symbols,'fcc')
            crystal.SimulationCell.__init__(self,atoms)
        elif structure_type == 'orthorhombic':
            atoms = ase.build.bulk(symbols,'fcc', orthorhombic=True)
            crystal.SimulationCell.__init__(self,atoms)
        elif structure_type == '111':
            # Reference for this approach is
            # https://wiki.fysik.dtu.dk/ase/ase/lattice.html
            x_axis = [1,-1,0]
            y_axis = [1,1,-2]
            z_axis = [1,1,1]
            
            # test for orthogonal vectors
            assert np.dot(np.array(x_axis),np.array(y_axis)) == 0
            assert np.dot(np.array(x_axis),np.array(z_axis)) == 0
            assert np.dot(np.array(y_axis),np.array(z_axis)) == 0

            atoms = ase.lattice.cubic.FaceCenteredCubic(
                        directions=[x_axis,y_axis,z_axis],
                        size=(1,1,1),
                        symbol=symbols,
                        pbc=(1,1,1))
            crystal.SimulationCell.__init__(self,atoms)
        elif structure_type == '110':
            x_axis = [-1,1,0]
            y_axis = [0,0,1]
            z_axis = [1,1,0]
            # test for orthogonal vectors
            assert np.dot(np.array(x_axis),np.array(y_axis)) == 0
            assert np.dot(np.array(x_axis),np.array(z_axis)) == 0
            assert np.dot(np.array(y_axis),np.array(z_axis)) == 0
            atoms = ase.lattice.cubic.FaceCenteredCubic(
                        directions=[x_axis,y_axis,z_axis],
                        size=(1,1,1),
                        symbol=symbols,
                        pbc=(1,1,1))
            crystal.SimulationCell.__init__(self,atoms)
        elif structure_type in self.supported_structure_types:
            raise NotImplementedError('structure_type not implemented')
        else:
            raise ValueError('structure_type not supported')
        self.structure_type = structure_type
        
    def get_interstitial_sites(self):
        if self.structure_type == 'cubic':
            self.sites = []
            self.sites.append([np.array([1/4,1/4,1/4]),'tetrahedral'])
            self.sites.append([np.array([1/4,1/4,3/4]),'tetrahedral'])
            self.sites.append([np.array([1/4,3/4,1/4]),'tetrahedral'])
            self.sites.append([np.array([1/4,3/4,3/4]),'tetrahedral'])
            self.sites.append([np.array([3/4,1/4,1/4]),'tetrahedral'])
            self.sites.append([np.array([3/4,1/4,3/4]),'tetrahedral'])
            self.sites.append([np.array([3/4,3/4,1/4]),'tetrahedral'])
            self.sites.append([np.array([1/2,1/2,1/2]),'octahederal'])
            self.sites.append([np.array([0,0,1/2]),'octahederal'])
            self.sites.append([np.array([0,1/2,0]),'octahederal'])
            self.sites.append([np.array([1/2,0,0]),'octahederal'])
            return sites
        else:
            raise NotImplementedError('interstitial_sites not implemented for {}'.format(self.structure_type))

    def get_kpath(self):
        if self.structure_type == 'primitive':
            # W. Setyawan and S. Curtarolo. Comp. Mat. Sci., 49.2 (2010): 299-312.
            self.kpath_points={
                    'G':[0,0,0],
                    'K':[3/8,3/8,3/8],
                    'L':[1/2,1/2,1/2],
                    'U':[5/8,1/4,5/8],
                    'W':[1/2,0,1/2]
                    }
            self.kpath(['G','X','W','K','L','U','W','L','K'],['U','X'])

class Bcc(crystal.SimulationCell):

    def __init__(self,symbols,structure_type='cubic'):

        self.supported_structure_types = ['cubic','orthorhombic','primitive']
        if structure_type == 'cubic':
            try:
                atoms = ase.build.bulk(symbols,'bcc', cubic=True)
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atomic_radius = get_atomic_radius(symbols)
                    lattice_parameter = 4*atomic_radius/(3**0.5)
                    atoms = ase.build.bulk(symbols,'bcc',\
                            a=lattice_parameter,cubic=True)
                else:
                    raise
            except:
                raise
        elif structure_type == 'primitive':
            try:
                atoms = ase.build.bulk(symbols,'bcc')
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atomic_radius = get_atomic_radius(symbols)
                    lattice_parameter = 4*atomic_radius/(3**0.5)
                    atoms = ase.build.bulk(
                                symbols,
                                'bcc',
                                a=lattice_parameter)
                else:
                    raise
            except:
                raise
        elif structure_type == 'orthorhombic':
            try:
                atoms = ase.build.bulk(symbols,'bcc',orthorhombic=True)
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atomic_radius = get_atomic_radius(symbols)
                    lattice_parameter = 4*atomic_radius/(3**0.5)
                    atoms = ase.build.bulk(\
                                symbols,
                                'bcc',
                                a=lattice_parameter,
                                orthorhombic=True)
                else:
                    raise
            except:
                raise
        elif structure_type in self.supported_structure_types:
            raise NotImplementedError('structure_type not implemented')
        else:
            raise ValueError('structure_type not supported')
        self.structure_type = structure_type
        
        crystal.SimulationCell.__init__(self,atoms)
        self.normalize_h_matrix()

class Sc(crystal.SimulationCell):
    def __init__(self,symbols,structure_type='cubic'):
        atomic_radius = get_atomic_radius(symbols)
        lattice_parameter = 2*atomic_radius

        self.supported_structure_types = ['cubic']
        if structure_type == 'cubic':
            try:
                atoms = ase.build.bulk(symbols,'sc',cubic=True)
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atoms = ase.build.bulk(symbols,'sc',
                                a=lattice_parameter,cubic=True)
                else:raise
            except:raise
        elif structure_type in self.supported_structure_types:
            raise NotImplementedError('structure_type not implemented')
        else:
            raise ValueError('structure_type not supported')
        self.structure_type = structure_type

        crystal.SimulationCell.__init__(self,atoms)
        self.normalize_h_matrix()

class Hcp(crystal.SimulationCell):

    def __init__(self,symbols,structure_type='cubic'):
        atomic_radius = get_atomic_radius(symbols)
        lattice_parameter = 2*atomic_radius

        self.supported_structure_types = ['cubic','orthorhombic','primitive']
        if structure_type == 'cubic':
            try:
                atoms = ase.build.bulk(symbols,'hcp', cubic=True)
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atoms = ase.build.bulk(symbols,'hcp',
                                a=lattice_parameter,cubic=True)
                else:raise
            except:raise
        elif structure_type == 'primitive':
            try:
                atoms = ase.build.bulk(symbols,'hcp')
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atoms = ase.build.bulk(symbols,'hcp',
                                a=lattice_parameter)
                else:raise
            except:raise
        elif structure_type == 'orthorhombic':
            try:
                atoms = ase.build.bulk(symbols,'hcp', orthorhombic=True)
            except ValueError as e:
                s = str(e)
                if s == 'You need to specify the lattice constant.':
                    atoms = ase.build.bulk(symbols,'hcp',
                                a=lattice_parameter,orthorhombic=True)
                else:raise
            except:raise
        elif structure_type in self.supported_structure_types:
            raise NotImplementedError('structure_type not implemented')
        else:
            raise ValueError('structure_type not supported')
        self.structure_type = structure_type
        
        crystal.SimulationCell.__init__(self,atoms)
        self.normalize_h_matrix()


class SurfaceSlab(object):
    pass


class Fcc100Slab(SurfaceSlab):

    def __init__(self,slab_thickness):
        pass

if __name__ == "__main__":

    print('making Ni fcc, cubic cell')
    fcc = Fcc('Ni',structure_type='cubic')
    fcc_vasp = vasp.Poscar(fcc)
    print('fcc,cubic,lattice:',fcc.a0)
    print('fcc,cubic,a1',fcc.H[0,:])
    print('fcc,cubic,a2',fcc.H[1,:])
    print('fcc,cubic,a3',fcc.H[2,:])
    fcc_vasp.write('Ni_fcc_cubic.vasp')

    print('making Ni fcc, primitive cell')
    fcc = Fcc('Ni',structure_type='primitive')
    fcc_vasp = vasp.Poscar(fcc)
    print('fcc,cubic,lattice:',fcc.a0)
    print('fcc,cubic,a1',fcc.H[0,:])
    print('fcc,cubic,a2',fcc.H[1,:])
    print('fcc,cubic,a3',fcc.H[2,:])
    fcc_vasp.write('Ni_fcc_prim.vasp')

    print('making Ni fcc, orthorhombic cell')
    fcc = Fcc('Ni',structure_type='orthorhombic')
    fcc_vasp = vasp.Poscar(fcc)
    fcc_vasp.write('Ni_fcc_ortho.vasp')

    print('making Ni bcc, cubic cell')
    bcc = Bcc('Ni',structure_type='cubic')
    bcc_vasp = vasp.Poscar(bcc)
    bcc_vasp.write('Ni_bcc_cubic.vasp')

    print('making Ni bcc, primitive cell')
    bcc = Bcc('Ni',structure_type='primitive')
    bcc_vasp = vasp.Poscar(bcc)
    bcc_vasp.write('Ni_bcc_prim.vasp')

    print('making Ni bcc, orthorhombic cell')
    bcc = Bcc('Ni',structure_type='orthorhombic')
    bcc_vasp = vasp.Poscar(bcc)
    bcc_vasp.write('Ni_bcc_ortho.vasp')

    print('making Ni hcp, cubic cell')
    hcp = Hcp('Ni',structure_type='cubic')
    hcp_vasp = vasp.Poscar(hcp)
    hcp_vasp.write('Ni_hcp_cubic.vasp')

    print('making Ni hcp, primitive cell')
    hcp = Hcp('Ni',structure_type='primitive')
    hcp_vasp = vasp.Poscar(hcp)
    hcp_vasp.write('Ni_hcp_prim.vasp')

    print('making Ni hcp, orthorhombic cell')
    hcp = Hcp('Ni',structure_type='orthorhombic')
    hcp_vasp = vasp.Poscar(hcp)
    hcp_vasp.write('Ni_hcp_ortho.vasp')
    
    print('making Ni sc, cubic_cell')
    sc = Sc('Ni',structure_type='cubic')
    sc_vasp = vasp.Poscar(sc)
    sc_vasp.write('Ni_sc_cubic.vasp')

    print('making Ni fcc 111 orientation, ortho cell')
    fcc_111 = Fcc('Ni',structure_type='111')
    fcc_111 = vasp.Poscar(fcc_111)
    fcc_111.write('Ni_fcc111_ortho.vasp')

    print('making Ni fcc 110 orientation, ortho cell')
    fcc_110 = Fcc('Ni',structure_type='110')
    fcc_110 = vasp.Poscar(fcc_110)
    fcc_110.write('Ni_fcc110_ortho.vasp')

    exit()
    print('b1:',fcc.b1)
    print('b2:',fcc.b2)
    print('b3:',fcc.b3)

    H = fcc.H
    G = np.zeros(shape=(3,3))
    G[0,:] = fcc.b1
    G[1,:] = fcc.b2
    G[2,:] = fcc.b3

    print(G)
    print(2*np.pi*inv(H.T))

