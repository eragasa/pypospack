import ase.build
import pypospack.crystallography as xtal

def build_bulk(symbols,sg,a,c=None,cubic=None,orthorhombic=None):
    if sg in ['fcc','bcc']:
        try:
            ase_cell = ase.build.bulk(symbols[0],sg,a=a,cubic,orthorhombic)
        except:
            raise
        else:
            pass
        finally:
            return ase
simulation_cells = {}
simulation_cells['Ni_fcc_unit']=
