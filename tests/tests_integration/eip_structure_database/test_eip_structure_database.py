import pypospack.potfit as potfit

if __name__ == '__main__':
    structure_db = potfit.StructureDatabase()
    structure_db.directory = "structure_db"
    structure_db.add_structure('MgO_NaCl','MgO_NaCl_unit.vasp','vasp')
    structure_db.add_structure('MgO_NaCl_fr_a','MgO_NaCl_333_fr_a.vasp','vasp')
    structure_db.add_structure('MgO_NaCl_fr_c','MgO_NaCl_333_fr_c.vasp','vasp')
    structure_db.add_structure('MgO_NaCl_sch','MgO_NaCl_333_sch.vasp','vasp')
    structure_db.add_structure('MgO_NaCl_001_s','MgO_NaCl_001_s.vasp','vasp')

    structure_db.write(fname='pypospack.structure.yaml')

    copy_structure_db = potfit.StructureDatabase()
    copy_structure_db.read('pypospack.structure.yaml')
    print(copy_structure_db.directory)
    for k,v in copy_structure_db.structures.items():
        print(k,v)
