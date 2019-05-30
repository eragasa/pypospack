import pypospack.potfit as potfit
import pypospack.qoi as qoi

if __name__ == "__main__":
    qoi_structures = ['NaCl_defect','NaCl_unit']
    qoi_name = 'NaCl.dfe'

    obj_qoi = qoi.DefectFormationEnergy(
            qoi_name = qoi_name,
            structures = qoi_structures)

    print(obj_qoi.get_required_simulations())