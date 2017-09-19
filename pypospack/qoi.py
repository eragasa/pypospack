import copy
import numpy as np

class QoiManager:
    """ Manager of Quantities of Interest 
    
    This class manages quantities of interest the simulation of multiple 
    quantities of interest will often require the results from the same
    simulation.  The purpose of this class is to identify the values required
    from different simulations, and then identify which simulations need to
    be done first.

    Args:

    Attributes:
        qoi_names (list): list of qoi_names
        qoi_targets (dict): key is the qoi_name, value is the reference value
        qoi_errors (dict): key is the qoi_name, value is the error
        qoi_dict (dict): 
    """ 
    def __init__(self):

        self._qoi_names = []
        self._qoi_targets = {}
        self._qoi_values = {}
        self._qoi_errors = {}
        self._qoi_dict = {}

        #TODO: this should be switched to a sorted list
        self._req_structure_names = None
        self._req_structures = {}

        #TODO: this should be switched to a sorted list
        self._req_variable_names = None
        self._req_variables = {}

    @property
    def qoi_names(self):
        return self._qoi_names

    @qoi_names.setter
    def qoi_names(self, qoi_names):
        assert type(qoi_names), list
        self._qoi_names = qoi_names 

    @property
    def qoi_values(self):
        return self._qoi_values

    @property
    def qoi_errors(self):
        return self._qoi_errors
    @property
    def qoi_definitions(self):
        return self._qoi_definitions

    @qoi_definitions.setter
    def qoi_definitions(self, qoi_defs):
        self._qoi_definitions = qoi_defs

    @property
    def qoi_targets(self):
        return self._qoi_targets

    @qoi_targets.setter
    def qoi_targets(self,qoi_targets):
        self._qoi_targets = qoi_targets

    @property
    def required_structure_names(self):
        if self._req_structure_names is None:
            self._set_required_structure_names()
 
    @property
    def required_variable_names(self):
        if self._req_variable_names is None:
            self._set_required_variable_names()
        return self._req_variable_names

        return self._req_structure_names

    @property
    def required_variables(self):
        return self._req_variables

    @required_variables.setter
    def required_variables(self, req_var_dict):
        self._req_variables = req_var_dict
 
    def process_qoi_definitions(self):
        for k in self._qoi_definitions:
            qoi = self._qoi_definitions[k]['variable']
            s   = self._qoi_definitions[k]['structure']
            if qoi in ['a0','a1','a2','a3','alpha','beta','gamma']:
                self.add_qoi(k, CrystalStructureGeometry(k,qoi,s))
            elif qoi in ['c11','c12','c44']:
                self.add_qoi(k, ElasticTensorComponent(k,qoi,s))
            elif qoi == 'bulk_modulus':
                self.add_qoi(k, BulkModulus(k,s))
            elif qoi == 'shear_modulus':
                self.add_qoi(k, ShearModulus(k,s))
            elif qoi == 'defect_energy':
                self.add_qoi(k, DefectFormationEnergy(k,s))
            elif qoi == 'surface_energy':
                self.add_qoi(k, SurfaceEnergy(k,s))
            elif qoi == 'stacking_fault':
                self.add_qoi(k, StackingFaultEnergy(k,s))
            else:
                err_msg = "can't add qoi [{}] unknown qoi_type [{}]"
                err_msg = err_msg.format(k,qoi)
                raise ValueError(err_msg)
            self._qoi_dict[k].reference_value = self._qoi_definitions[k]['target']

    def add_qoi(self,qoi_name, qoi_obj):
        self._qoi_names.append(qoi_name)
        self._qoi_dict[qoi_name] = copy.deepcopy(qoi_obj)

    def calculate_qois(self,var_dict):
        self._req_variables = var_dict
        self._qoi_values = {}
        self._qoi_errors = {}
        for k in self._qoi_dict.keys():
            for v_name in self._qoi_dict[k].required_variable_names:
                v_value = self._req_variables[v_name]
                self._qoi_dict[k].set_variable(v_name,v_value)
            self._qoi_dict[k].calculate_qoi()
            self._qoi_values[k] = self._qoi_dict[k].predicted_value
            self._qoi_errors[k] = self._qoi_dict[k].error

    def _set_required_structure_names(self):
        self._req_structure_names = []
        for k,v in self._qoi_dict.items():
            for s in v.required_structure_names:
                if s not in self._req_structure_names:
                    self._req_structure_names.append(s)
        self._req_structure_names.sort()
        return self._req_structure_names 

    def _set_required_variable_names(self):
        self._req_variable_names = []
        for k,v in self._qoi_dict.items():
            for var in v.required_variable_names:
                if var not in self._req_variable_names:
                    self._req_variable_names.append(var)
        self._req_variable_names.sort()
        return self._req_variable_names

class Qoi:
    def __init__(self, qoi_name, qoi_type, req_structures_names):
        self.qoi_name = qoi_name
        self.qoi_type = qoi_type

        self.ref_value = None
        self.predicted_value = None
        self.error = None

        # set required structure names and set up the dictionary
        self.req_structure_names = list(req_structures_names)
        self.req_structures = {}
        for sn  in self._req_structure_names:
            self._req_structures[sn] = None

        # set required variable names and set up the dictionary
        self.req_var_names = [qoi_name]
        self.req_vars = {}
        for var in self._req_var_names:
            self.req_vars[var] = None

    def calculate_qoi(self):
        raise NotImplementedError

    @property
    def reference_value(self):
        return self._ref_value

    @reference_value.setter
    def reference_value(self, value):
        assert type(value), float
        self._ref_value = value

    @property
    def predicted_value(self):
        return self._predicted_value 

    @predicted_value.setter
    def predicted_value(self, qhat):
        self._predicted_value = qhat
    def set_variable(self, v_name, v_value):
        self._req_vars[v_name] = v_value

class CrystalStructureGeometry(Qoi):
    def __init__(self,qoi_name,qoi_type,structures):
        Qoi.__init__(self,qoi_name,qoi_type,structures)

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        self._predicted_value = self._req_vars[self._qoi_name]
        return self._predicted_value

class ElasticTensorComponent(Qoi):
    def __init__(self,qoi_name,qoi_type,structures):
        Qoi.__init__(self,qoi_name,qoi_type,structures)

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        self._predicted_value = self._req_vars[self._qoi_name]
        return self._predicted_value

class CohesiveEnergy(Qoi):
    def __init__(self,qoi_name,structures):
        qoi_type = 'E_coh'
        Qoi.__init__(self,qoi_name,qoi_type,structures)

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        raise NotImplementedError

class BulkModulus(Qoi):
    def __init__(self,qoi_name, structures):
        qoi_type = 'bulk_modulus'
        assert len(structures), 1
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self._req_var_names = []
        self._req_var_names.append("{}.{}".format(structures[0],'c11'))
        self._req_var_names.append("{}.{}".format(structures[0],'c12'))
        self._req_var_names.append("{}.{}".format(structures[0],'c44'))

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        s_name = self._req_structure_names[0]
        c11 = self._req_vars["{}.c11".format(s_name)]
        c12 = self._req_vars["{}.c12".format(s_name)]
        c44 = self._req_vars["{}.c44".format(s_name)]
        self._predicted_value = (c11+2*c12)/3.
        return self._predicted_value


class ShearModulus(Qoi):
    def __init__(self,qoi_name, structures):
        qoi_type = 'shear_modulus'
        assert len(structures), 1
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self._req_var_names = []
        self._req_var_names.append("{}.{}".format(structures[0],'c11'))
        self._req_var_names.append("{}.{}".format(structures[0],'c12'))
        self._req_var_names.append("{}.{}".format(structures[0],'c44'))

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        s_name = self._req_structure_names[0]
        c11 = self._req_vars["{}.c11".format(s_name)]
        c12 = self._req_vars["{}.c12".format(s_name)]
        c44 = self._req_vars["{}.c44".format(s_name)]
        self.predicted_value = (c11-c12)/2.
        return self.predicted_value

class DefectFormationEnergy(Qoi):
    def __init__(self,qoi_name, structures):
        qoi_type = 'defect_energy'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self._req_var_names = []
        self._req_var_names.append("{}.{}".format(structures[0],'E_min_pos'))
        self._req_var_names.append("{}.{}".format(structures[0],'n_atoms'))
        self._req_var_names.append("{}.{}".format(structures[1],'E_min'))
        self._req_var_names.append("{}.{}".format(structures[1],'n_atoms'))

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        s_name_defect = self._req_structure_names[0]
        s_name_bulk   = self._req_structure_names[1]
        e_defect = self._req_vars["{}.E_min_pos".format(s_name_defect)]
        e_bulk   = self._req_vars["{}.E_min".format(s_name_bulk)]
        n_atoms_defect = self._req_vars["{}.n_atoms".format(s_name_defect)]
        n_atoms_bulk   = self._req_vars["{}.n_atoms".format(s_name_bulk)]
        e_f = e_defect - n_atoms_defect/n_atoms_bulk*e_bulk
        self._predicted_value = e_f
        return self._predicted_value

class SurfaceEnergy(Qoi):
    def __init__(self, qoi_name, structures):
        qoi_type = 'surface_energy'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self._req_var_names = []
        self._req_var_names.append("{}.{}".format(structures[0],'E_min_pos'))
        self._req_var_names.append("{}.{}".format(structures[0],'a1_min_pos'))
        self._req_var_names.append("{}.{}".format(structures[0],'a2_min_pos'))
        self._req_var_names.append("{}.{}".format(structures[0],'n_atoms'))
        self._req_var_names.append("{}.{}".format(structures[1],'E_min'))
        self._req_var_names.append("{}.{}".format(structures[1],'n_atoms'))

        self._req_vars = {}
        for var in self._req_var_names:
            self._req_vars[var] = None

    def calculate_qoi(self):
        s_name_slab = self._req_structure_names[0]
        s_name_bulk = self._req_structure_names[1]
        e_slab = self._req_vars["{}.E_min_pos".format(s_name_slab)]
        e_bulk = self._req_vars["{}.E_min".format(s_name_bulk)]
        n_atoms_slab = self._req_vars["{}.n_atoms".format(s_name_slab)]
        n_atoms_bulk = self._req_vars["{}.n_atoms".format(s_name_bulk)]
        a1 = self._req_vars["{}.a1_min_pos".format(s_name_slab)]
        a2 = self._req_vars["{}.a2_min_pos".format(s_name_slab)]
        e_surf = (e_slab - n_atoms_slab/n_atoms_bulk*e_bulk)/(2*a1*a2)
        self._predicted_value = e_surf
        return self._predicted_value

class StackingFaultEnergy(Qoi):

