
class TersoffPotential as object:
    # Refernce:
    # http://lammps.sandia.gov/doc/pair_tersoff.html
    def __init__(self):
        self.elements = None
        self.param_dict = None

    def _create_parameter_dictionary(self, elements = None):
        if elements is not None:
            self.elements = list(elements)
        self.param_dict = []
        for el1 in self.elements:
            for el2 in self.elements:
                for el3 in self.elements:
                    self._create_param_dict_entry(e1,e2,e3)

    def _create_parameter_dict_entry(self, el1, el2, el3):
        name = "{}.{}.{}".format(el1,el2,el3)
        self.param_dict[name] = {}
        self.param_dict[name]['element1'] = el1 # center atom
        self.param_dict[name]['element2'] = el2 # atom bonded to el1
        self.param_dict[name]['element3'] = el3 # atom influencing 1-2 bond in bond order sense
        self.param_dict[name]['m'] = None
        self.param_dict[name]['gamma'] = None
        self.param_dict[name]['lambda3'] = None
        self.param_dict[name]['c'] = None
        self.param_dict[name]['d'] = None
        self.param_dict[name]['costheta0'] = None
        self.param_dict[name]['n'] = None
        self.param_dict[name]['beta'] = None
        self.param_dict[name]['lambda2'] = None
        self.param_dict[name]['B'] = None
        self.param_dict[name]['R'] = None
        self.param_dict[name]['D'] = None
        self.param_dict[name]['lambda1'] = None
        self.param_dict[name]['A'] = None

    def _check_parameters(self):
        # -1 < costheta0 < 1
