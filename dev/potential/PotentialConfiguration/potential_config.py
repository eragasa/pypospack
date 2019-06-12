from collections import OrderedDict

class PotentialConfiguration(object):

    def __init__(self):

    def read(self,filename):
        self.filename = filename
        self.configuration = None
        with open(filename,'r') as f:
            self.configuration = yaml.load(f,OrderedDictYamlLoader)

def validate_potential_configuration(d_config):
    validate_pair_potential(d_config['pair'])
    validate_density_function(d_config['density'])

def validate_potential_configuration(d_config):
    if d_config['formalism'] == 'eam':
        validate_eam_potential_configuration(d_config)
if __name__ == "__main__":
    potential = OrderedDict()
    potential['symbols']=['Ni']
    potential['formalism'] = 'eam'
    potential['pair'] = OrderedDict()
    potential['pair']['NiNi']['formalism'] = None
    potential['pair']['NiNi']['param'] = None
    potential['density']['Ni'] = OrderedDict()
    potential['density']['Ni']['formalism'] = None
    potential['density']['Ni']['param'] = None
    potential['embedding']['Ni']['formalism'] = None
    potential['embedding']['Ni']['param']
