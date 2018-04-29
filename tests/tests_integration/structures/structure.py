import json
from collections import OrderedDict

class Structure(object):
    def __init__(self, symbols,stoichiometry,sg):
        self.symbols = list(symbols)
        self.stoichiometriy = OrderedDict(stoichiometry)
        self.sg = sg

    def from_json(json_serial_obj):
        self.__dict__ = json.loads(json_serial_obj)
        return structure

if __name__ == "__main__":
    json_serial_str
