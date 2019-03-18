import json
from collections import OrderedDict


class PyposmatReferenceDataFile(object):
    """ Used to interface with data from publications for validation purposes.
    
    Args:
        filename(str, optional): the filename of the data file to read (must be JSON).
    """

    def __init__(self, filename=None):
        self.filename = filename
        self._bibtex = None
        self._qois = None

    @property
    def bibtex(self):
        if self._bibtex is None:
            raise RuntimeError("bibtex value not set or not found in json")
        else:
            return self._bibtex

    @property
    def qois(self):
        if self._qois is None:
            raise RuntimeError("qois value not set")
        else:
            return self._qois

    def read(self, filename=None):
        # filename uses self.filename if not set
        if filename is None:
            if self.filename is None:
                raise ValueError("must provide a filename to read")
            else:
                filename = self.filename
        # filename must refer to a JSON file
        if not filename.endswith(".json"):
            raise ValueError("filename must be a json file")
        # read JSON into a dict
        with open(filename) as f:
            raw_json = json.load(f)
        # set the private attributes which the properties pull from
        self._bibtex = raw_json.pop("bibtex", None)
        qois =  raw_json.pop("qois", None)
        self._qois = OrderedDict(qois)