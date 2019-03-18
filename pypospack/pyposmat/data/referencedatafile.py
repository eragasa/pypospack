import json


class PyposmatReferenceDataFile(object):

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
        if filename is None:
            if self.filename is None:
                raise ValueError("must provide a filename to read")
            else:
                filename = self.filename

        if not filename.endswith(".json"):
            raise ValueError("filename must be a json file")
        
        with open(filename) as f:
            raw_json = json.load(f)

        self._bibtex = raw_json.pop("bibtex", None)
        self._qois =  raw_json
