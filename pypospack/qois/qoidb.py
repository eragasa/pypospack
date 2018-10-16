import copy
from collections import OrderedDict
from pypospack.io.filesystem import OrderedDictYAMLLoader


class QoiDatabase(object):
    """ Qoi Database 
    
        Attributes:
            filename(str): file to read/write to yaml file
    """
        
    def __init__(self,qoidb_OrderedDict=None):
        assert any([
            isinstance(qoidb_OrderedDict,OrderedDict),
            type(qoidb_OrderedDict) in [type(None)],
            ])

        self.filename = 'pypospack.qoi.yaml'
        self.qois = None
        self.qoi_names = None

        if qoidb_OrderedDict is not None:
            self.__init_from_OrderedDict(qoidb=qoidb_OrderedDict)

    def __init_from_OrderedDict(self,qoidb):
        for k_qoi,v_qoi in qoidb.items():
            _qoi_name = k_qoi
            _qoi_type = v_qoi['qoi_type']
            _structures = v_qoi['structures']
            _target = v_qoi['target']

            try:
                _qoi_options = v_qoi['qoi_options']
            except KeyError as e:
                _qoi_options = None
            self.add_qoi(
                    qoi_name=_qoi_name,
                    qoi_type=_qoi_type,
                    structures=_structures,
                    target=_target,
                    qoi_options=_qoi_options)
        
    def add_qoi(self,
            qoi_name,
            qoi_type,
            structures,
            target,
            qoi_options=None):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """
        assert isinstance(qoi_name,str)
        assert isinstance(qoi_type,str)
        assert any([
            isinstance(structures,str),
            isinstance(structures,list),
            isinstance(structures,dict),
            ])
        assert isinstance(target,float)
        
        _structures = None
        if isinstance(structures,list):
            _structures = OrderedDict()
            _structures['ideal'] = structures[0]
        elif isinstance(structures,OrderedDict):
            _structures = copy.deepcopy(structures)
        else:
            raise ValueError

        #<--------- initialize internal atributes if not already set
        if self.qois is None: self.qois = OrderedDict()
        if self.qoi_names is None: self.qoi_names = []

        #<--------- create the dictionary entry for this qoi
        self.qois[qoi_name] = OrderedDict()
        self.qois[qoi_name]['qoi_type'] = qoi_type
        self.qois[qoi_name]['structures'] = copy.deepcopy(_structures)
        self.qois[qoi_name]['target'] = target

        if qoi_options is not None:
            self.qois[qoi_name]['qoi_options']= qoi_options 
        #<--------- let's add the value for qoi names
        self.qoi_names.append(qoi_name)

    
    def read(self,filename=None): 
        """ read qoi configuration from yaml file
        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """
        assert isinstance(filename,str)

        # set the attribute if not none
        if filename is not None:
            self.filename = filename
        
        try:
            with open(self.filename,'r') as f:
                self.qois = yaml.load(f, OrderedDictYAMLLoader)
        except:
            raise

        # <------------------ 
        self.qoi_names = [k for k in self.qois]

    def write(self,filename=None):
        """ write qoi configuration to yaml file

        Args:
            fname(str): file to write yaml from from.  If no argument is passed
               then use the filename attribute.  If the filename is set, then 
               the filename attribute is also set.
        """

        # set the attribute if not none
        if filename is not None:
            self.filename = filename

        # marshall attributes into a dictionary
        _qoidb = copy.deepcopy(self.qois)

        # dump to yaml file
        with open(self.filename,'w') as f:
            yaml.dump(_qoidb,f, default_flow_style=False)

    def to_string(self):
        s = 80*'-'+"\n"
        s += '{:^80}\n'.format('QUANTITIES OF INTEREST')
        s += 80*'-'+"\n"
        for k,v in self.qois.items():
            s += '{}\n'.format(k)
            
        
        return s
#------------------------------------------------------------------------------


