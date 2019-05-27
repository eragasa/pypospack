
class TaskDictionary(dict):
    def __init__(self,*args,**kwargs):
        super(TaskDictionary,self).__init__(*args,**kwargs)
        self.itemlist = super(TaskDictionary,self).keys()

    def __setitem__(self,key,value):
        self.itemlist.append(key)
        super(TaskDictionary,self).__setitem__(key,value)

    def __iter__(self):
        return iter(self.itemlist)

    def __keys__(self):
        return self.itemlist

    def __values__(self):
        return [self[key] for key in self]

    def __itervalues__(self):
        return (self[key]) for key in self)

