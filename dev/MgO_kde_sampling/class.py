class myclass:

    self.myproperty = None
    self._myproperty = something

    def myfunction(self):
        if self.myproperty is None:
            self.myproperty = self.calculate_my_properties()

    @property
    def parameters(self):
        return self.qoi.configuration.parameters

    def property(self,myproperty:
        self._myproperty=myproperty
        self.do_something_with_my_property(myproperty)
