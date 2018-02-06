
    def initialize_kde_sampler(self,fname_in):
        f = open(fname_in,'r')
        lines = f.readlines()
        f.close()

        self._kde_names = lines[0].strip().split(',')
        self._kde_names = [str(v.strip()) for v in self._kde_names]

        self._kde_types = lines[1].strip().split(',')
        self._kde_types = [str(v).strip() for v in self._kde_types]

        datas = []
        for i in range(2,len(lines)):
            line = lines[i]
            line = line.strip().split(',')
            line = [float(v) for v in line]
            datas.append(line)
        datas = np.array(datas)

        # construct free parameter list
        free_param_list = self.get_free_parameter_list()
        self._kde_free_param_indx = []
        for i,v in enumerate(self._kde_names):
            if v in free_param_list:
                self._kde_free_param_indx.append(i)
        # DEBUG
        free_params = datas[:,self._kde_free_param_indx]
        self._kde_kernel = scipy.stats.gaussian_kde(free_params.transpose())
    
    def get_kde_sample(self):
        if self._kde_kernel is None:
            self.initialize_kde_sampler(self._fname_results_in)

        # initialize param dict
        param_dict = {}
        for pn in self._param_names:
            param_dict[pn] = None

        # construct free parameter list
        free_param_list = self.get_free_parameter_list()

        is_good = False
        while not is_good:
            # sample free parameters from kde
            free_params = self._kde_kernel.resample(size=1)
            for i,pn in enumerate(free_param_list):
                param_dict[pn] = free_params[i,0]

            # static variables
            for pn in self._param_names:
                if self._param_info[pn]['type'] == 'static':
                    param_dict[pn] = self._param_info[pn]['info'][0]

            # constrained variables
            for pn in self._param_names:
                if self._param_info[pn]['type'] == 'equals':
                    info = self._param_info[pn]['info'][0]
                    for p in self._param_names:
                        if p in info:
                            info = info.replace(p,"{}".format(param_dict[p]))
                    param_dict[pn] = eval(info)

            # check parameter constraints
            if param_dict['MgMg_rho'] < 0.:
                is_good = False
            elif param_dict['MgO_rho'] < 0.:
                is_good = False
            elif param_dict['OO_rho'] < 0.:
                is_good = False
            else:
                is_good =True

        return param_dict
