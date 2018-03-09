import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def str__qoi_constraints(qoi_constraints):
    s = []

    for k,v in qoi_constraints.items():
        if k == 'qoi_constraints':
            s.append(k)
            for qoic_n,qoic_v in v.items():
                s.append("{:5} {:>20} {:^3} {:<10}".format("",qoic_n,qoic_v[0],qoic_v[1]))
        else:
            s.append(",".join([str(k),str(v)]))
    return "\n".join(s)
def print__qoi_constraints(qoi_constraints):
    s = str__qoi_constraints(qoi_constraints)
    print(s)
class Dev__PyposmatDataAnalyzer(PyposmatDataAnalyzer):

    def filter_with__qoi_constraints(self,kde_df,qoi_constraints):
        _df = copy.deepcopy(kde_df)

        for qn in self.datafile.qoi_names:
            aen = "{}.abserr".format(qn)
            en = "{}.err".format(qn)
            _df[aen] = _df[en].abs()

        for qoic_n,qoic_v in qoi_constraints.items():
            nr0,nc0 =_df.shape
            if qoic_v[0] == '<':
                _df = _df[_df[qoic_n] < qoic_v[1]]
            elif qoic_v[0] == '>':
                _df = _df[_df[qoic_n] > qoic_v[1]]
            elif qoic_v[0] == '=':
                _df = _df[_df[qoic_n] == qoic_v[1]]
            else:
                raise ValueError('unknown operator, {}'.format(k))
            nr1,nc0 =_df.shape
            s = "{:>20} {:^3} {:<10} {:10} {:10} {:10}".format(qoic_n,qoic_v[0],qoic_v[1],nr1,nr0,nr0-nr1)
            print(s)
        return _df

    def write_kde_file(self,filename):
        _qoi_constraints = self.configuration.qoi_constraints
        print__qoi_constraints(_qoi_constraints)
        kde_df = copy.deepcopy(self._df)
        for k,v in _qoi_constraints.items():
            if k == 'qoi_constraints':
                kde_df = self.filter_with__qoi_constraints(kde_df,v)
                kde_df = kde_df.reset_index(drop=True)
            elif k == 'filter_by_qoi_error':
                kde_df = self.filter_performance_requirements(kde_df,v)
                kde_df = kde_df.loc[kde_df['is_survive'] == 1]
                kde_df = kde_df.reset_index(drop=True)
            elif k == 'filter_by_phase_order':
                kde_df = self.filter_phase_order(kde_df,v)
            elif k == 'filter_by_pareto' or k == 'select_pareto_only':
                if v == True:
                    kde_df = self.calculate_pareto_set(kde_df,v)
                    kde_df = kde_df.loc[kde_df['is_pareto'] == 1]
                    kde_df = kde_df.reset_index(drop=True)
            elif k == 'filter_by_dmetric':
                    (nr,nc) = kde_df.shape
                    kde_df = self.calculate_d_metric(kde_df)
                    (nr,nc) = kde_df.shape
                    if v[1] == 'pct':
                        pct_to_keep = v[0]/100
                        n = int((nr * pct_to_keep) //1)
                    else:
                        n = min(v[0],nr)
                    kde_df = kde_df.nsmallest(n,'d_metric')
                    for en in self.error_names:
                        print(en,kde_df[en].max())
                    (nr,nc) = kde_df.shape
            else:
                raise ValueError("unknown qoi_constraint method {}".format(k))
            (nr,nc) = kde_df.shape
            print('after {}: {} remainings'.format(k,nr))
        names = ['sim_id'] \
                + self.parameter_names \
                + self.qoi_names\
                + self.error_names
        types = ['sim_id'] \
                + len(self.parameter_names)*['param']\
                + len(self.qoi_names)*['qoi']\
                + len(self.error_names)*['err']
        str_list = []
        str_list.append(','.join(names))
        str_list.append(','.join(types))
        for i_row, row in kde_df[names].iterrows():
            str_list.append(','.join([str(v) for v in row.values.tolist()]))

        with open(filename,'w') as f:
            f.write("\n".join(str_list))
if __name__ == "__main__":
    _fn_config=os.path.join("resources","pyposmat.config.in")
    _fn_data=os.path.join("resources","pyposmat.results.0a.out")
    _fn_pareto_out=os.path.join("pyposmat.pareto.out")
    _fn_kde_out=os.path.join("pyposmat.kde.out")

    pda = Dev__PyposmatDataAnalyzer(
        fn_config=_fn_config,
        fn_data=_fn_data)
    pda.write_kde_file(filename=_fn_kde_out)
