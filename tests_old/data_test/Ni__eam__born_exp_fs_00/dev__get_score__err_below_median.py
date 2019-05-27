import collection as OrderedDict
import pandas as df

def get_closest_to_orgin(df,N):
    assert isinstance(df,pd.DataFrame)

def get_pareto_set(df,N):
    assert isinstance(df,pd.DataFrame)

def get_error_constraints_by_max_percentage_error(qoi_ref,qoi_pct_err):
    assert isinstance(qoi_ref,OrderedDict)
    
    _qoi_pct_err = OrderedDict()
    if isinstance(qoi_pct_err,OrderedDict):
        _qoi_pct_err = OrderedDict(qoi_pct_err)
    elif type(qoi_pct_err) in [int,float]:
        for qn,qv in qoi_ref.items():
            _qoi_pct_err[qn] = qoi_pct_err
    else:
        err_msg = "qoi_pct_err must either be OrderedDict or a numeric value"
        raise ValueError(err_msg)

    _error_constraints = OrderedDict()
    for qn,qv in qoi_ref.items():
        _qce = _qoi_pct_err[qn]
        _error_constraints[qn] = ['<',abs(qv) * _qce]

    return _error_constraints

def apply_error_constraints(df,
        error_constraints):
    for ec in error_constraints.iterrows

def apply_parameter_constraints(df,
        parameter_constraints):

def get_score_by_simulation(df,
        name_list, 
        err_constraints=None,
        parameter_constraints=None):
    """
    df(pandas.DataFrame):
    err_constraints(OrderedDict)
    """
    assert isinstance(df,pd.DataFrame)
    assert isinstance(error_constraints)
    avg_val_dict = {}
    row_score_list = []
    for names in name_list:
        mean = df[names].median()
        avg_val_dict[names] = mean
    for i in range(0, len(df)):
        row = df.iloc[i]
        row_score = 0
        for names in name_list:
            value = row[names]
            if value < avg_val_dict[names]:
                row_score += 1
        row_score_list.append(row_score)

    df['score'] = row_score_list
    return row_score_list
