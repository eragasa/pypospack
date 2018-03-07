import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
def str__sampling_configuration(sampling_config):
    _n_iterations = sampling_config['n_iterations']
    _mc_seed = sampling_config['n_iterations']
    _sampling_table = [sampling_config[i] for i in range(_n_iterations)]

    s = []
    s.append(80*"-")
    s.append('{:^80}'.format('SAMPLING INFORMATION'))
    s.append(80*"-")
    s.append('{:^3} {:^10} {:^10} {:^50}'.format('i','type','n_samples','filename'))
    s.append('{} {} {} {}'.format(2*'-',10*'-',10*'-',50*'-'))
    for i,v in enumerate(_sampling_table):
        _type=v['type']
        _nsamples = v['n_samples']
        
        _filename=''
        if 'filename' in v:
            _filename =v['filename']
        s.append('{:^2} {:^10} {:^10} {:^50}'.format(i,_type,_nsamples,_filename))
    s.append('{} {} {} {}'.format(2*'-',10*'-',10*'-',50*'-'))
    s.append("n_iterations: {}".format(sampling_config['n_iterations']))
    s.append("mc_seed:      {}".format(sampling_config['mc_seed']))
    return "\n".join(s)

def str__qoi_information(qois):
    QOI_LINE_FORMAT="{:^14} {:^14} {:^14} {:^14} {:^14}"

    s = []
    s.append(80*'-')
    s.append('{:^80}'.format('QOI INFORMATION'))
    s.append(80*'-')
    s.append(QOI_LINE_FORMAT.format(
        'qoi_name','qoi_type','target','structure','structure_name'))
    s.append(QOI_LINE_FORMAT.format(14*'-',14*'-',14*'-',14*'-',14*'-'))
    for k,v in qois.items():
        qoi_name=k
        qoi_type=v['qoi_type']
        qoi_target=v['target']
        qoi_s = [ks for ks in v['structures']]
        qoi_n = [kv for ks,kv in v['structures'].items()]
        for i in range(len(qoi_s)):
            if i == 0:
                s.append(QOI_LINE_FORMAT.format(
                    qoi_name,qoi_type,qoi_target,qoi_s[i],qoi_n[i]))
            else:
                s.append(QOI_LINE_FORMAT.format(
                    '','','',qoi_s[i],qoi_n[i]))
    s.append(QOI_LINE_FORMAT.format(14*'-',14*'-',14*'-',14*'-',14*'-'))
    return '\n'.join(s)

def analyze_pypospack_iterative_sampler(data_dir):
    filename_config=os.path.join(
            data_dir,'pyposmat.config.in')
    configuration=PyposmatConfigurationFile()
    configuration.read(filename=filename_config)

    n_iterations = configuration.sampling_type['n_iterations']
    print(str__sampling_configuration(configuration.sampling_type))
    print(str__qoi_information(configuration.qois))
    

    data_results=[]
    for i in range(n_iterations):
        filename_results=os.path.join(
                data_dir,
                'pyposmat.results.{}.out'.format(i))
        print('reading {}...'.format(filename_results))
        data_results.append(PyposmatDataAnalyzer())
        data_results[i].read_configuration_file(
                filename=filename_config)
        data_results[i].read_data_file(
                filename=filename_results)
        data_results[i].write_kde_file(filename='temp.out')  
    data_kde=[]
    for i in range(n_iterations):
        filename_kde=os.path.join(
                data_dir,
                'pyposmat.kde.{}.out'.format(i))
        print('reading {}...'.format(filename_kde))
        data_kde.append(PyposmatDataAnalyzer())
        data_kde[i].read_configuration_file(
                filename=filename_config)
        try:
            data_kde[i].read_data_file(
                    filename=filename_kde)
        except FileNotFoundError as e:
            data_kde[i] = None
    
    n_simulation_info = OrderedDict()
    for i in range(n_iterations):
        _total=configuration.sampling_type[i]['n_samples']
        try:
            _results,ncols=data_results[i].df.shape
            _kde,ncols=data_kde[i].df.shape
        except AttributeError as e:
            _kde=0
        
        n_simulation_info[i] = OrderedDict()
        n_simulation_info[i]['total']=_total
        n_simulation_info[i]['results']=_results
        n_simulation_info[i]['kde']=_kde

    s = []
    s.append('{:^5} {:^10} {:^10} {:^10}'.format(
        'i','total','results','kde'))
    s.append('{:^5} {:^10} {:^10} {:^10}'.format(
        5*'-',10*'-',10*'-',10*'-'))
    for i in range(n_iterations):
        s.append('{:^5} {:>10} {:>10} {:>10}'.format(
            i,
            n_simulation_info[i]['total'],
            n_simulation_info[i]['results'],
            n_simulation_info[i]['kde']))
    print('\n'.join(s))
    qoi_targets = configuration.qoi_targets


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir",
            action='store',
            dest='data_dir',
            type=str,
            help="location of the data directory")
    parse_args = parser.parse_args()

    data_dir=parse_args.data_dir

    analyze_pypospack_iterative_sampler(data_dir)
    exit()
     
    _filename_pyposmat_data = 'data/pyposmat.kde.10.out'
    _n_potentials = 30
    datafile=PyposmatDataFile(filename=_filename_pyposmat_data)
    datafile.read()
    datafile.qoi_references = OrderedDict()
    datafile.qoi_references['TARGET'] = copy.deepcopy(qoi_targets)
    datafile.score_by_d_metric(scaling_factors='TARGET')
    datafile.subselect_by_score(
            score_name='d_metric',
            n=_n_potentials)
    _filename_subselect = datafile.write_subselect()

    datafile=PyposmatDataFile(filename=_filename_subselect)
    datafile.read()
    error_names = datafile.error_names
    qoi_names = datafile.qoi_names
    for i_error,n_error in enumerate(error_names):
        _qoi_name = qoi_names[i_error]
        datafile.df[n_error] = datafile.df[n_error]/qoi_targets[_qoi_name]

    (_nrows,_ncols) = datafile.df.shape
    import matplotlib.pyplot as plt
    print(datafile.df[error_names])
    fig, ax = plt.subplots()
    for i_error,n_error in enumerate(error_names):
        _yloc = [i_error+1]
        ax.scatter(
                datafile.df[n_error],
                _nrows*[i_error+1],
                marker='|',
                s=100.,
                color='k')
    plt.sca(ax)
    plt.yticks(range(len(qoi_names)+1),['']+qoi_names)
    fig.savefig('rugplots.png')

