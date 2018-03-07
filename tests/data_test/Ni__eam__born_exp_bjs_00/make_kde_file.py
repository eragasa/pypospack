import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename_config",
            action='store',
            dest='filename_config',
            type=str,
            default='data/pyposmat.config.in')
    parser.add_argument("--filename_in",
            action='store',
            dest='filename_in',
            type=str,
            default='data/pyposmat.kde.10.out')
    parser.add_argument("--filename_out",
            action='store',
            dest='filename_out',
            type=str,
            default='pyposmat.kde.0.out')
    parser.add_argument("--n_potentials",
            action='store',
            dest='n_potentials',
            type=int,
            default=1000)
    parse_args = parser.parse_args()

    _filename_config=parse_args.filename_config
    _filename_in=parse_args.filename_in
    _filename_out=parse_args.filename_out
    _n_potentials=parse_args.n_potentials

    config=PyposmatConfigurationFile()
    config.read(filename=_filename_config)
    qoi_targets=config.qoi_targets

    datafile=PyposmatDataFile(filename=_filename_in)
    datafile.read()
    datafile.qoi_references=OrderedDict()
    datafile.qoi_references['TARGET']=copy.deepcopy(qoi_targets)
    datafile.score_by_d_metric(scaling_factors='TARGET')
    datafile.subselect_by_score(score_name='d_metric',n=_n_potentials)
    datafile.write_subselect(filename=_filename_out)
