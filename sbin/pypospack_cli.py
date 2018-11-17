import argparse
from collections import OrderedDict


parsers = OrderedDict()
parsers['obj_parser'] = argparse.ArgumentParser()
parsers['obj_subparsers'] = parsers['obj_parser'].add_subparsers()

def do_visualization(args):
    print('visualization')
    config_fn = args.configuration
    data_fn = args.data
    print('configuration:{}'.format(config_fn))
    print('data:{}'.format(data_fn))

    from pypospack.pyposmat.visualization.bokeh import PyposmatBokehVisualizer

    o = PyposmatBokehVisualizer()
    o.read_configuration(filename=config_fn)
    o.read_data(filename=data_fn)
    o.start_bokeh_server()
    o.setup_bokeh_frame()
    
    # from pypospack.pyposmat.visualization import start_bokeh_visualization
    # start_bokeh_visualization(config_fn=config_fn,data_fn=data_fn)

parsers['visualization'] = parsers['obj_subparsers'].add_parser('visualization')
parsers['visualization'].add_argument('--configuration','--config')
parsers['visualization'].add_argument('--data')
parsers['visualization'].set_defaults(func=do_visualization)

def check_pyposmat_configuration(args):
    print('checking pyposmat configuration file')
    config_fn = args.configuration
    print('pyposmat_configuration_file:{}'.format(config_fn))

#parsers['check_pyposmat_configuration'] = parsers['obj_subparsers'].add_parser('check_pyposmat_configuration')
#parsers['check_pyposmat_configuration'].add_argument('--configuration')
#parsers['check_pyposmat_configuration'].set_defaults(func=check_pyposmat_configuration)
if __name__ == "__main__":
    #args = parsers['obj_parser'].parse_args()
    #args.func(args)
    parsers['args'] = parsers['obj_parser'].parse_args()
    parsers['args'].func(parsers['args'])

    # Examples
    # python ~/repos/pypospack/sbin/pypospack_cli.py visualization --configuration=data/pyposmat.config.in --data=data/pypospack.kde.3.out

