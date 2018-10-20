
import copy,os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class PyposmatParetoRugplot(object):

    def __init__(self):
        """

        Attributes:
            parameter_names (list):  A list of parameter names taken from the configuration file
            qoi_names(list): A list of qoi names taken from the configuration file
            error_names(list): A list of error names taken from the configuration file
            qoi_v_names(list): A list of testing qoi names taken from the configuration file
            error_v_names(list): A list of testing error names taken from the configuration file

        """
        self.data_directory = None
        self.configuration_fn = None
        self.data_fn = None
        self.plot_fn = None

        self.qoi_reference_data = None
        self.logger = None
        
        self.parameter_names = None   
        self.qoi_names = None
        self.error_names = None
        self.qoi_validation_names = None
        self.error_validation_names = None

        self.qoi_excluded_names = None

        self.df = None
        self.data_nrows = None
        self.data_ncols = None

        self.reference_data_marker_size = 300,
        self.reference_data_marker_type = '|',
        self.reference_data_colors = None
        self.reference_data_colors_default = ['red','blue','green']        
        
        self.datum_line_width = 1
        self.datum_line_style = '-'
        self.datum_line_color = 'grey'

        self.qoi_v_separator_width=1
        self.qoi_v_separator_style='-'
        self.qoi_v_separator_color='k'

        self.plot_legend_location='upper left'
        self.plot_format = 'svg'
        self.plot_dpi = 1200
    def _log(self,msg):
        if self.logger is None:
            print(msg)
        else:
            self.logger.log(msg)

    def read_configuration(self,filename=None):
        
        if filename is not None: self.configuration_fn = filename
        _filename = self.configuration_fn
           
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(_filename)

        self.parameter_names = list(self.configuration.parameter_names)
        self.qoi_names = list(self.configuration.qoi_names)
        self.error_names = list(self.configuration.error_names)
        self.qoi_validation_names = list(self.configuration.qoi_validation_names)
        self.error_validation_names = list(self.configuration.error_validation_names)

        self.qoi_targets = self.configuration.qoi_targets
        self.qoi_validation_targets = self.configuration.qoi_validation_targets

    def read_datafile(self,filename=None):

        if filename is not None: self.data_fn = filename
        _filename = self.data_fn

        self.data = PyposmatDataFile()
        self.data.read(_filename)

        self.df = copy.deepcopy(self.data.df)
        
        (_nrows,_ncols) = self.df.shape
        self.data_nrows = _nrows
        self.data_ncols = _ncols

    def make_plot(self,filename,include_qois=True,include_qois_v=False,qoi_excluded_names=None):
        """
        Args:
            include_qois (bool): Including fitting qois in the rugplot
            include_qois_v (bool): Including the testing qois in the rugplot
        """

        self.plot_fn = filename

        if qoi_excluded_names is not None:
            self.qoi_excluded_names = list(qoi_excluded_names)
        else:
            self.qoi_excluded_names = []

        _qoi_names = []
        if include_qois:
            _qoi_names += [q for q in self.qoi_names if q not in self.qoi_excluded_names]
        if include_qois_v:
            _qoi_names += [q for q in self.qoi_validation_names if q not in self.qoi_excluded_names]
        
        self.calculate_normed_errors(qoi_names=_qoi_names)

        try:
            _reference_potential_names = [k for k in self.configuration.reference_potentials]
        except KeyError as e:
            _reference_potentials_names = None

        if _reference_potential_names is not None:
            self.reference_df = self.df.loc[self.df['sim_id'].isin(_reference_potential_names)]
            self.results_df = self.df.loc[~self.df['sim_id'].isin(_reference_potential_names)]
        else:
            self.reference_df = None
            self.results_df = self.df
    
        _fig,_ax = plt.subplots()
        
        self.add_vertical_datum_line_to_plot(ax=_ax)
        self.add_results_to_plot(ax=_ax,qoi_names=_qoi_names)
        self.add_reference_potentials_to_plot(ax=_ax,qoi_names=_qoi_names)
        self.add_horizontal_separator(ax=_ax)
        self.add_x_axis_labels(ax=_ax)
        self.add_y_axis_labels(ax=_ax,qoi_names=_qoi_names)
        
        # set limits
        _ax.set_xlim([-1,1])
        _ax.set_ylim([0,len(_qoi_names)+1])
        
        _ref_data_name = [k for k in self.configuration.reference_potentials]
        _ref_data_colors = self.reference_data_colors_default
        _best_data_name = "BEST"
        _best_data_color = "grey"
        _legend_location = self.plot_legend_location
        _handles = []
        for i,v in enumerate(_ref_data_name):
            _label = v
            _color = _ref_data_colors[i]
            _handles.append(mpatches.Patch(color=_color,label=_label))


        _handles.append(mpatches.Patch(color=_best_data_color,label=_best_data_name))

        _ax.legend(handles=_handles,loc=_legend_location)
        
        _fig.tight_layout()

        _plot_fn = self.plot_fn
        _plot_dpi = self.plot_dpi
        if _plot_fn[-3:] in ['png','eps','svg']:
            _plot_format = _plot_fn[-3:]
        else:
            _plot_fn = plot_fn + ".png"
            self.plot_fn = _plot_fn

            _plot_format = ".png"
       
        #s = ("{}\n  plot_format:{}\n  plot_dpi:{}".format(_plot_fn,_plot_format,_plot_dpi))
        #print(s)
        _fig.savefig(_plot_fn,format=_plot_format,dpi=_plot_dpi)
    
    def calculate_normed_errors(self,df=None,qoi_names=None):
        """

        If a pandas.DataFrame is passed to df, then it will be set as the df attribute 
        for this class.  It willl
        Args:
            df (pandas.DataFrame)
            qoi_names (list) - a list of qoi names.  Default behavior will use all the 
                qois specified in the configuration object

        """
        if df is not None:
            self.df = copy.deepcopy(df)

        if qoi_names is not None:
            _qoi_names = list(qoi_names)
        else:
            _qoi_names = list(self.qoi_names)

        self.normed_error_names = []
        self.normed_error_validation_names = []
        for qn in _qoi_names:
            
            if qn in self.qoi_names:

                en = "{}.err".format(qn)
                nen = "{}.nerr".format(qn)
                self.normed_error_names.append(nen)
                
                
                q = self.qoi_targets[qn]

            elif qn in self.qoi_validation_names:
                en = "{}.err_v".format(qn)
                nen = "{}.nerr_v".format(qn)
                self.normed_error_validation_names.append(nen)

                q = self.qoi_validation_targets[qn]
            else:
                s = 80*"-"+"\n"
                s += "{:^80}\n".format('debugging information')
                s += 80*"-"+"\n"
                s += "qoi_name:{}\n".format(qn)
                s += "qoi_names\n"
                s += "\n".join(["  {}".format(v) for v in _qoi_names])+"\n"
                s += 80*"-"+"\n"
                s += "{:^80}\n".format('debugging information')
                s += "qoi_names\n"
                print(s)
                raise ValueError()

            self.df[nen] = self.df[qn]/q-1
               
    def add_results_to_plot(self,ax,qoi_names):
        for iqn,qn in enumerate(qoi_names):
            (nrows,ncols) = self.results_df.shape
            if qn in self.qoi_names:
                x = self.results_df['{}.nerr'.format(qn)]
            elif qn in self.qoi_validation_names:
                x = self.results_df['{}.nerr_v'.format(qn)]
            else:
                raise ValueError()
            y = nrows*[len(qoi_names)-iqn]
            ax.scatter(x,y,marker='|',color='grey')

    def add_reference_potentials_to_plot(self,
            ax,
            qoi_names,
            reference_data_marker_size = 300,
            reference_data_marker_type = '|',
            reference_data_colors = None):

        _default_data_colors = self.reference_data_colors_default
        
        _reference_potential_names = [k for k in self.configuration.reference_potentials]

        if reference_data_marker_size is not None: 
            self.reference_data_marker_size=reference_data_marker_size
        if reference_data_marker_type is not None:
            self.reference_data_marker_type = reference_data_marker_type
        if reference_data_colors is not None:
            if isinstance(reference_data_colors,dict):
                self.reference_data_colors = copy.deepcopy(reference_data_colors)
            elif isinstance(reference_data_colors,list):
                _default_data_colors = list(reference_data_colors)
                self.reference_data_colors = OrderedDict()
                for i,v in enumerate(_reference_potential_names):
                    self.reference_data_colors[v] = _default_data_colors[i]
        else:
            self.reference_data_colors = OrderedDict()
            for i,v in enumerate(_reference_potential_names):
                self.reference_data_colors[v] = _default_data_colors[i]

        _marker_size = self.reference_data_marker_size
        _marker_type = self.reference_data_marker_type
        _marker_color = self.reference_data_colors
       
        for _,row in self.reference_df.iterrows():
            for iqn,qn in enumerate(qoi_names):
                (nrows,ncols) = self.results_df.shape
                if qn in self.qoi_names:
                    x = row['{}.nerr'.format(qn)]
                elif qn in self.qoi_validation_names:
                    x = row['{}.nerr_v'.format(qn)]
                else:
                    raise ValueError()
                y = [len(qoi_names)-iqn]
                ax.scatter(x,y,
                    s=_marker_size,
                    marker=_marker_type,
                    color=_marker_color[row['sim_id']]
                    )

    def add_x_axis_labels(self,ax):
        # set numbers to percentages
        vals = ax.get_xticks()
        ax.set_xticklabels(['{:,.0%}'.format(x) for x in vals])
    
    def add_y_axis_labels(self,ax,qoi_names):
        try:
            _latex_labels = self.configuration.latex_labels 
            _y_ticks_labels = [_latex_labels[qn]['name'] for qn in qoi_names]

            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = 'Helvetica'
            plt.rcParams['text.usetex'] = True
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
            
        except TypeError as e:
            # use text labels
            _y_ticks_labels = qoi_names

        plt.sca(ax)
        plt.yticks(
                list(range(1,len(qoi_names)+1)),
                list(reversed(_y_ticks_labels))
            )

    def add_horizontal_separator(self,ax,
            qoi_v_separator_width=None,
            qoi_v_separator_style=None,
            qoi_v_separator_color=None):

        if qoi_v_separator_width is not None: self.qoi_v_separator_width = qoi_v_separator_width
        if qoi_v_separator_style is not None: self.qoi_v_separator_style = qoi_v_separator_style
        if qoi_v_separator_color is not None: self.qoi_v_separator_color = qoi_v_separator_color

        _qoi_v_separator_width = self.qoi_v_separator_width
        _qoi_v_separator_style = self.qoi_v_separator_style
        _qoi_v_separator_color = self.qoi_v_separator_color

        _qoi_v_separator_location = len(self.qoi_validation_names) + 0.5
        ax.axhline(
                _qoi_v_separator_location,
                color = _qoi_v_separator_color,
                linestyle = _qoi_v_separator_style,
                linewidth = _qoi_v_separator_width)

    def add_vertical_datum_line_to_plot(self,
            ax,
            datum_line_width=None,
            datum_line_style=None,
            datum_line_color=None):

        if datum_line_width is not None: self.datum_line_width = datum_line_width
        if datum_line_style is not None: self.datum_line_style = datum_line_style
        if datum_line_color is not None : self.datum_line_color = datum_line_color

        _datum_line_width = self.datum_line_width
        _datum_line_style = self.datum_line_style
        _datum_line_color = self.datum_line_color

        ax.axvline(0,
                color=_datum_line_color,
                linestyle=_datum_line_style,
                linewidth=_datum_line_width
        )
