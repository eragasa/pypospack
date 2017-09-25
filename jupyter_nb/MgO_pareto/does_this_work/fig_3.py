# this script produces figure 3 from the PRL article
import os,sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pandas.tools.plotting import parallel_coordinates

import pyflamestk.pyposmat as pyposmat
import pyflamestk.pareto as pareto
import pyflamestk.pyposmatpost as post

# THINGS TO DO:
# 1.  Move some of these functions into the pyflamestk module
#     a.  qoi_ref_values, should be read in from the configuration file
#     b.  absolute_scaling factor should be part of SimulationResults
# 2.  Move the reference data into input/data
# description of global variables

absolute_scale = None    # to use absolute scale

free_params = ['chrg_Mg','MgO_A', 'MgO_rho', 'OO_A', 'OO_rho', 'OO_C']
x_axis_err_name = 'MgO_NaCl.c12.err'
y_axis_err_name = 'MgO_NaCl.c44.err'
x_axis_latex_name = '$|$error in $C_{12}|$ [GPa]'
y_axis_latex_name = '$|$error in $C_{44}|$ [GPa]'

qoi_ref_dft    = {'MgO_NaCl.a0': 4.246,
                  'MgO_NaCl.c11': 277.00031,
                  'MgO_NaCl.c12': 91.67016,
                  'MgO_NaCl.c44': 144.00722,
                  'MgO_NaCl.B': 153.4468767,
                  'MgO_NaCl.G': 92.665075,
                  'MgO_NaCl.fr_a': 10.9781666,
                  'MgO_NaCl.fr_c': 8.98642095,
                  'MgO_NaCl.sch':5.067179685,
                  'MgO_NaCl.001s': 0.055950069}


qoi_ref_names = ['chrg_Mg','chrg_O',\
                 'MgMg_A','MgMg_rho','MgMg_C',
                 'MgO_A','MgO_rho','MgO_C',
                 'OO_A','OO_rho','OO_C',
                 'MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44',
                 'MgO_NaCl.B','MgO_NaCl.G',
                 'MgO_NaCl.fr_a','MgO_NaCl.fr_c','MgO_NaCl.sch',
                 'MgO_NaCl.001s',
                 'MgO_NaCl.a0.err','MgO_NaCl.c11.err','MgO_NaCl.c12.err',
                 'MgO_NaCl.c44.err','MgO_NaCl.B.err','MgO_NaCl.G.err',
                 'MgO_NaCl.fr_a.err','MgO_NaCl.fr_c.err','MgO_NaCl.sch.err',
                 'MgO_NaCl.001s.err']

qoi_ref_LC = [2.0,-2.0,0.0,0.5,0.0,821.6,0.3242,0.0,22764.0,0.149,27.88,4.21078276561128,307.5718095128,171.135602331774,168.168424521017,216.61433805878266,68.21810359051298,9.68024989411606,9.810715180656189,5.797051474639375,0.06783817649966246,-0.035217234388720264,30.571809512799973,79.46560233177401,24.158424521016997,63.164338058782675,-24.441896409487015,-1.2977501058839405,0.8247151806561881,0.7300514746393745,0.011888176499662464]
qoi_ref_BG1 = [2.0,-2.0,0.0,0.5,0.0,1279.69,0.29969,0.0,9547.96,0.21916,32.0,4.20923604431415,383.274119165401,169.434215310753,179.601185701851,240.71418326230233,106.91995192732399,12.419259511088967,11.869114175328832,7.198887069605007,0.08070791160146304,-0.036763955685850114,106.27411916540098,77.76421531075299,35.591185701851,87.26418326230234,14.259951927323996,1.4412595110889672,2.8831141753288314,2.131887069605007,0.02475791160146304]
qoi_ref_BG2 = [1.7,-1.7,0.0,0.5,0.0,929.69,0.29909,0.0,4870,0.2679,77.0,4.222448,301.315822490901,150.827961179874,142.471471673523,200.990581616883,75.2439306555135,10.434727086962994,8.526633932683126,5.509135247188169,0.0692527749868838,-0.02355200000000046,24.315822490900985,59.15796117987399,-1.5385283264769782,47.540581616883,-17.416069344486502,-0.543272913037006,-0.4593660673168749,0.442135247188169,0.0133027749868838]
# calculated using the parameter set of Lewis and Catlow
qoi_ref_lc    = {}
# calculated using the parameter set of Ball and Grimes #1 as referenced by Henkleman
qoi_ref_bg1   = {}
# calculated using the parameter set of Ball and Grimes #2 as references by Henkleman
qoi_ref_bg2   = {}
for i in range(len(qoi_ref_names)):
    v = qoi_ref_names[i]
    qoi_ref_lc[v] = qoi_ref_LC[i]
    qoi_ref_bg1[v] = qoi_ref_BG1[i]
    qoi_ref_bg2[v] = qoi_ref_BG2[i]

qoi_ref_latex = {'MgO_NaCl.a0':r'$a_0$',
                 'MgO_NaCl.c11':r'$c_{11}$',
                 'MgO_NaCl.c12':r'$c_{12}$',
                 'MgO_NaCl.c44':r'$c_{44}$',
                 'MgO_NaCl.B':r'$B$',
                 'MgO_NaCl.G':r'$G$',
                 'MgO_NaCl.fr_a':r'$E_{fr,a}$',
                 'MgO_NaCl.fr_c':r'$E_{fr,c}$',
                 'MgO_NaCl.sch':r'$E_{sch}$',
                 'MgO_NaCl.001s':r'$\gamma_{001}$'}
qoi_ref_values = qoi_ref_dft
qoi_ref_dict = {}
qoi_ref_dict['DFT'] = qoi_ref_dft
qoi_ref_dict['LC'] = qoi_ref_bg1
qoi_ref_dict['BG1'] = qoi_ref_bg2
qoi_ref_dict['BG2'] = qoi_ref_lc
# ---- DATA DIRECTORY ----
data = '../data/output'
dpi = 800
fig_fname = 'fig_3.png'

# ---- CHANGE HERE TO CONTROL COLOR ---
cmap_name = 'Blues' # This uses a color gradient map scheme defined in matplotlib
cmap_min = 0
cmap_max = 100
c_results = 30
c_pareto = 70
c_culled = 100
rgb_results = None # required in the script
rgb_pareto = None
rgb_culled = None

# ---- CHANGE HERE TO CONTROL FRONT -----
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}


#class MplColorHelper(object):
#    """
#    Req:
#       import matplotlib as mpl
#       import matplotlib.pyplot as plt
#    Ref:
#        http://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
#    """
#    def __init__(self, cmap_name, start_val, stop_val):
#        self.cmap_name = cmap_name
#        self.cmap = plt.get_cmap(cmap_name)
#        self.norm = mpl.colors.Normalize(vmin=start_val,vmax=stop_val)
#        self.scalar_map = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
#
#    def get_rgb(self,val):
#        return self.scalar_map.to_rgba(val)

#def get_rgb_from_cmap(c_val,cmap_name,start_val,stop_val):
#    """
#    Req:
#       import pyflamestk as pftk
#
#    Ref:
#        http://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
#    """
#    mpl_color_helper = MplColorHelper(cmap_name,start_val,stop_val)
#    return mpl_color_helper.get_rgb(c_val)

#rgb_results = get_rgb_from_cmap(c_results,cmap_name,cmap_min,cmap_max)
#rgb_pareto = get_rgb_from_cmap(c_pareto,cmap_name,cmap_min,cmap_max)
#rgb_culled = get_rgb_from_cmap(c_culled,cmap_name,cmap_min,cmap_max)

#print('rgb_results = {}'.format(rgb_results))
#print('rgb_pareto = {}'.format(rgb_pareto))
#print('rgb_culled = {}'.format(rgb_culled))

# --- read in values

def rescale_errors(sr, scaling_factors, results = 'culled'):
    """calculates the closest value to the utopia point

    Args:
        sr (pyflamestk.pyposmat.SimulationResults): provides the 
            information required from the Pareto fitting process.
        scale (dict of str:float): the key is the error name, the value is
            a scalar value from which the errors will be divided for the 
            purposes of scaling
        n (int): this will find the closest points to the utopia point
        result (str,optional): should be either results, pareto or culled.
            Default is culled.

    Raises:
        ValueError: if results is of the wrong type
    """
    
    # make a copy of the values, which we will then mutate to find the 
    # solutions
    rescaled_values = None
    if results == 'results':
        rescaled_values = np.copy(sr.results)
    elif results == 'pareto':
        rescaled_values = np.copy(sr.pareto)
    elif results == 'culled':
        rescaled_values = np.copy(sr.culled)
    else:
        str_out = "results must be results, pareto, or culled. {} sent"
        str_out = str_out.format(results)
        raise ValueError(str_out)


    # subselect just the error_columns, error columns
    err_names = sr.err_names # local copy
    names = ['sim_id'] + err_names
    names_idx = [i for i,v in enumerate(sr._names) if v in names]
    rescaled_values = rescaled_values[:,names_idx]    

    # rescale the error values
    for v in err_names:
        i = names.index(v)
        eta = scaling_factors[v]
        rescaled_values[:,i] = rescaled_values[:,i]/eta

    # calculate distance metric
    n_rows, n_columns = rescaled_values.shape
    err_col_idx = [names.index(qen) for qen in err_names]
    distance_metric_vector = np.zeros(n_rows)
    for i in range(n_rows):
        err_vector = rescaled_values[i,err_col_idx]
        dist = np.sqrt(err_vector.dot(err_vector))
        distance_metric_vector[i] = dist

   # calculate the distance metrics

    print('--- REPACKING ERROR INFORMATION ---')
    head_format = '{:10} {:10} {:25} {:10}'
    line_format = '{:10} {:10} {:25} {:10.6f}'
    print(head_format.format('new_col_i', 'old_col_i','col_name','scaling_factor'))
    for i,v in enumerate(names):
        new_col_i = i
        old_col_i = names_idx[i]
        col_name = v

        if v == 'sim_id':
            print(head_format.format(new_col_i, old_col_i, col_name, 'N/A'))
        else:
            scaling_factor = scaling_factors[col_name]
            print(line_format.format(new_col_i, old_col_i, col_name, scaling_factor))


    return distance_metric_vector, rescaled_values

def get_closest_to_utopia_point(d_vector, n = 1):
    """ return an index of the dvector closest to the utopia point
    """
    assert isinstance(d_vector, np.ndarray)
    assert isinstance(n,int)

    idx_shortest = np.where(d_vector.min())
    idx_shortest_pop = np.argpartition(d_vector,n)[:n]

    return idx_shortest, idx_shortest_pop

def get_absolute_scaling_factors(sr):
    """Determines scaling factor

    If q is a quantity of interest and eta is the scaling factor, then the 
    scaled error should be q/eta ~= 1
    Returns:
        sf - (list of float) a vector of scaling factors
    """
    sf = {} 
    for k,v in sr.qoi_ref.items():
        err_name = k + ".err"
        sf[err_name] = v

    return sf

# read in the files from simulation
fname_results = os.path.join(data,'results_009.out')
fname_pareto = os.path.join(data,'pareto_009.out')
fname_culled = os.path.join(data,'culled_009.out')

# not currently used
fname_pyposmat_config = os.path.join(data,'pyposmat.config')
fname_pyposmat_potential = os.path.join(data,'pyposmat.potential')
fname_pyposmat_qoi = os.path.join(data,'pyposmat.qoi')


sr = pareto.SimulationResults()
sr.read_simulation_results(fname_results, fname_pareto, fname_culled)
#sr.read_configuration_files(pyposmat_config = fname_pyposmat_config:
#                            pyposmat_potential = fname_pyposmat_potential,
#                            pyposmat_qoi = fname_pyposmat_qoi)
sr.qoi_ref = qoi_ref_values # need to fix, should be done automatically

# SCALING FACTORS
# - scaling factors renormalize errors so that errors are of approximately the
# same magnitude, the scaling factors chosen here are by division.
absolute_sf = get_absolute_scaling_factors(sr) # get absolute scaling factors
print('---SCALING FACTORS---')
row_format = "{:25}{:15.4f}"
for k,v in absolute_sf.items():
    print(row_format.format(k,v))

# RESCALE ERRORS
distance_metrics, rescaled_values = rescale_errors(sr,absolute_sf,'culled') 
names = ['sim_id'] + sr.err_names
err_names_idx = [names.index(v) for v in sr.err_names]
print(len(names), rescaled_values.shape)
#bp determine the point closest to the utopia point by the distance metric vector
n_potentials = [100]
pot_idx = {}
utopia_sim_id = {}
for n in n_potentials:
    temp, pot_idx[n] = get_closest_to_utopia_point(distance_metrics,n)
    utopia_sim_id[n] = sr.culled[pot_idx[n],0]

n=100
data = rescaled_values[np.ix_(pot_idx[n],err_names_idx)]
xmin = -100
xmax =  100
ymin = -11 
ymax = 0

plt.rc('text', usetex=True)
fig = plt.figure(1)
ax = plt.subplot(111)
my_yticks = [-i-1 for i in range(len(sr.err_names))]
my_yticklabels = []
for v in sr.err_names:
    qoi_name = v.strip('.err')
    qoi_latex = qoi_ref_latex[qoi_name]
    my_yticklabels.append(qoi_latex)
plt.setp([ax], yticks=my_yticks, yticklabels=my_yticklabels)
n_errors = len(sr.err_names)
n_samples = len(pot_idx[n])
ax.axvline(x=0, color='k')
for i,v in enumerate(sr.err_names):
    ax.plot(100*data[:,i],[-i-1]*n_samples,'|', color='k')
    print(v,qoi_ref_dft[v.strip('.err')],
          qoi_ref_lc[v]/qoi_ref_dft[v.strip('.err')],
          qoi_ref_bg1[v]/qoi_ref_dft[v.strip('.err')],
          qoi_ref_bg2[v]/qoi_ref_dft[v.strip('.err')])
    ax.plot(100*qoi_ref_lc[v]/qoi_ref_dft[v.strip('.err')],
            [-i-1],'|', color='g', markersize = 20, mew = 5,
            label='LC')
    ax.plot(100*qoi_ref_bg1[v]/qoi_ref_dft[v.strip('.err')],
            [-i-1],'|', color='b', markersize = 20, mew = 5,
            label='BG+2.0')
    ax.plot(100*qoi_ref_bg2[v]/qoi_ref_dft[v.strip('.err')],
            [-i-1],'|', color='r', markersize = 20, mew = 5,
            label='BG+1.7')
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_xlabel('Percent Difference from Reference Value')

ax.legend(handles=[mpatches.Patch(color='g',label='LC'),
                   mpatches.Patch(color='b',label='BG+2.0'),
                   mpatches.Patch(color='r',label='BG+1.7')],
          loc='upper left')

fig.savefig('fig_3.png',png=800)
#ax.violinplot(data,showmeans = False, showmedians=True)
#fig.savefig('fig_3.png',dpi=800)

# make pandas dataframe of just errors
import pandas as pd

#culled_df = pd.DataFrame(data=data,columns=sr.err_names)
#culled_df[:10].plot(kind='barh',rot=0)

def get_pandas_dataframe(sr, df_type = 'err'):

    col_names =  None
    if df_type is 'err':
        col_names = ['sim_id'] + sr.err_names
    elif df_type is 'param':
        col_names = ['sim_id'] + sr.param_names
    else:
        raise AttributeError('dt_type was {}, must be err, or param'.format(dt_type))
    col_idx = [i for i,v in enumerate(sr._names) if v in col_names]
    np_pareto = np.copy(sr.pareto[:,col_idx])
    return  pd.DataFrame(data = sr.pareto[:,col_idx],
                         columns = col_names)

err_df = get_pandas_dataframe(sr, df_type='err')

#n_qois = len(sr.qoi_names)
#fig_violin, ax_violin = plt.subplots(nrows = 1, ncols = 1)

#fig, axes = plt.s
#fig.savefig(fig_fname,dpi=800)
