# this script produces figure 1 from the PRL article
import os,sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyflamestk.pyposmat as pyposmat
import pyflamestk.pareto as pareto
import pyflamestk.pyposmatpost as post
# import pyflamestk.stats as stats
x_axis_err_name = 'MgO_NaCl.c12.err'
y_axis_err_name = 'MgO_NaCl.c44.err'
x_axis_latex_name = '$|$error in $C_{12}|$ [GPa]'
y_axis_latex_name = '$|$error in $C_{44}|$ [GPa]'


err_names_pftk = [x_axis_err_name,y_axis_err_name]
err_names_latex = [x_axis_latex_name,y_axis_latex_name]
err_range_plot = {}
err_range_plot['MgO_NaCl.c12.err'] = [0,100]  # min,max of the range to plot
err_range_plot['MgO_NaCl.c44.err'] = [0,100] # min,max of the range to plot

# the directory from which we get data
data = '../data/output'
dpi = 800
fig_fname = 'fig_2.png'

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
mpl.rc('font', **font)

class MplColorHelper(object):
    """
    Req:
       import matplotlib as mpl
       import matplotlib.pyplot as plt
    Ref:
        http://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
    """
    def __init__(self, cmap_name, start_val, stop_val):
        self.cmap_name = cmap_name
        self.cmap = plt.get_cmap(cmap_name)
        self.norm = mpl.colors.Normalize(vmin=start_val,vmax=stop_val)
        self.scalar_map = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self,val):
        return self.scalar_map.to_rgba(val)

def get_rgb_from_cmap(c_val,cmap_name,start_val,stop_val):
    """
    Req:
       import pyflamestk as pftk

    Ref:
        http://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
    """
    mpl_color_helper = MplColorHelper(cmap_name,start_val,stop_val)
    return mpl_color_helper.get_rgb(c_val)

rgb_results = get_rgb_from_cmap(c_results,cmap_name,cmap_min,cmap_max)
rgb_pareto = get_rgb_from_cmap(c_pareto,cmap_name,cmap_min,cmap_max)
rgb_culled = get_rgb_from_cmap(c_culled,cmap_name,cmap_min,cmap_max)

print('rgb_results = {}'.format(rgb_results))
print('rgb_pareto = {}'.format(rgb_pareto))
print('rgb_culled = {}'.format(rgb_culled))

free_params = ['chrg_Mg','MgO_A', 'MgO_rho', 'OO_A', 'OO_rho', 'OO_C']
qoi_ref_values = {'MgO_NaCl.a0': 4.246,
                  'MgO_NaCl.c11': 277.00031,
                  'MgO_NaCl.c12': 91.67016,
                  'MgO_NaCl.c44': 144.00722,
                  'MgO_NaCl.B': 153.4468767,
                  'MgO_NaCl.G': 92.665075,
                  'MgO_NaCl.fr_a': 10.9781666,
                  'MgO_NaCl.fr_c': 8.98642095,
                  'MgO_NaCl.sch':5.067179685,
                  'MgO_NaCl.001s': 0.055950069}
qoi_ref_latex = {'MgO_NaCl.a0':'$a_0$',
                 'MgO_NaCl.c11':'$c_{11}$',
                 'MgO_NaCl.c12':'$c_{12}$',
                 'MgO_NaCl.c44':'$c_{44}$',
                 'MgO_NaCl.B':'$B$',
                 'MgO_NaCl.G':'$G$',
                 'MgO_NaCl.fr_a':'$E_{fr,a}$',
                 'MgO_NaCl.fr_c':'$E_{fr,c}$',
                 'MgO_NaCl.sch':'$E_{sch}$',
                 'MgO_NaCl.001s':'$\gamma_{001}$'}

n_iterations = 10

# these are the formats of the files generated by pyposmat
fname_results_format = "results_{:03d}.out"
fname_pareto_format = "pareto_{:03d}.out"
fname_culled_format = "culled_{:03d}.out"

print('creating iteration results....')
sr = []
for i in range(n_iterations):
    fname_results = os.path.join(data,
                                 fname_results_format.format(i))
    fname_pareto = os.path.join(data,
                                fname_pareto_format.format(i))
    fname_culled = os.path.join(data,
                                fname_culled_format.format(i))
    print('reading_iteration: ({}) '.format(i),
          '{}'.format(fname_results),
          '{}'.format(fname_pareto),
          '{}'.format(fname_culled))
    sr.append(pareto.SimulationResults())
    sr[i].read_simulation_results(fname_results,
                                  fname_pareto,
                                  fname_culled)

print('creating post-processor engine')
# create post processor engine
srpp = []
for i in range(n_iterations):
    srpp.append(post.SimulationResultsPostProcessor(sr[i]))
    
for i in range(n_iterations):
    srpp[i].free_parameter_names = free_params

for i in range(n_iterations):
    srpp[i].qoi_reference_values = qoi_ref_values


print('create {}'.format(fig_fname))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 8}
mpl.rc('font', **font)

pair_names = err_names_pftk

# plotting all the results
# = plt.subplot(0,0,1)
fig = plt.figure()
ax1 = plt.scatter(srpp[0].get_data_frame(dt='abs_err',ds='results')[pair_names[0]],
                          srpp[0].get_data_frame(dt='abs_err',ds='results')[pair_names[1]],
                          s=1, color=rgb_results,label='results')
ax2 = plt.scatter(srpp[0].get_data_frame(dt='abs_err',ds='pareto')[pair_names[0]],
                          srpp[0].get_data_frame(dt='abs_err',ds='pareto')[pair_names[1]],
                          s=1, color=rgb_pareto,label='pareto')
ax3 = plt.scatter(srpp[0].get_data_frame(dt='abs_err',ds='culled')[pair_names[0]],
                          srpp[0].get_data_frame(dt='abs_err',ds='culled')[pair_names[1]],
                          s=1, color=rgb_culled,label='culled')
#fig.legend(handles=ax)
#ax.legend(handles=[ax1,ax2,ax3])
plt.legend(loc='upper right')
plt.xlim(err_range_plot[err_names_pftk[0]][0],
         err_range_plot[err_names_pftk[0]][1])
plt.ylim(err_range_plot[err_names_pftk[1]][0],
         err_range_plot[err_names_pftk[1]][1])
plt.xlabel(err_names_latex[0])
plt.ylabel(err_names_latex[1])

fig.savefig(fig_fname,dpi=800)