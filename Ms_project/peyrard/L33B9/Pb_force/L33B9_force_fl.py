import numpy as np

import matplotlib                 # module for making pretty figures
matplotlib.use('Agg')             # use 'Agg' backend - do not show plot, just save figure

from matplotlib import rcParams   # class to modify figure properties
#rcParams.update({'figure.autolayout': True}) # autolayout
rcParams['legend.numpoints'] = 1  # set 1 marker in legend (instead of two)

import matplotlib.pyplot as plt   # module for making pretty figures

##################################################################################################

# Figure of l(F), f(F) for L60B36

# set figure properties
fig = plt.figure(figsize=(3.2,3))   # makes a 3.3in x 3in figure
ax = fig.add_subplot(111)

# load data for l
data = np.loadtxt('l_force_py.dat')
l, = plt.plot(data[:,0], data[:,1], ms = 2.5, marker = 'o', color='Blue', mec='Blue', mfc='Blue',
              mew=0.5, lw=0.5, label=r' $l$', clip_on=False, zorder=100)

# load data for f
data = np.loadtxt('f_force_py.dat')
f, = plt.plot(data[:,0], data[:,1], ms = 2.5, marker = 'o', color='DarkCyan', mec='DarkCyan',
              mfc='White', mew=0.5, lw=0.5, label=r'$f$', clip_on=False, zorder=100)

plt.legend(handles=[l, f], bbox_to_anchor=(0.32,0.90), loc=1, fontsize=10,
           handletextpad=0.4, borderpad=0.5, handlelength=1.5, labelspacing=0.5)
plt.title('$\mathrm{L33B9}$ \n',fontsize=10)                 # text in mathmode using \mathrm{}
                                                              # to produce Times New Roman font 
plt.xlabel('$\mathrm{F}_{\,}(\mathrm{pN})$',fontsize=10) 
plt.xlim((0, 40))
plt.ylim((0,1))
ax.tick_params(which='major', length=6, labelsize=8, right=True, top=True, direction='in')
plt.tight_layout()
plt.savefig('L33B9_force_lf.eps', format='eps', dpi=300)      # save the figure

plt.close()
