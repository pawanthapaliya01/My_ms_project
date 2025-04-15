#!/usr/bin/env python

import numpy as np

import matplotlib                 # module for making pretty figures
matplotlib.use('Agg')             # use 'Agg' backend - do not show plot, just save figure

from matplotlib import rcParams   # class to modify figure properties
#rcParams.update({'figure.autolayout': True}) # autolayout
rcParams['legend.numpoints'] = 1  # set 1 marker in legend (instead of two)

import matplotlib.pyplot as plt   # module for making pretty figures

##################################################################################################

# Figure of nbub(F), bubsize(F) for L60B36

# set figure properties

fig = plt.figure(figsize=(3.8,3))   # makes a 3.8in x 3in figure
ax1 = fig.add_subplot(111)

# load data for nbub(T)
data = np.loadtxt('nbub_force_py.dat')
nbub, = ax1.plot(data[:,0], data[:,1], ms = 2.5, marker = 'o', color='Blue', mec='Blue', mfc='Blue',
              mew=0.5, lw=0.5, label=r'$n_b$', clip_on=False, zorder=100)

ax2 = ax1.twinx()

# load data for theta
data = np.loadtxt('bubsize_force_py.dat')
bubsize, = ax2.plot(data[:,0], data[:,1], ms = 2.5, marker = 'o', color='DarkCyan', mec='DarkCyan', mfc='White',
                  mew=0.5, lw=0.5, label=r'$l_b$', clip_on=False, zorder=100)


plt.legend(handles=[nbub, bubsize], bbox_to_anchor=(0.35,0.6), loc=1, fontsize=10,
           handletextpad=0.4, borderpad=0.5, handlelength=1.5, labelspacing=0.5)
plt.title('$\mathrm{L33B9}$ \n',fontsize=10)         # text in mathmode using \mathrm{}
                                                      # to produce Times New Roman font 

ax1.set_xlabel('$\mathrm{F}_{\,}(\mathrm{pN})$',fontsize=10)  # x-axis
ax1.set_xlim((0, 40))

ax1.set_ylabel('$n_b$',fontsize=10)           # left y-axis
ax1.set_ylim((0.0,6.8))

ax2.set_ylabel('$l_b\,(\mathrm{bp})$',fontsize=10)           # right y-axis
ax2.set_ylim((0,60))

ax1.tick_params(which='major', length=6, labelsize=8, top=True, direction='in')
ax2.tick_params(which='major', length=6, labelsize=8, top=True, direction='in')

plt.tight_layout()
plt.savefig('L33B9_force_bub.eps', format='eps', dpi=300)      # save the figure

plt.close()
