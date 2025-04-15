#!/usr/bin/env python

import numpy as np

import matplotlib                 # module for making pretty figures
matplotlib.use('Agg')             # use 'Agg' backend - do not show plot, just save figure

from matplotlib import rcParams   # class to modify figure properties
#rcParams.update({'figure.autolayout': True}) # autolayout
rcParams['legend.numpoints'] = 1  # set 1 marker in legend (instead of two)

import matplotlib.pyplot as plt   # module for making pretty figures

##################################################################################################

# Figure of h(F), theta(F) for L33B9

# set figure properties

fig = plt.figure(figsize=(3.8,3))   # makes a 3.8in x 3in figure
ax1 = fig.add_subplot(111)

# load data for h
data = np.loadtxt('h_force_py.dat')
h, = ax1.plot(data[:,0], data[:,1], ms = 2.5, marker = 'o', color='Blue', mec='Blue', mfc='Blue',
              mew=0.5, lw=0.5, label=r'$h$', clip_on=False, zorder=100)

ax2 = ax1.twinx()

# load data for theta
data = np.loadtxt('theta_force_py.dat')
theta, = ax2.plot(data[:,0], data[:,1], ms = 2.5, marker = 'o', color='red', mec='red', mfc='White',
                  mew=0.5, lw=0.5, label=r'$\theta$', clip_on=False, zorder=100)


plt.legend(handles=[h, theta], bbox_to_anchor=(0.35,0.7), loc=1, fontsize=10,
           handletextpad=0.4, borderpad=0.5, handlelength=1.5, labelspacing=0.5)
plt.title('$\mathrm{L33B9}$ \n',fontsize=10)         # text in mathmode using \mathrm{}
                                                      # to produce Times New Roman font 

ax1.set_xlabel('$\mathrm{F}_{\,}(\mathrm{pN})$',fontsize=10)  # x-axis
ax1.set_xlim((0, 40))

ax1.set_ylabel('$\mathrm{h}\,(\AA)$',fontsize=10)               # left y-axis
ax1.set_ylim((3.3,6.4))

ax2.set_ylabel('$\\theta\,(\mathrm{rad})$',fontsize=10)         # right y-axis
ax2.set_ylim((0.0,0.65))

ax1.tick_params(which='major', length=6, labelsize=8, top=True, direction='in')
ax2.tick_params(which='major', length=6, labelsize=8, top=True, direction='in')

plt.tight_layout()
plt.savefig('L33B9_force_h_theta.eps', format='eps', dpi=300)      # save the figure

plt.close()
