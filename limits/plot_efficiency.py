import numpy as np
import sys
import matplotlib.pyplot as plt
from pylab import setp
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
from scipy.interpolate import splrep, splev
import pandas as pd
from limits import LimitFigure, ara_energies, ara_available, ara_projected, efficiencies

import constants as const
import tools as tool
import plots as plotter

###################################
# Just Plot our efficiencies
###################################

eff_ebins = np.array([17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25])
eff_c1 = np.array([0.34, 0.42, 0.48, 0.56, 0.58, 0.66, 0.71])
eff_c2 = np.array([0.25, 0.34, 0.40, 0.44, 0.48, 0.60, 0.64])
eff_c3 = np.array([0.33, 0.43, 0.49, 0.58, 0.63, 0.67, 0.87])
eff_c4 = np.array([0.29, 0.38, 0.46, 0.55, 0.58, 0.65, 0.69])
eff_c5 = np.array([0.31, 0.40, 0.46, 0.54, 0.60, 0.67, 0.80])

frac_uptime = np.array([0.165, 0.132, 0.087, 0.375, 0.242])
eff_ave = (frac_uptime[0]*eff_c1) + (frac_uptime[1]*eff_c2) + (frac_uptime[2]*eff_c3) + (frac_uptime[3]*eff_c4) + (frac_uptime[4]*eff_c5)

fig = plt.figure(figsize=(8.5,8.5))
ax = fig.add_subplot(1,1,1)

ax.plot(np.power(10.,eff_ebins), eff_c1, '-^', color='blue', label='Config 1',alpha=0.5)
ax.plot(np.power(10.,eff_ebins), eff_c2, '-v', color='orange', label='Config 2',alpha=0.5)
ax.plot(np.power(10.,eff_ebins), eff_c3, '->', color='green', label='Config 3',alpha=0.5)
ax.plot(np.power(10.,eff_ebins), eff_c4, '-<', color='cyan', label='Config 4',alpha=0.5)
ax.plot(np.power(10.,eff_ebins), eff_c5, '-s', color='magenta', label='Config 5',alpha=0.5)
ax.plot(np.power(10.,eff_ebins), eff_ave, '-o', color='black', label='Average', linewidth=3)

sizer=20
xlow =  1.e17 #the lower x limit
xup = 1e21 #the uppper x limit
ylow = 0 #the lower y limit
yup = 1.05 #the upper y limit
ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
ax.set_ylabel('Efficiency',size=sizer)
# ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(labelsize=sizer)
ax.set_xlim([xlow,xup]) #set the x limits of the plot
ax.set_ylim([ylow,yup]) #set the y limits of the plot
ax.grid()
this_legend = ax.legend(loc='upper left')
setp(this_legend.get_texts(), fontsize=20)
setp(this_legend.get_title(), fontsize=20)
fig.savefig("final/efficiencies.png",edgecolor='none',bbox_inches="tight", dpi=300) #save the figure

###################################
# And compare them to others
###################################

data_tb = np.genfromtxt("testbed_efficiency.csv",delimiter=',',skip_header=1,names=['energy_logev','eff'])
tb_logeV = data_tb['energy_logev']
tb_eff = data_tb['eff']

thomas_energy = np.array([1.015313e+17, 3.114584e+17, 1.000000e+18, 3.174680e+18, 1.019295e+19, 3.235192e+19, 1.015313e+20, 3.149923e+20])
thomas_eff =np.array([0.476333, 0.492164, 0.571537, 0.567818, 0.609883, 0.681175, 0.690122, 0.741392])

fig2 = plt.figure(figsize=(8.5,8.5))
ax2 = fig2.add_subplot(1,1,1)

ax2.plot(thomas_energy, thomas_eff, '-', color='green', label='Previous A23 Analysis', linewidth=3)
# ax2.plot(thomas_energy, thomas_eff, '-', color='green', label='Previous A23 Analysis [10.1103/93.082003]', linewidth=3)
ax2.plot(np.power(10.,eff_ebins), eff_ave, '-o', color='black', label='This Analysis', linewidth=3)
# ax2.plot(np.power(10.,tb_logeV[2:-6]), tb_eff[2:-6], '-', color='blue', label='Testbed Analysis [10.1016/2015.04.006]', linewidth=3)
ax2.plot(np.power(10.,tb_logeV[2:-6]), tb_eff[2:-6], '-', color='blue', label='Testbed Analysis', linewidth=3)


sizer=20
xlow =  1.e17 #the lower x limit
xup = 1e21 #the uppper x limit
ylow = 0 #the lower y limit
yup = 1.05 #the upper y limit
ax2.set_xlabel('Energy [eV]',size=sizer) #give it a title
ax2.set_ylabel('Efficiency',size=sizer)
# ax.set_yscale('log')
ax2.set_xscale('log')
ax2.tick_params(labelsize=sizer)
ax2.set_xlim([xlow,xup]) #set the x limits of the plot
ax2.set_ylim([ylow,yup]) #set the y limits of the plot
ax2.grid()
this_legend2 = ax2.legend(loc='upper left')
setp(this_legend2.get_texts(), fontsize=20)
setp(this_legend2.get_title(), fontsize=20)
fig2.savefig("final/compare_efficiencies.png",edgecolor='none',bbox_inches="tight",dpi=300) #save the figure    
