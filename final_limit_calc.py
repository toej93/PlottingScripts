import numpy as np
import sys
import matplotlib.pyplot as plt
from pylab import setp
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
from scipy.interpolate import splrep, splev

import constants as const
import tools as tool
import plots as plotter

#import plots as plotter
# the livetime from each configuration, this we know by counting up how many good seconds we had
livetime = np.array([15518810,12404200,8160666,35390429,22788969])

# energy bins we have in units of log10(energy)
energy_bins = np.array([17,18,19,20])

# the effective volume (in cubic kilometers * steradians), which we know from Jorge calculating them for us
veff_steradian_config1 = np.array([0.45, 5.06, 23.64, 59.62])
veff_steradian_config2 = np.array([0.39, 5.84, 23.28, 59.27])
veff_steradian_config3 = np.array([0.43, 5.97, 23.47, 59.44])
veff_steradian_config4 = np.array([0.41, 5.06, 23.64, 59.90])
veff_steradian_config5 = np.array([0.41, 5.06, 23.64, 59.90]) #assume same as c5 until Jorge can figure out his problem w/ c5 sim

# convert to km^3 to cm^3
veff_steradian_config1*=1e15
veff_steradian_config2*=1e15
veff_steradian_config3*=1e15
veff_steradian_config4*=1e15
veff_steradian_config5*=1e15

print(veff_steradian_config1)

# convert to effective area by dividing by interaction length
aeff_steradian_config1=veff_steradian_config1/tool.get_Lint(np.power(10.,energy_bins))
aeff_steradian_config2=veff_steradian_config1/tool.get_Lint(np.power(10.,energy_bins))
aeff_steradian_config3=veff_steradian_config1/tool.get_Lint(np.power(10.,energy_bins))
aeff_steradian_config4=veff_steradian_config1/tool.get_Lint(np.power(10.,energy_bins))
aeff_steradian_config5=veff_steradian_config1/tool.get_Lint(np.power(10.,energy_bins))

aeff_interpolator = splrep(energy_bins, np.log10(aeff_steradian_config1))

eff_energy_bins = np.array([17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25])
eff_config1 = np.array([0.34, 0.42, 0.47, 0.56, 0.58, 0.65, 0.69])
eff_interpolator = splrep(eff_energy_bins, eff_config1)

counts_koteramax=[]
energy_bins_counts=[]

bins = np.arange(17,20,0.5)
for bin in bins:
	print(bin)
	temp_logev = np.arange(bin,bin+0.5,0.1)
	temp_energy = np.power(10.,temp_logev)
	temp_aeff = np.power(10.,splev(temp_logev, aeff_interpolator))
	temp_eff = splev(temp_logev,eff_interpolator)
	temp_koteramax = tool.get_flux('kotera_max',temp_logev)
	temp_counts_koteramax = np.trapz(temp_koteramax*temp_aeff*livetime[0]*temp_eff,temp_energy)
	counts_koteramax.append(temp_counts_koteramax)
	energy_bins_counts.append(np.power(10.,bin))

fig = plt.figure(figsize=(2.*11,8.5))
ax_veff = fig.add_subplot(1,2,1)
ax_counts = fig.add_subplot(1,2,2)
ax_veff.plot(pow(10.,energy_bins),veff_steradian_config1,'-o',linewidth=2.0,color='blue',label=r'Config 1')

n_koteramax, bins, patches= ax_counts.hist(energy_bins_counts,
									bins=np.power(10.,np.arange(15,22,0.5)),
									weights=counts_koteramax,
									label=r'Kotera Max: %.3f'%np.sum(counts_koteramax),
									fill=False, 
									stacked=True, 
									histtype='step', 
									edgecolor='blue',
									linewidth=4)

plotter.beautify_veff(ax_veff)
plotter.beautify_counts(ax_counts)
ax_counts.set_ylabel('Events',size=17) #modify this axis title from plotter default
fig.savefig("example.png",edgecolor='none',bbox_inches="tight") #save the figure
