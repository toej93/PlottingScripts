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

def get_veff_plot():
    data = np.genfromtxt("veff.csv",delimiter=',',skip_header=1,names=['energy_logev','veffc1','veffc1err','veffc2','veffc2err','veffc3','veffc3err','veffc4','veffc4err','veffc5','veffc5err'])
    logeV = data['energy_logev']
    veff_c1 = data['veffc1']
    veff_c1_err = data['veffc1err']
    veff_c2 = data['veffc2']
    veff_c2_err = data['veffc2err']
    veff_c3 = data['veffc3']
    veff_c3_err = data['veffc3err']
    veff_c4 = data['veffc4']
    veff_c4_err = data['veffc4err']
    veff_c5 = data['veffc5']
    veff_c5_err = data['veffc5err']


    veff_c1*=1e15
    veff_c1_err*=1e15
    veff_c2*=1e15
    veff_c2_err*=1e15
    veff_c3*=1e15
    veff_c3_err*=1e15
    veff_c4*=1e15
    veff_c4_err*=1e15
    veff_c5*=1e15
    veff_c5_err*=1e15

    fig = plt.figure(figsize=(8.5,8.5))
    ax_veff = fig.add_subplot(1,1,1)
    ax_veff.errorbar(pow(10.,logeV),veff_c1,yerr=veff_c1_err,label=r'Config 1',color='blue')
    ax_veff.errorbar(pow(10.,logeV),veff_c2,yerr=veff_c2_err,label=r'Config 2',color='red')
    ax_veff.errorbar(pow(10.,logeV),veff_c3,yerr=veff_c3_err,label=r'Config 3',color='green')
    ax_veff.errorbar(pow(10.,logeV),veff_c4,yerr=veff_c4_err,label=r'Config 4',color='orange')
    ax_veff.errorbar(pow(10.,logeV),veff_c5,yerr=veff_c5_err,label=r'Config 5',color='gray')


    # ax_veff.errorbar(pow(10.,logeV),veff_c1,yerr=veff_c1_err,'-o',linewidth=2.0,color='blue',label=r'Config 1')
    plotter.beautify_veff(ax_veff)
    fig.savefig("veff.png",edgecolor='none',bbox_inches="tight") #save the figure


def plot_final_limit(make_counts_plot=False): #if true, then make counts plot
    #import plots as plotter
    # the livetime from each configuration, this we know by counting up how many good seconds we had
    livetime = np.array([15518810,12404200,8160666,35390429,22788969])

    # energy bins we have in units of log10(energy)
    energy_bins = np.array([1E17,1E18,1E19,1E20])
    
    # the effective volume (in cubic kilometers * steradians), which we know from Jorge calculating them for us
    config={}
    config["veff_steradian_config1"] = np.array([0.45, 5.06, 23.64, 59.62])
    config["veff_steradian_config2"] = np.array([0.39, 5.84, 23.28, 59.27])
    config["veff_steradian_config3"] = np.array([0.43, 5.97, 23.47, 59.44])
    config["veff_steradian_config4"] = np.array([0.41, 5.06, 23.64, 59.90])
    config["veff_steradian_config5"] = np.array([0.40,4.96, 23.53, 59.20])

    # convert to km^3 to cm^3
    # veff_steradian_config1*=1e15

    #The eficiencies for each config
    # config["eff_config1"] = np.array([0.34, 0.42, 0.47, 0.56, 0.58, 0.65, 0.69])
    #config["eff_config2"] = np.array([0.25,0.34,0.40,0.44,0.48,0.60,0.64])
    #config["eff_config3"] = np.array([0.33,0.43,0.49,0.58,0.63,0.67,0.87])
    #config["eff_config4"] = np.array([0.29,0.38,0.46,0.55,0.58,0.65,0.69])
    #config["eff_config5"] = np.array([0.31,0.40,0.46,0.54,0.60,0.67,0.80])


    #Efficiencies if only integer log10(energy) 
    config["eff_config1"] = np.array([0.34, 0.47, 0.58, 0.69])
    config["eff_config2"] = np.array([0.25,0.40,0.48,0.64])
    config["eff_config3"] = np.array([0.33,0.49,0.63,0.87])
    config["eff_config4"] = np.array([0.29,0.46,0.58,0.69])
    config["eff_config5"] = np.array([0.31,0.46,0.60,0.80])

    combined_aeff=0
    combined_limit=0
    for x in range(1,6):
        config["veff_steradian_config{0}".format(x)]*=1e15
        #config["veff_steradian_config{0}".format(x)]*=config["eff_config{0}".format(x)]#change to config["eff_config{0}".format(x)] whenever you add the rest of the efficiencies
        config["aeff_steradian_config{0}".format(x)]=config["veff_steradian_config{0}".format(x)]/tool.get_Lint(energy_bins)
       # print(tool.get_Lint(energy_bins))
        config["EFE_config{0}".format(x)]=1/(np.log(10)*0.5*livetime[x-1]*config["aeff_steradian_config{0}".format(x)]*config["eff_config{0}".format(x)])
        combined_aeff+=config["aeff_steradian_config{0}".format(x)]*livetime[x-1]*config["eff_config{0}".format(x)]
        
        
    combined_limit=2.44/(np.log(10)*0.5*combined_aeff)
    #aeff_interpolator = splrep(energy_bins, np.log10(aeff_steradian_config1))

    #plot shit

    if __name__=="__main__":
        figure = LimitFigure(e_power=1, xlims=(1e6, 1e11), ylims=(1e-19, 2e-13), font_size=16, tick_size=14)
        figure.build_base_plot('ara')
        ara_analysis_level = pd.read_csv("ara_analysis_level.csv")
        plt.plot(ara_analysis_level["energy"][-8:]/1e9,ara_analysis_level["limit"][-8:], color= "#2288AA")
        plt.annotate('ARA (2x1yr)',xy=(2.0e8, 1.6e-14), xycoords='data',horizontalalignment='left', color='#2288AA', rotation=-40, fontsize=12)
        #for x in range(1,6):
           #print(config["EFE_config{0}".format(x)])
           # plt.plot(energy_bins/1e9,config["EFE_config{0}".format(x)], label="config{0}".format(x))
            #plt.yscale('log')
        # plt.plot(energy_bins/1e9,combined_limit, label="combined limit", color='#2288AA')
        plt.plot(energy_bins/1e9,combined_limit, label="combined limit", color='black', linewidth=3)
        plt.legend(loc = 'best')
    #plotter.beautify_limit(ax_veff)
    #plt.tight_layout()
    #plt.show()
    #fig.savefig("final_limits.png",edgecolor='none',bbox_inches="tight") #save the figure
    figure.show(legend_size=10, save_name='current_limits.png',dpi=300)
    #plt.show()

    if(make_counts_plot):
        counts_koteramax_total=np.zeros(6)
        #counts_koteramax_total=np.array(counts_koteramax_total)
        for x in range(1,6):
            aeff_interpolator = splrep(np.log10(energy_bins), np.log10(config["aeff_steradian_config{0}".format(x)]))	
            
            eff_energy_bins = np.array([17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25])
            config["eff_config1"] = np.array([0.34, 0.42, 0.47, 0.56, 0.58, 0.65, 0.69])
            config["eff_config2"] = np.array([0.25,0.34,0.40,0.44,0.48,0.60,0.64])
            config["eff_config3"] = np.array([0.33,0.43,0.49,0.58,0.63,0.67,0.87])
            config["eff_config4"] = np.array([0.29,0.38,0.46,0.55,0.58,0.65,0.69])
            config["eff_config5"] = np.array([0.31,0.40,0.46,0.54,0.60,0.67,0.80])

            eff_interpolator = splrep(eff_energy_bins, config["eff_config{0}".format(x)])
	    
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
                temp_counts_koteramax = np.trapz(temp_koteramax*temp_aeff*livetime[x-1]*temp_eff,temp_energy)
                counts_koteramax.append(temp_counts_koteramax)
                energy_bins_counts.append(np.power(10.,bin))
	        
            fig = plt.figure(figsize=(2.*11,8.5))
            ax_veff = fig.add_subplot(1,2,1)
            ax_counts = fig.add_subplot(1,2,2)
            ax_veff.plot(energy_bins,config["veff_steradian_config{0}".format(x)],'-o',linewidth=2.0,color='blue',label=r'Config{0}'.format(x))
            print(counts_koteramax)
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
            fig.savefig("counts_config{0}.png".format(x),edgecolor='none',bbox_inches="tight") #save the figure
            counts_koteramax=np.array(counts_koteramax)
            counts_koteramax_total+=counts_koteramax
            print(counts_koteramax_total)
            #plt.show()


        fig = plt.figure(figsize=(11,8.5))
        ax_counts = fig.add_subplot(1,2,2)
        n_koteramax, bins, patches= ax_counts.hist(energy_bins_counts,
                                                   bins=np.power(10.,np.arange(15,22,0.5)),
                                                   weights=counts_koteramax_total,
                                                   label=r'Kotera Max: %.3f'%np.sum(counts_koteramax_total),
						   fill=False, 
						   stacked=True, 
						   histtype='step', 
						   edgecolor='blue',
						   linewidth=4)
        
        plotter.beautify_counts(ax_counts)
        ax_counts.set_ylabel('Events',size=17) #modify this axis title from plotter default
        fig.savefig("counts_total.png",edgecolor='none',bbox_inches="tight") #save the figure
"""
Uncomment the fuction, depending of what you want to get.
"""
#get_veff_plot()
plot_final_limit(make_counts_plot=True)

    
