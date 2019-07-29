import numpy as np #import numpy
import math
from scipy.interpolate import splrep, splev
import constants as const

def get_Lint(energy_eV, which_sigma=None):
	"""
	get_Lint
	This is a parameterization of the neutrino effective length.
	By default, it will use a cross section 
	measurement from Ghandi et. al. (10.1103/PhysRevD.58.093009)
	Parameters
	----------
	logev (float): energy in eV
		energy of your neutrino in electron volts
	Returns
	-------
	Lint (float): interaction length
		interaction length of a neutrino in centimeters
	"""
	Lint=0.

	if(which_sigma==None):
		#by default, assuem Ghandi et. al. (10.1103/PhysRevD.58.093009)
                sigma = 7.84e-36 * ((energy_eV/1.e9)**0.363)
                Lint = const.NucleonMass / (const.EarthDensity * sigma)
	if(which_sigma==1):
                print("Using Connolly et al. Xsec")
                sigma = nuCrsScnCTW(energy_eV)*1e4 #Converting to cm^2
                Lint = const.NucleonMass / (const.EarthDensity * sigma)
        
	return Lint

def get_flux(resource_name, energy_vals_logeV):
	"""
	get_flux
	Return the flux of neutrinos as a function of energy
	Parameters
	----------
	resource_name: name of the flux you want
		name of the flux you want
	energy_vals_logeV:
		the energies you want the flux evaluted at, in units of log10(eV)
	Returns
	-------
	flux: ndarray
		the flux prediction in units of 1/cm^2/s/sr
	"""
	flux = np.array([])

	if(resource_name=='kotera_max'):
		data = np.genfromtxt("data/kotera_max_PyREx.csv",delimiter=',',skip_header=1,names=['energy','flux'])
		data_energy_logev = data['energy']
		data_energy = np.power(10.,data_energy_logev)
		data_limit_log = data['flux']
		data_limit = np.power(10.,data_limit_log)
		data_interpolator = splrep(data_energy_logev, data_limit_log,k=4)

		energy_vals = np.power(10.,energy_vals_logeV)
		flux = np.power(10.,splev(energy_vals_logeV, data_interpolator))/energy_vals
		return flux

	return flux


def nuCrsScnCTW(energy):								    
    e = np.log10(energy*1e-9)
    cCC = (-1.826, -17.31, -6.406, 1.431, -17.91) # nu, cCC
    l = np.log(e - cCC[0])
    xCC = cCC[1] + (cCC[2] * l) + (cCC[3] * l*l) + (cCC[4]/l)
    sCC  = 10.0** xCC
    cNC = (-1.826, -17.31, -6.448, 1.431, -18.61) # nu, NC
    xNC = cNC[1] + (cNC[2] * l) + (cNC[3] * l*l) + (cNC[4]/l)
    sNC  = 10.0 ** xNC
    s = sCC+sNC
    s *= 1e-4 # cm^2 -> m^2
    return s

