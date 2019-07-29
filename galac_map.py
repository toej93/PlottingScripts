import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.coordinates as coord
from astropy.io import ascii

data = ascii.read("analysis_window_clean_altaz_FOV_2018.04.11.csv")
	
ra = coord.Angle(data['RA']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(data['Decl']*u.degree)
#ax.set_yticklabels([])
#ax.set_xticklabels([])

	    
	#also want to plot the galactic plane
galactic_longitudes = np.arange(start=0, stop=360, step=0.1)
galactic_latitudes = [0] * len(galactic_longitudes)
icrs = SkyCoord(galactic_longitudes, galactic_latitudes, unit="deg", frame="galactic").icrs
gal_ra = icrs.ra
gal_ra = gal_ra.wrap_at(180*u.degree)
gal_dec = icrs.dec
	
gal_center_long = 180.
gal_center_lat = 0.
cena_long = 309.51589568
cena_lat = 19.41727350
icrs2 = SkyCoord(gal_center_long, gal_center_lat, unit="deg", frame="galactic").icrs
gal_center_ra = icrs2.ra
gal_center_ra = gal_center_ra.wrap_at(180*u.degree)
gal_center_dec = icrs2.dec

icrs3 = SkyCoord(cena_long, cena_lat, unit="deg", frame="galactic").icrs
#cena_ra = icrs3.ra
#cena_ra = cena_ra.wrap_at(180*u.degree)
#cena_dec = icrs3.dec
ra_gc=299.3*u.degree
dec_gc=-28.72* u.degree
ra_cena = 201.3625*u.degree
dec_cena = -43.0192*u.degree
c2 = SkyCoord(ra=ra_cena, dec=dec_cena, frame='icrs')
cena_ra = c2.ra.wrap_at(180 * u.deg).radian
cena_dec = c2.dec.radian


c = SkyCoord(ra=ra_gc, dec=dec_gc, frame='icrs')
fig = plt.figure()
ax = fig.add_subplot(111, projection="hammer")
ra_rad = c.ra.wrap_at(180 * u.deg).radian
dec_rad = c.dec.radian
ax.grid(color='k', linestyle='solid', linewidth=0.5)
r = 90
#theta = np.arange(0,2*np.pi,0.1)
x = np.array([-np.pi,np.pi,np.pi,-np.pi,-np.pi])
y = np.array([0.05,0.05,-0.36,-0.36,0.05])
y2 = np.array([0.05,0.05,-0.52,-0.52,0.05])
y3 = np.array([-0.52,-0.52,-0.64,-0.64,-0.52])
#ax.plot(x,y)

ax.fill(x,y2,label='60 m depth', facecolor='#8b9dc3',alpha=1)
ax.fill(x,y3,label='Added reach at 100 m depth',alpha=1,facecolor='#3b5998')

#ax.fill(x, y,label='15 m', facecolor='#d0e1f9',alpha=1 )


ax.scatter(ra_rad,dec_rad)
ax.plot(ra_rad, dec_rad, '*',color='gold',markersize=11,mec='firebrick',label='Galactic Center')
#ax.plot(cena_ra, cena_dec,'^',markersize=9,color='m',label='Centaurus A')
ax.plot(gal_ra.radian[0:2970], -gal_dec.radian[0:2970],color='firebrick',linewidth=2,label='Galactic Plane',zorder=1)
ax.plot(gal_ra.radian[2980:], -gal_dec.radian[2980:],color='firebrick',linewidth=2,zorder=2)
#ax.plot(x+np.pi/6,y+np.pi/6)
#ax.plot(x-np.pi/2,y-np.pi/2);
ax.plot([-np.pi, np.pi], [-0.36, -0.36], '--',color='lightgreen', alpha=0.7)#15 m depth

legend = ax.legend(loc='upper right')
legend.get_frame().set_facecolor('#ffe4c4')
ax.legend( bbox_to_anchor=(.66,.8))
ax.set_ylabel('Declination (deg)') #give it a title
ax.set_xlabel('Right Ascension (deg)',labelpad=20) #give it a title
ax.axes.get_xaxis().set_ticks([-np.pi/3, -2*np.pi/3, -np.pi,0,np.pi/3, 2*np.pi/3, np.pi])
ax.annotate('15 m depth',xy=(np.pi, -0.36), xycoords='data',xytext=(.87, .3), textcoords='figure fraction',arrowprops=dict(arrowstyle="->", color='royalblue'),horizontalalignment='right', verticalalignment='top')
#ax.annotate('annotate', xy=(3, 0), xytext=(11, 4))
#plt.gca().set_aspect('0.8', adjustable='box')
plt.tight_layout()
fig.savefig("Sky_map.pdf")
