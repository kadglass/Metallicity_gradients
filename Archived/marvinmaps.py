
### Import Libraries ###

import matplotlib.pyplot as plt

from numpy import log10, pi

from marvin.tools.maps import Maps
from marvin import config
config.setDR("DR15") # config.download = True

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import astropy.constants as const


### Read Maps ###

ID = "7443-12705"
maps = Maps(ID)  # download = True

HaF_map = maps["emline_gflux_ha_6564"]
HbF_map = maps["emline_gflux_hb_4862"]
OII_map = maps["emline_gflux_oii_3727"]
OIII_map = maps["emline_gflux_oiii_5008"]
NII_map = maps["emline_gflux_nii_6585"]


### Plot Emission Lines ###

HaF_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/HaF_map.pdf", overwrite = True)
HbF_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/HbF_map.pdf", overwrite = True)
OII_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/OII_map.pdf", overwrite = True)
OIII_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/OIII_map.pdf", overwrite = True)
NII_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/NII_map.pdf", overwrite = True)


### Plot Emission Ratios ###

N2_map = NII_map / HaF_map
O3N2_map = OIII_map / HbF_map / N2_map
N2O2_map = NII_map / OII_map

N2_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/N2_map.pdf", overwrite = True)
O3N2_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/O3N2_map.pdf", overwrite = True)
N2O2_map.plot()
plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/N2O2_map.pdf", overwrite = True)


### Calculate Metallicity ###

c = const.c.to("km/s")  # [speed of light] = [astropy.constants].[c].to[desired units]
H0 = 70 * (u.km/u.s/u.Mpc) # [astropy.units].[desired units]

# drpall has useful information about each of the galaxies in marvin
# drpall is a fits file
#   Format: [variable name][z].[y][x]
drpall = fits.open("drpall-v2_4_3.fits")  #drpall[1].header
drpall_table = Table()
drpall_table["plateifu"] = drpall[1].data["plateifu"]
drpall_table["z"] = drpall[1].data["z"]
drpall.close()

# given galaxy ID, find its redshift
index = drpall_table["plateifu"] == ID
z = drpall_table["z"][index]



### Nifty Equations ###

Dl = ( z * c ) / H0  # Distance Luminosity
Dl = Dl.to("cm")

AHa_map = 5.91 * log10( HaF_map / HbF_map ) -2.70  # Dust Attentuation Correction for H 

# Reference Equations:
#   fHa_map = HaF_map * ( 10 ** ( AHa_map / 2.5) )  # Corrected H Alpha
#   LHa_map = 4 * pi * ( Dl ** 2 ) * fHa_map  # H Alpha Luminosity
#   SFR_map = 7.9 * ( 10 ** -41.28 ) * ( LHa_map )  # Star Formation Rate: H Alpha Method

# marvin doesn't want to take a number to the power of a map
# logfHa_map = log10(HaF_map) + (AHa_map.value / 2.5)
# logLHa_map = log10(4*pi*Dl.value*Dl.value) + logfHa_map
# logSFR_map = -41.28 + log10(7.9) + logLHa_map

logSFR_map = -41.28 + log10(7.9) + log10(4*pi*Dl.value*Dl.value) + log10(HaF_map.value) + (AHa_map.value / 2.5)

logSFR_map.plot()
plt.show()
#plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/SFR_map.pdf", overwrite = True)





