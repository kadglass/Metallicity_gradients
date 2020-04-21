
### import python libraries ###

import matplotlib.pyplot as plt

from numpy import log10, pi

from getspectro import *

from marvin.tools.maps import Maps
from marvin import config
config.setDR("DR15") # config.download = True

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import astropy.constants as const



### read maps ###
plate = "7443"
fiber = "12705"

ID = plate + "-" + fiber

maps = Maps(ID)  # download = True

    # observed

HaF_map = maps["emline_gflux_ha_6564"]
HbF_map = maps["emline_gflux_hb_4862"]
OII_map = maps["emline_gflux_oii_3727"]
OIII_map = maps["emline_gflux_oiii_5008"]
NII_map = maps["emline_gflux_nii_6585"]

    # calculated

N2_map = NII_map / HaF_map
O3N2_map = OIII_map / HbF_map / N2_map
N2O2_map = NII_map / OII_map

    # masks

nocov_HaF = HaF_map.pixmask.get_mask("NOCOV")
nocov_HbF = HbF_map.pixmask.get_mask("NOCOV")
mask = nocov_HaF | nocov_HbF




## calculate metallicity ###

#constants

c = const.c.to("km/s")  # [speed of light] = [astropy.constants].[c].to[desired units]
H0 = 70 * (u.km/u.s/u.Mpc) # [astropy.units].[desired units]

    # retrieve redshift & mass

drpall = fits.open("drpall-v2_4_3.fits")  #drpall[1].header
drpall_table = Table()  # Format: [variable name][z].[y][x]
drpall_table["plateifu"] = drpall[1].data["plateifu"]
drpall_table["z"] = drpall[1].data["z"]
drpall.close()

index = drpall_table["plateifu"] == ID
z = drpall_table["z"][index]


mass = getspectro(plate, fiber)

# equations from moorman:
#   fHa_map = HaF_map * ( 10 ** ( AHa_map / 2.5) ) = corrected H Alpha
#   LHa_map = 4 * pi * ( Dl ** 2 ) * fHa_map = H Alpha luminosity
#   SFR_map = 7.9 * ( 10 ** -41.28 ) * ( LHa_map ) = star formation rate: H Alpha Method

Dl = (( z * c ) / H0).to("cm")  # Distance Luminosity

AHa_map = 5.91 * log10( HaF_map.value / HbF_map.value ) -2.70  # Dust Attentuation Correction for H Alpha

logSFR_map = -41.28 + log10(7.9) + log10(4*pi*Dl.value*Dl.value) + log10(HaF_map.value) + (AHa_map / 2.5)  # Change to logs and combine


    # retrieve mass

logsSFR_map = log10(logSFR_map) - log10(mass)
aveSSFR = 283.728 - ( 116.265 * mass) + ( 17.4403 * (mass ** 2) ) - ( 1.17146 * (mass ** 3) ) + ( 0.0296526 * (mass ** 4) )
d_logSSFR_map = logsSFR_map - aveSSFR

    # N2 Method
N2_metallicity = 9.12  + ( 0.58 * log10(N2_map) ) - ( 0.19 * d_logSSFR_map )
N2_metallicity.plot()

    # O3N2 Method
O3N2_metallicity = 8.98 - ( 0.32 * log10(O3N2_map) ) - ( 0.18 * d_logSSFR_map )
O3N2_metallicity.plot()

    # N2O2 Method
N2O2_metallicity = 9.20 + ( 0.54 * log10(N2O2_map) ) - ( 0.36 * d_logSSFR_map )
N2O2_metallicity.plot()

plt.show()

#plt.savefig(f"/Users/leilani/Desktop/SHELS/marvin/plots/SFR_map.pdf", overwrite = True)





