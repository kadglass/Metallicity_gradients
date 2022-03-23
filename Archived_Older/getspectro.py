import os

from astropy.io import fits
from astropy.table import table

def getspectro(plate, fiber):

    # des
    file = "manga-{0}-{1}.Pipe3D.cube.fits.gz".format(plate, fiber)

    os.system("rsync -avz --no-motd rsync://data.sdss.org/dr16/manga/spectro/pipe3d/v2_4_3/2.4.3/{0}/{1} /Users/leilani/Desktop/Metallicity_gradients".format(plate, file))

    cube = fits.open(file)

    index = cube[1].header.index("FILE_18")

    if cube[1].header["FILE_18"] == "map.CS.manga-{0}-{1}_Mass_ssp.fits.gz".format(plate, fiber):
        pass

    else:
        print(":(")


    mass = cube[1].data[18]

    return(mass)