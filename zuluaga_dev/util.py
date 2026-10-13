import textwrap
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

def print_list(lst, width=80):
    print(textwrap.fill(", ".join(lst), width=width))

def ecl_gal(npoints=1000):
    # Ecliptic
    lon_ecl = np.linspace(0, 360, npoints) * u.deg
    ecl = SkyCoord(lon=lon_ecl, lat=0 * u.deg, frame='geocentrictrueecliptic').icrs
    idx_ecl = np.argsort(ecl.ra.deg)
    ecl_ra = ecl.ra.deg[idx_ecl]
    ecl_dec = ecl.dec.deg[idx_ecl]
    # Galactic 
    l_gal = np.linspace(0, 360, npoints) * u.deg
    gal = SkyCoord(l=l_gal, b=0 * u.deg, frame='galactic').icrs
    idx_gal = np.argsort(gal.ra.deg)
    gal_ra = gal.ra.deg[idx_gal]
    gal_dec = gal.dec.deg[idx_gal]

    return [[ecl_ra, ecl_dec], [gal_ra, gal_dec]]