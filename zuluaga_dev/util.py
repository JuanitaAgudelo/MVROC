import textwrap
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

def print_list(lst, width=80):
    print(textwrap.fill(", ".join(lst), width=width))

DEG = np.pi / 180.0
RAD = 1.0 / DEG

# Convert longitudes [0, 360) to Mollweide angles: 0° at center, East (0° -> 180°) left, West (360° -> 180°) right
f2m = lambda lon: -((lon + 180.0) % 360.0 - 180.0)
def split_mollweide(ra_deg, dec_deg):
    """Convert RA and Dec to Mollweide radians and split lines across the 180° boundary."""
    x = f2m(ra_deg) * DEG
    y = np.deg2rad(dec_deg)
    jumps = np.where(np.abs(np.diff(x)) > np.pi)[0]
    return np.insert(x, jumps + 1, np.nan), np.insert(y, jumps + 1, np.nan)

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
