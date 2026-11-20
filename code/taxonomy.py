"""Reflectance, gri slope, i-z and taxonomy classification of asteroid colours.

Method of Carruba et al. (2026): each object is placed on the plane
(gri slope [%/100nm], i-z [mag]) and assigned the class whose box contains it.

All slope functions accept scalars, numpy arrays or pandas Series, so they work
on one object or on a whole catalogue.

Conventions
-----------
- Colours are asteroid colours in magnitudes: g-r, r-i, i-z (g-i = g-r + r-i).
- Reflectance normalised to g, for a band f:
      R_fg = 10^(-0.4 [(m_f - m_g) - (M_f_sun - M_g_sun)])
  so R_gg = 1. Note that m_f - m_g = -(g - f), e.g. m_r - m_g = -(g-r).
- gri slope (spectral gradient between r and i), as used in the paper:
      S = 1e4 (R_rg - R_ig) / |lam_i - lam_r|     [%/100nm]
- The paper's slopes are reproduced with gri_slope_carruba (M_SUN_SDSS, SDSS wavelengths).
"""
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
# Effective wavelengths (nm)
LAM_LSST = {'u': 372.4, 'g': 480.7, 'r': 622.1, 'i': 755.9, 'z': 868.0, 'y': 975.3}
LAM_SDSS = {'g': 467.3, 'r': 614.2, 'i': 745.9, 'z': 892.5}

# Solar absolute magnitudes (solar colours g-r 0.70, r-i 0.34, i-z 0.18)
M_SUN_SDSS = {'g': -26.34, 'r': -27.04, 'i': -27.38, 'z': -27.56}  #Vega
M_SUN_LSST = {'g': -26.52, 'r': -26.93, 'i': -27.05, 'z': -27.07}  #AB

# Solar colours (g-r, r-i, i-z) in mag
SUN_COLORS = {
    'notebook': {'g-r': M_SUN_LSST['g'] - M_SUN_LSST['r'],
                 'r-i': M_SUN_LSST['r'] - M_SUN_LSST['i'],
                 'i-z': M_SUN_LSST['i'] - M_SUN_LSST['z']},
    # catalogue pipeline (util.py / catalogue_pipeline.py); does NOT reproduce the paper
    'pipeline': {'g-r': 0.45, 'r-i': 0.12, 'i-z': 0.04},
}

_K = 0.4 * np.log(10)  # d(10^(-0.4 x))/dx = -_K * 10^(-0.4 x)


def _sun(sun):
    """Accept a preset name or a dict with 'g-r', 'r-i', 'i-z'."""
    return SUN_COLORS[sun] if isinstance(sun, str) else sun


# --------------------------------------------------------------------------
# Reflectance and slopes
# --------------------------------------------------------------------------
def reflectance(Mfg, Mf_sun, Mg_sun):
    """Reflectance of band f normalised to g.

    Mfg    : m_f - m_g of the asteroid (e.g. -(g-r) for f = r, -(g-i) for f = i)
    Mf_sun : solar absolute magnitude in band f (e.g. M_SUN_LSST['r'])
    Mg_sun : solar absolute magnitude in g, from the same set as Mf_sun
    """
    Mfg = np.asarray(Mfg, dtype=float)
    dM_sun = Mf_sun - Mg_sun              # (M_f - M_g) of the Sun
    return 10 ** (-0.4 * (Mfg - dM_sun))


def spectral_gradient_carruba(Rrg, Rig, lamr, lami):
    """S = 1e4 (Rrg - Rig) / |lamr - lami|, in %/100nm (R normalised to 1 at g)."""
    return 1e4 * (np.asarray(Rrg) - np.asarray(Rig)) / abs(lamr - lami)


def gri_slope_carruba(g_r, r_i, lam=LAM_SDSS, m_sun=M_SUN_SDSS):
    """gri slope [%/100nm] from the colours g-r and r-i.

    lam = LAM_SDSS reproduces the paper; LAM_LSST uses the Rubin filters.
    m_sun : solar absolute magnitudes (default M_SUN_SDSS; e.g. M_SUN_LSST).
    """

    g_i = np.asarray(g_r, dtype=float) + np.asarray(r_i, dtype=float)
    Rrg = reflectance(-np.asarray(g_r, dtype=float), m_sun['r'], m_sun['g'])
    Rig = reflectance(-g_i, m_sun['i'], m_sun['g'])
    return spectral_gradient_carruba(Rrg, Rig, lam['r'], lam['i'])


def gri_slope_err_carruba(g_r, r_i, g_r_err, r_i_err, lam=LAM_SDSS, m_sun=M_SUN_SDSS):
    """1-sigma error of gri_slope by linear error propagation.

    With R_r = R_r(g-r) and R_i = R_i(g-r + r-i), and g-r, r-i independent:
        dS/d(g-r) = k (R_r - R_i),   dS/d(r-i) = -k R_i,   k = 1e4 * 0.4 ln10 / dlam
    """
    g_r, r_i = np.asarray(g_r, dtype=float), np.asarray(r_i, dtype=float)
    Rrg = reflectance(-g_r, m_sun['r'], m_sun['g'])
    Rig = reflectance(-(g_r + r_i), m_sun['i'], m_sun['g'])
    k = 1e4 * _K / abs(lam['r'] - lam['i'])
    return k * np.sqrt(((Rrg - Rig) * np.asarray(g_r_err)) ** 2 + (Rig * np.asarray(r_i_err)) ** 2)


def _fit_weights(lam):
    """Least-squares weights w_k so that the fitted slope of R(lam_g, lam_r, lam_i) is sum(w_k R_k)."""
    x = np.array([lam['g'], lam['r'], lam['i']], dtype=float)
    d = x - x.mean()
    return d / (d ** 2).sum()


def gri_slope_fit_colors(g_r, r_i, lam=LAM_LSST, m_sun=M_SUN_LSST):
    """gri slope [%/100nm] from an unweighted linear fit of the reflectance vs wavelength in g, r and i,
    computed from the colours (closed form, no errors used; see gri_slope_fit for the weighted fit).

    The reflectances are normalised to g (R_g = 1, R_r, R_i from reflectance()); the slope of the
    least-squares line R(lam) through the three points is returned in %/100nm (1e4 * dR/dlam[nm]).
    Positive for objects redder than the Sun.
    """

    g_r = np.asarray(g_r, dtype=float)
    g_i = g_r + np.asarray(r_i, dtype=float)
    Rr = reflectance(-g_r, m_sun['r'], m_sun['g'])
    Ri = reflectance(-g_i, m_sun['i'], m_sun['g'])
    w = _fit_weights(lam)
    return 1e4 * (w[0] * 1.0 + w[1] * Rr + w[2] * Ri)


def gri_slope_fit_colors_err(g_r, r_i, g_r_err, r_i_err, lam=LAM_LSST, m_sun=M_SUN_LSST):
    """1-sigma error of gri_slope_fit_colors by linear error propagation.

    dR_r/d(g-r) = k R_r, dR_i/d(g-r) = dR_i/d(r-i) = k R_i with k = 0.4 ln10, so
        dS/d(g-r) = 1e4 k (w_r R_r + w_i R_i),   dS/d(r-i) = 1e4 k w_i R_i
    """
    g_r = np.asarray(g_r, dtype=float)
    g_i = g_r + np.asarray(r_i, dtype=float)
    Rr = reflectance(-g_r, m_sun['r'], m_sun['g'])
    Ri = reflectance(-g_i, m_sun['i'], m_sun['g'])
    w = _fit_weights(lam)
    d_gr = 1e4 * _K * (w[1] * Rr + w[2] * Ri)
    d_ri = 1e4 * _K * w[2] * Ri
    return np.sqrt((d_gr * np.asarray(g_r_err)) ** 2 + (d_ri * np.asarray(r_i_err)) ** 2)


def reflectance_err(R, m_err):
    """1-sigma error of a reflectance R = 10^(-0.4 m) by error propagation: sigma_R = 0.4 ln10 R sigma_m."""
    return _K * np.abs(np.asarray(R, dtype=float)) * np.asarray(m_err, dtype=float)


def gri_reflectances(m_g, m_r, m_i, m_g_err, m_r_err, m_i_err, m_sun=M_SUN_LSST):
    """Reflectances in g, r, i and their 1-sigma errors from the magnitude of each band (e.g. H per band).

        R_f = 10^(-0.4 [(m_f - m_g) - (M_f_sun - M_g_sun)])      (R_g = 1)

    The measured m_g only sets the normalisation (a constant), so every band keeps its own error:
        sigma_R_f = 0.4 ln10 R_f sigma_m_f
    and the three points are independent (sigma_R_g = 0.4 ln10 sigma_m_g).
    Returns a DataFrame with R_g, R_r, R_i, R_g_err, R_r_err, R_i_err.
    """
    index = getattr(m_g, 'index', None)
    m = {'g': np.atleast_1d(np.asarray(m_g, dtype=float)), 'r': np.atleast_1d(np.asarray(m_r, dtype=float)),
         'i': np.atleast_1d(np.asarray(m_i, dtype=float))}
    e = {'g': np.atleast_1d(np.asarray(m_g_err, dtype=float)), 'r': np.atleast_1d(np.asarray(m_r_err, dtype=float)),
         'i': np.atleast_1d(np.asarray(m_i_err, dtype=float))}
    out = {}
    for b in 'gri':
        out[f'R_{b}'] = reflectance(m[b] - m['g'], m_sun[b], m_sun['g'])
    for b in 'gri':
        out[f'R_{b}_err'] = reflectance_err(out[f'R_{b}'], e[b])
    return pd.DataFrame(out, index=index)


def _line(x, a, b):
    return a + b * x


def gri_slope_fit(R, R_err, lam=LAM_LSST, bands=('g', 'r', 'i')):
    """gri slope from a weighted linear fit R(lam) = a + b lam of the reflectances, with scipy curve_fit.

    R, R_err : arrays of shape (n, 3) or (3,) (or DataFrames) with the reflectances and their 1-sigma
               errors, columns in the order of `bands`.
    The wavelengths are taken as exact (no error); sigma = R_err with absolute_sigma=True, so the errors
    of a and b come only from the reflectance errors (covariance matrix of the fit).

    Returns a DataFrame, one row per object:
        slope, slope_err          1e4 b [%/100nm]  (R normalised to g)
        intercept, intercept_err  a (reflectance at lam = 0)
        cov_ab                    covariance of a and b (b in %/100nm)
        chi2, ndof                sum(((R - a - b lam) / R_err)^2), n_points - 2
        p_value                   P(chi2_ndof >= chi2): probability of a chi2 this large if the line is right
                                  and the errors are correct (small p = poor fit; p < 0.05 <=> chi2 > 3.84)
    Rows with a missing value or a non-positive error get NaN.
    """
    from scipy.optimize import curve_fit
    from scipy.stats import chi2 as chi2_dist

    x = np.array([lam[b] for b in bands], dtype=float)
    R = np.atleast_2d(np.asarray(R, dtype=float))
    E = np.atleast_2d(np.asarray(R_err, dtype=float))
    nan = dict(slope=np.nan, slope_err=np.nan, intercept=np.nan, intercept_err=np.nan,
               cov_ab=np.nan, chi2=np.nan, ndof=len(x) - 2, p_value=np.nan)
    rows = []
    for y, e in zip(R, E):
        if not (np.all(np.isfinite(y)) and np.all(np.isfinite(e)) and np.all(e > 0)):
            rows.append(nan)
            continue
        p0 = np.polyfit(x, y, 1)[::-1]
        popt, pcov = curve_fit(_line, x, y, p0=p0, sigma=e, absolute_sigma=True)
        chi2 = float(np.sum(((y - _line(x, *popt)) / e) ** 2))
        ndof = len(x) - 2
        rows.append(dict(slope=1e4 * popt[1], slope_err=1e4 * np.sqrt(pcov[1, 1]),
                         intercept=popt[0], intercept_err=np.sqrt(pcov[0, 0]),
                         cov_ab=1e4 * pcov[0, 1], chi2=chi2, ndof=ndof,
                         p_value=float(chi2_dist.sf(chi2, ndof))))
    return pd.DataFrame(rows)


def iz_color(i_z, sun='notebook', solar_corrected=False):
    """i-z for the diagram.

    solar_corrected=False returns the raw colour (what the paper plots).
    solar_corrected=True returns (i-z) - (i-z)_sun = 2.5 log10(R_z / R_i).
    """
    i_z = np.asarray(i_z, dtype=float)
    return i_z - _sun(sun)['i-z'] if solar_corrected else i_z


# --------------------------------------------------------------------------
# Apply to a catalogue
# --------------------------------------------------------------------------
def taxonomy_params(df, g_r='g-r', r_i='r-i', i_z='i-z',
                    g_r_err='g-r_err', r_i_err='r-i_err', i_z_err='i-z_err'):
    """Slopes (SDSS and LSST wavelengths), their errors and i-z for every row of df.

    Arguments are column names in df; error columns may be None.
    """
    out = pd.DataFrame(index=df.index)
    for name, lam in [('sdss', LAM_SDSS), ('lsst', LAM_LSST)]:
        out[f'slope_{name}'] = gri_slope_carruba(df[g_r], df[r_i], lam=lam)
        if g_r_err and r_i_err:
            out[f'slope_{name}_err'] = gri_slope_err_carruba(df[g_r], df[r_i], df[g_r_err],
                                                             df[r_i_err], lam=lam)
    out['i_z'] = iz_color(df[i_z])
    if i_z_err:
        out['i_z_err'] = df[i_z_err].astype(float)
    return out


def add_h_colors(df):
    """Add colours from the absolute magnitudes: gr_H, ri_H, iz_H."""
    df = df.copy()
    df['gr_H'] = df['H_g'] - df['H_r']
    df['ri_H'] = df['H_r'] - df['H_i']
    df['iz_H'] = df['H_i'] - df['H_z']
    return df


# --------------------------------------------------------------------------
# Taxonomy boxes and classification
# --------------------------------------------------------------------------
# Class boxes in the (gri slope [%/100nm], i-z [mag]) plane: (x_min, x_max, y_min, y_max, colour)
# colour of each class (boxes, regions and class colours in the plots)
PAL = {'C': '#0072B2', 'B': '#E6AB02', 'X': '#009E73', 'D': '#8C564B', 'L': '#CC79A7',
       'S': '#D55E00', 'Q': '#999933', 'A': '#117733', 'V': '#882255'}
TAXO_BOXES = {
    'A': (21.5, 28.0, -0.265, -0.115, PAL['A']),
    'B': (-5.0,  0.0, -0.200,  0.000, PAL['B']),
    'C': (-5.0,  6.0, -0.200,  0.180, PAL['C']),
    'D': ( 6.0, 25.0,  0.080,  0.330, PAL['D']),
    'L': ( 9.0, 25.0,  0.000,  0.085, PAL['L']),
    'Q': ( 5.0,  9.5, -0.265, -0.165, PAL['Q']),
    'S': ( 6.0, 25.0, -0.265,  0.000, PAL['S']),
    'V': ( 5.0, 25.0, -0.665, -0.265, PAL['V']),
    'X': ( 2.5,  9.0,  0.000,  0.180, PAL['X']),
}
LABEL_POS = {'A': (21.8, -0.17), 'B': (-4.5, -0.06), 'C': (-4.5, 0.12),
             'D': (6.5, 0.27), 'L': (10.0, 0.01), 'Q': (8.0, -0.24),
             'S': (6.5, -0.06), 'V': (5.5, -0.33), 'X': (3.0, 0.12)}

# Probabilistic classification parameters
SYS_SLOPE = 0.25   # %/100nm, systematic floor on the slope (RMS of our slopes vs the paper)
SYS_IZ = 0.01      # mag, systematic floor on i-z (precision of the published colours)
P_MIN = 0.60       # minimum probability of the most likely class to assign it
N_MC = 20000       # random draws per object

INDETERMINATE = 'indeterminate'
UNCLASSIFIED = 'unclassified'   # most of the probability falls outside every box
NO_DATA = 'no data'             # slope or i-z missing


# Classification regions. Several boxes overlap; each overlap belongs to one class, so every
# point of the plane has exactly one class (or none). The first box in this list that contains
# the point wins:
#   B over C (B lies inside C)        X over C and D (C/X and D/X overlaps are X)
#   Q over S and C                    L over D (thin D/L strip at i-z 0.080-0.085)
#   A over S
REGION_PRIORITY = ['B', 'Q', 'A', 'X', 'L', 'C', 'D', 'S', 'V']


def _contains(box, x, y):
    x0, x1, y0, y1 = box[:4]
    return x0 <= x <= x1 and y0 <= y <= y1


def region_index(x, y):
    """Index into REGION_PRIORITY of the region containing each point (-1 = outside every box).

    x, y: arrays (or scalars) of gri slope and i-z.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    out = np.full(np.broadcast(x, y).shape, -1, dtype=int)
    for k in range(len(REGION_PRIORITY) - 1, -1, -1):   # lowest priority first, overwritten by higher
        x0, x1, y0, y1 = TAXO_BOXES[REGION_PRIORITY[k]][:4]
        out[(x0 <= x) & (x <= x1) & (y0 <= y) & (y <= y1)] = k
    return out


def classify(slope, iz):
    """Class of the region that contains the point ('' if none). No errors."""
    if not (np.isfinite(slope) and np.isfinite(iz)):
        return ''
    k = int(region_index(slope, iz))
    return REGION_PRIORITY[k] if k >= 0 else ''


def class_probabilities(slope, iz, slope_err=0.0, iz_err=0.0, sys_slope=SYS_SLOPE,
                        sys_iz=SYS_IZ, n_mc=N_MC, seed=0):
    """Probability that the object belongs to each class.

    The true position is drawn n_mc times from a Gaussian centred on (slope, iz) with
        sigma_slope = sqrt(slope_err^2 + sys_slope^2),  sigma_iz = sqrt(iz_err^2 + sys_iz^2)
    and each draw votes for the class of the region that contains it (see REGION_PRIORITY:
    overlaps of boxes belong to one class), or for 'none' if it is outside every box.

    Returns (dict {class: probability}, probability of 'none'). The values add up to 1.
    """
    rng = np.random.default_rng(seed)
    se = 0.0 if not np.isfinite(slope_err) else float(slope_err)
    ze = 0.0 if not np.isfinite(iz_err) else float(iz_err)
    x = slope + np.hypot(se, sys_slope) * rng.standard_normal(n_mc)
    y = iz + np.hypot(ze, sys_iz) * rng.standard_normal(n_mc)

    k = region_index(x, y)
    counts = np.bincount(k[k >= 0], minlength=len(REGION_PRIORITY)) / n_mc
    P = {c: float(counts[REGION_PRIORITY.index(c)]) for c in TAXO_BOXES}
    return P, float((k < 0).mean())


def classify_prob(slope, iz, slope_err=0.0, iz_err=0.0, p_min=P_MIN, **kw):
    """Classify one object from its class probabilities.

        most likely class has P >= p_min       -> that class
        P(outside every box) is the largest
          and >= p_min                         -> 'unclassified'
        otherwise                              -> 'indeterminate'

    Returns (taxonomy, p_class, candidates, probabilities):
      p_class    : probability of the most likely class
      candidates : classes with P >= 0.05, most likely first, e.g. 'L 0.55 / S 0.40'
    """
    if not (np.isfinite(slope) and np.isfinite(iz)):
        return NO_DATA, np.nan, '', {}
    P, p_none = class_probabilities(slope, iz, slope_err, iz_err, **kw)
    best = max(P, key=P.get)
    ranked = sorted(((p, c) for c, p in P.items() if p >= 0.05), reverse=True)
    candidates = ' / '.join(f'{c} {p:.2f}' for p, c in ranked)
    if p_none > P[best]:
        tax = UNCLASSIFIED if p_none >= p_min else INDETERMINATE
        return tax, p_none, candidates, P
    tax = best if P[best] >= p_min else INDETERMINATE
    return tax, P[best], candidates, P


def classify_catalogue(slope, iz, slope_err=None, iz_err=None, **kw):
    """classify_prob on every row.

    Returns a DataFrame with taxonomy, p_class, candidates and one P_<class> column per class.
    """
    slope = pd.Series(slope)
    idx = slope.index
    n = len(idx)
    iz = np.asarray(iz, dtype=float)
    se = np.zeros(n) if slope_err is None else np.asarray(slope_err, dtype=float)
    ze = np.zeros(n) if iz_err is None else np.asarray(iz_err, dtype=float)
    rows = []
    for s, z, a, b in zip(slope, iz, se, ze):
        tax, p, cand, P = classify_prob(s, z, a, b, **kw)
        rows.append({'taxonomy': tax, 'p_class': p, 'candidates': cand,
                     **{f'P_{c}': P.get(c, np.nan) for c in TAXO_BOXES}})
    return pd.DataFrame(rows, index=idx)


# --------------------------------------------------------------------------
# Taxonomy diagram
# --------------------------------------------------------------------------
def plot_taxonomy(points, title=None, ax=None, legend=True, auto_limits=True,
                  fill_regions=True, region_alpha=0.22):
    """Plot objects on the gri-slope vs i-z diagram.

    fill_regions shades each class region (as used by the classification) with its colour.

    points: list of dicts with keys
        label, slope, iz, color, and optional slope_err, iz_err,
        name_text (text box next to the star), text_offset (points), marker, size
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5.2))

    if fill_regions:   # shade the classification regions (overlaps take the colour of their class)
        from matplotlib.colors import to_rgba
        xs = np.linspace(-5.0, 28.0, 1321)
        ys = np.linspace(-0.665, 0.330, 996)
        k = region_index(xs[None, :], ys[:, None])
        rgba = np.zeros(k.shape + (4,))
        for i, cls in enumerate(REGION_PRIORITY):
            rgba[k == i] = to_rgba(TAXO_BOXES[cls][4], alpha=region_alpha)
        ax.imshow(rgba, extent=(xs[0], xs[-1], ys[0], ys[-1]), origin='lower', aspect='auto',
                  interpolation='nearest', zorder=0)

    for cls, (x0, x1, y0, y1, c) in TAXO_BOXES.items():
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=c, lw=2.5))
        lx, ly = LABEL_POS[cls]
        ax.text(lx, ly, cls, color=c, fontsize=15, fontweight='bold')

    for p in points:
        se, ye = p.get('slope_err', 0) or 0, p.get('iz_err', 0) or 0
        se = 0 if not np.isfinite(se) else se
        ye = 0 if not np.isfinite(ye) else ye
        if se or ye:  # 1-sigma error bars with short caps
            ax.errorbar(p['slope'], p['iz'], xerr=se, yerr=ye, fmt='none', ecolor=p['color'],
                        elinewidth=1.8, capsize=5, capthick=1.8, zorder=4)
        ax.scatter(p['slope'], p['iz'], marker=p.get('marker', '*'), s=p.get('size', 350),
                   c=p['color'], edgecolors='white', linewidths=1, zorder=5, label=p.get('label'))
        if p.get('name_text'):
            ax.annotate(p['name_text'], (p['slope'], p['iz']),
                        xytext=p.get('text_offset', (10, 10)), textcoords='offset points',
                        fontsize=13, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.6'))

    # default limits, widened if a point falls outside them
    xs = [p['slope'] for p in points] if auto_limits else []
    ys = [p['iz'] for p in points] if auto_limits else []
    ax.set_xlim(min([-10] + [x - 3 for x in xs]), max([30] + [x + 3 for x in xs]))
    ax.set_ylim(min([-0.8] + [y - 0.1 for y in ys]), max([0.4] + [y + 0.1 for y in ys]))
    ax.set_xlabel('gri slope [%/100nm]', fontsize=15)
    ax.set_ylabel('i-z [mag]', fontsize=15)
    ax.tick_params(labelsize=13)
    ax.minorticks_on()
    ax.grid(True, which='major', alpha=0.4)
    ax.grid(True, which='minor', alpha=0.15)
    if title:
        ax.set_title(title, fontsize=17)
    if legend:
        ax.legend(loc='lower left', fontsize=12, framealpha=0.9)
    return ax


def plot_object(row, name, out_dir=None, dpi=120):
    """Taxonomy diagram for one object of the classified table.

    row needs: slope_sdss, slope_sdss_err, slope_lsst, slope_lsst_err, i_z, i_z_err,
    taxonomy, p_class, candidates. Saves <out_dir>/<name with '_'>_taxonomy.png if out_dir is given.
    """
    points = [
        dict(label=f"λ SDSS: {row['slope_sdss']:.2f} ± {row['slope_sdss_err']:.2f}",
             slope=row['slope_sdss'], slope_err=row['slope_sdss_err'],
             iz=row['i_z'], iz_err=row['i_z_err'], color='#d62728', name_text=name),
        dict(label=f"λ LSST: {row['slope_lsst']:.2f} ± {row['slope_lsst_err']:.2f}",
             slope=row['slope_lsst'], slope_err=row['slope_lsst_err'],
             iz=row['i_z'], iz_err=row['i_z_err'], color='#2ca02c', marker='o', size=120),
    ]
    tax = row['taxonomy']
    if tax == INDETERMINATE:
        tax = f"indeterminate ({row['candidates']})"
    elif tax in TAXO_BOXES:
        tax = f"{tax}  (P = {row['p_class']:.2f})"
    fig, ax = plt.subplots(figsize=(9, 5.2))
    plot_taxonomy(points, title=f"{name}  —  {tax}", ax=ax)
    ax.text(0.99, 0.02, f"i-z = {row['i_z']:.3f} ± {row['i_z_err']:.3f}",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=12,
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.6'))
    fig.tight_layout()
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        fig.savefig(os.path.join(out_dir, f"{name.replace(' ', '_')}_taxonomy.png"), dpi=dpi)
    return fig, ax


# --------------------------------------------------------------------------
# Literature taxonomy from SsODNet (rocks package)
# --------------------------------------------------------------------------
def _use_certifi():
    """Point SSL to certifi's certificates if no bundle is set (fixes SSL errors on Windows)."""
    if 'SSL_CERT_FILE' not in os.environ:
        try:
            import certifi
            os.environ['SSL_CERT_FILE'] = certifi.where()
        except ImportError:
            pass


def _first(lst, key):
    return lst[0].get(key) if lst else None


def query_rocks(names, cache_csv=None, refresh=False):
    """Taxonomy, diameter and albedo from the SsODNet ssoCard of each object (rocks.Rock).

    names     : designations (e.g. provid)
    cache_csv : if given, results are saved there and re-used on the next call
                (only the missing names are queried); refresh=True queries everything again.

    Returns a DataFrame indexed by the input name with columns:
    ssod_name, ssod_number, rocks_class, rocks_complex, rocks_scheme, rocks_technique,
    rocks_waverange, rocks_source, rocks_method, rocks_diameter, rocks_albedo, rocks_error.
    """
    names = list(dict.fromkeys(names))
    cached = pd.DataFrame()
    if cache_csv and os.path.exists(cache_csv) and not refresh:
        cached = pd.read_csv(cache_csv, index_col='name', keep_default_na=False,
                             na_values=['', 'nan', 'NaN'])
    todo = [n for n in names if n not in cached.index]

    rows = []
    if todo:
        _use_certifi()
        from rocks import config
        config.PATH_CACHE.mkdir(parents=True, exist_ok=True)  # skip the first-run prompt
        import rocks
        for n in todo:
            row = {'name': n}
            try:
                r = rocks.Rock(n)
                t = r.taxonomy.model_dump()
                row.update(
                    ssod_name=r.name, ssod_number=r.number,
                    rocks_class=t['class_']['value'] or None,
                    rocks_complex=t['complex']['value'] or None,
                    rocks_scheme=t['scheme']['value'] or None,
                    rocks_technique=t['technique']['value'] or None,
                    rocks_waverange=t['waverange']['value'] or None,
                    rocks_source=_first(t['bibref'], 'shortbib'),
                    rocks_method=_first(t['method'], 'shortbib'),
                    rocks_diameter=r.diameter.value, rocks_albedo=r.albedo.value,
                    rocks_error=None)
            except Exception as e:  # object not resolved, network error, ...
                row['rocks_error'] = f'{type(e).__name__}: {e}'[:200]
            rows.append(row)

    frames = [f for f in [cached, pd.DataFrame(rows).set_index('name') if rows else None]
              if f is not None and len(f)]
    out = pd.concat(frames) if frames else pd.DataFrame()
    out = out[~out.index.duplicated(keep='last')]
    if cache_csv and rows:
        out.to_csv(cache_csv, index_label='name')
    return out.reindex(names)
