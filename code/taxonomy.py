"""Taxonomy of asteroid colours with LSST parameters: reflectances, weighted gri slope and class boxes.

Each object is placed on the plane (gri slope [%/100nm], i-z [mag]) and assigned the class of the box
that contains it (boxes of Carruba et al. 2026).

Method
------
- Reflectance normalised to g for a band f, with the LSST AB solar magnitudes:
      R_f = 10^(-0.4 [(m_f - m_g) - (M_f_sun - M_g_sun)])      (R_g = 1)
  and its error by propagation, sigma_R_f = 0.4 ln10 R_f sigma_m_f.
- gri slope: weighted straight-line fit R(lam) = a + b lam through R_g, R_r, R_i at the LSST effective
  wavelengths (scipy curve_fit, sigma = reflectance errors); slope = 1e4 b [%/100nm], with its error,
  chi2 and p-value.
- Class: probability of each box from Monte Carlo draws of (slope, i-z) with their errors.

The comparison code (SDSS parameters, Carruba slope formula, unweighted fit, SsODNet queries) lives in
carruba_comparison/taxonomy_comparison.py, which is not part of the repository.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
# Effective wavelengths (nm)
LAM_LSST = {'u': 372.4, 'g': 480.7, 'r': 622.1, 'i': 755.9, 'z': 868.0, 'y': 975.3}

# Solar absolute magnitudes, LSST AB
M_SUN_LSST = {'g': -26.52, 'r': -26.93, 'i': -27.05, 'z': -27.07}

_K = 0.4 * np.log(10)  # d(10^(-0.4 x))/dx = -_K * 10^(-0.4 x)


# --------------------------------------------------------------------------
# Reflectances and weighted gri slope
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
SYS_SLOPE = 0.25   # %/100nm, systematic floor on the slope
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
