# Column Descriptions for `asteroid_analysis_results.csv`

This document describes the definition, physical meaning, and calculation method for each column in `asteroid_analysis_results.csv`, as written by `plot_rotation_figures.py`.

Each row is one asteroid. Every file starts with the metadata columns. The other blocks are written only when that analysis is selected. `select_provid_medellin`, the selection used for the Medellín catalogue, writes the voted period, the joint multiband fit, the shape fits on the synthetic channel \(L_k\), and the reflectance. It does not write the single-band periodogram or the single-band shape fit.

The light curve is taken to be double-peaked, so a rotational period is twice the photometric period. A photometric frequency \(f\) in cycles per day corresponds to \(P_{\mathrm{rot}} = 48/f\) hours and \(P_{\mathrm{phot}} = 24/f\) hours.

---

## 1. Mandatory columns

Present at the start of every row:

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `provid` | `str` | Asteroid provisional designation from the Minor Planet Center (MPC). | Object identifier (e.g., `2002 TH52`). |
| `epoch_utc_ini` | `str` | UTC timestamp of the earliest photometric measurement. | Earliest observation time in the data, formatted as `YYYY-MM-DD HH:MM:SS`. |
| `epoch_utc_end` | `str` | UTC timestamp of the latest photometric measurement. | Latest observation time in the data, formatted as `YYYY-MM-DD HH:MM:SS`. |
| `n_obs` | `int` | Total number of photometric measurements of the asteroid, in all filters. | Number of measurements in the data for that asteroid. |
| `n_obs_u` | `int` | Number of photometric measurements in the Rubin \(u\) filter. | Number of measurements in the \(u\) band. |
| `n_obs_g` | `int` | Number of photometric measurements in the Rubin \(g\) filter. | Number of measurements in the \(g\) band. |
| `n_obs_r` | `int` | Number of photometric measurements in the Rubin \(r\) filter. | Number of measurements in the \(r\) band. |
| `n_obs_i` | `int` | Number of photometric measurements in the Rubin \(i\) filter. | Number of measurements in the \(i\) band. |
| `n_obs_z` | `int` | Number of photometric measurements in the Rubin \(z\) filter. | Number of measurements in the \(z\) band. |
| `n_obs_y` | `int` | Number of photometric measurements in the Rubin \(y\) filter. | Number of measurements in the \(y\) band. |

---

## 2. Single-band periodogram

Written when `periodogram` is selected. These columns use one filter, `band_periodogram`. It is \(r\) (`Lr`) unless the selection asks for another band. The request `Lx` means the filter with the most measurements, excluding \(u\) and \(y\); \(r\) is kept when it ties with another filter.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `band_periodogram` | `str` | Filter used for the periodogram, the phase curve, the phase model, the single-band shape fit, and the nightly time model. | For example `Lr` or `Li`. |
| `periods_h` | `str` | Candidate rotation periods in hours, separated by commas. | Lomb–Scargle peaks with power of at least 80% of the strongest peak. The strongest period is listed first. Each stored period is twice the photometric period of that peak. |
| `power` | `float` | Lomb–Scargle power of the strongest peak. | Power of the first period in `periods_h`. |
| `amp_lsm` | `float` | Peak-to-peak amplitude of the Lomb–Scargle model in `band_periodogram` (mag). | \(\max - \min\) of the model over one cycle at the adopted frequency. |
| `amp_photo` | `float` | Photometric amplitude in `band_periodogram` (mag). | 95th minus 5th percentile of the magnitudes, \(P_{95} - P_{5}\). |

---

## 3. Voted period

Written when `periodogram_democratic` is selected. The phase curve and any later fit that needs a starting period use this value. The vote uses the \(N_b\) filters with the most measurements and at least five points. \(u\) and \(y\) are excluded. The default is \(N_b = 3\), and the Medellín selection uses that default.

Each of those filters casts one vote with its own Lomb–Scargle maximum. Each pair of consecutive filters in that list casts one joint multiband vote, and, when more than two filters are used, all of them together cast one more. For \(N_b = 3\) that is six series. Periods that agree to 3%, or whose photometric frequencies differ by \(1 \pm 0.1\) or \(2 \pm 0.1\) cycles per day, count as the same vote. A unanimous vote reports the period of the joint series that includes the most filters. Otherwise the reported period is the Lomb–Scargle period of the filter with the most measurements among the single-filter series in the winning group.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_dem_h` | `float` | Rotation period chosen by the vote (h). | \(P_{\mathrm{rot}} = 48/f\) at the Lomb–Scargle frequency of the series selected above. |
| `votes_dem` | `int` | Number of series that voted for `p_rot_dem_h`. | Size of the winning group of agreeing periods. |
| `n_series_dem` | `int` | Number of series in the vote. | The \(N_b\) filters, their consecutive pairs, and the joint series of all of them when \(N_b > 2\). |
| `bands_dem` | `str` | Filters that entered the vote, in order of the number of measurements. | For example `Lr+Lg+Li`. This is the filter list, not the subset that won. |
| `periods_dem_h` | `str` | Period adopted by every series. | Written as `Lr=3.3296,Lg=3.3296,Lr+Lg=3.3296`, with four digits after the decimal point. |

---

## 4. Single-band shape fit

Written when `shape_fit` is selected. The definitions are the same as the \(L_k\) columns in section 6, without the `_k` suffix. The measurements are those in `band_periodogram`, folded at the voted period when that vote is part of the analysis and at the single-band periodogram otherwise.

- Touching spheres: `chi_r_sph`, `q_sph`, `q_sph_err`, `phi_sph_deg`, `phi_sph_err_deg`, `amp_sph`, `m_mid_sph`, `m_sep`, `p_rot_sph_h`.
- Equator-on ellipsoid: `chi_r_ell`, `b_over_a`, `b_over_a_err`, `amp_ell`, `m_mid_ell`, `m_ab`, `p_rot_ell_h`.

---

## 5. Multiband photometry and synthetic channel \(L_k\)

Written when `multiband_photometry` is selected. The filters \(g\), \(r\), \(i\), \(z\), and \(y\) share one harmonic model. \(u\) is left out of that model and is compared with it afterwards. The reference filter is \(r\) when \(r\) has measurements; otherwise it is the filter with the most measurements.

The starting photometric frequency \(f_{\mathrm{ph0}}\) is the frequency passed in by the periodogram or by the voted period. The rotational frequency is then refined: for each harmonic order \(K\) from 2 to 10, the fit searches a grid of step \(0.0002\,\mathrm{d}^{-1}\) within \(\pm 0.02\,\mathrm{d}^{-1}\) of \(f_{\mathrm{ph0}}/2\). The adopted pair \((f_{\mathrm{rot}}, K)\) is the one with the lowest Bayesian information criterion. At each epoch the filters are shifted onto that common curve and their fluxes are added. \(L_k\) is the magnitude of the sum.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_multi_used_h` | `float` | Rotation period, in hours, used to start the multiband fit. | \(P_{\mathrm{rot}} = 48/f_{\mathrm{ph0}}\). This is the voted period when that vote is part of the analysis, and the single-band periodogram period otherwise. |
| `p_rot_multi_h` | `float` | Rotational period adopted by the multiband fit (h). | \(P_{\mathrm{rot}} = 24/f_{\mathrm{rot}}\), where \(f_{\mathrm{rot}}\) is the refined rotational frequency. |
| `f_rot_multi` | `float` | Rotational frequency adopted by the multiband fit (\(\mathrm{d}^{-1}\)). | The grid frequency, near \(f_{\mathrm{ph0}}/2\), that minimizes the BIC together with \(K\). |
| `k_harm_multi` | `int` | Number of Fourier harmonics in the joint fit. | The value of \(K\) between 2 and 10 that minimizes the BIC. |
| `bic` | `float` | Bayesian information criterion of the adopted harmonic model. | \(\chi^2 + k \ln N\). Here \(k = 2K + N_{\mathrm{filters}}\) and \(N\) is the number of points kept by the preliminary \(5\sigma\) clip at \(f_{\mathrm{ph0}}\) with \(K = 3\). \(\chi^2\) is the weighted chi-squared at the chosen grid frequency. |
| `amp_lsm_k` | `float` | Peak-to-peak amplitude of the Lomb–Scargle model on \(L_k\) (mag). | \(\max(y_{\mathrm{model}}) - \min(y_{\mathrm{model}})\) over one rotation, using \(K\) harmonics of \(f_{\mathrm{rot}}\). |
| `amp_photo_k` | `float` | Empirical amplitude of the \(L_k\) magnitudes after rejecting outliers (mag). | Residuals to that model are clipped at \(3 \times \mathrm{MAD}\). The amplitude is \(P_{99}(y_{\mathrm{clean}}) - P_{1}(y_{\mathrm{clean}})\). If fewer than five points remain, all \(L_k\) points are used. |
| `rms_multi` | `float` | Root-mean-square residual of the joint harmonic fit (mag). | \(\sqrt{\langle r^2 \rangle}\) of the magnitude residuals of the points kept by the clip. The residuals are not divided by their uncertainties. |
| `n_outliers_multi` | `int` | Number of measurements rejected in the joint fit. | Points in \(g\), \(r\), \(i\), \(z\), and \(y\) removed by the iterative clip at \(5 \times \mathrm{MAD}\). Measurements in \(u\) are not included. |
| `u-g`, `u-g_err` | `float` | \(u-g\) color and its 1-sigma uncertainty (mag). | \((u-r) - (g-r)\) when \(u\) is measured against the model and \(g\) is in the fit. The two uncertainties are added in quadrature. |
| `g-r`, `g-r_err` | `float` | \(g-r\) color and its 1-sigma uncertainty (mag). | Difference of the \(g\) and \(r\) magnitude offsets. When \(r\) is the reference, this is the fitted offset of \(g\). |
| `r-i`, `r-i_err` | `float` | \(r-i\) color and its 1-sigma uncertainty (mag). | Difference of the \(r\) and \(i\) magnitude offsets. When \(r\) is the reference, this is minus the fitted offset of \(i\). |
| `i-z`, `i-z_err` | `float` | \(i-z\) color and its 1-sigma uncertainty (mag). | Difference of the \(i\) and \(z\) offsets. When both offsets are free, \(\sigma = \sqrt{\mathrm{Var}(c_i)+\mathrm{Var}(c_z)-2\mathrm{Cov}(c_i,c_z)}\). |
| `z-y`, `z-y_err` | `float` | \(z-y\) color and its 1-sigma uncertainty (mag). | Difference of the \(z\) and \(y\) offsets, with the covariance when both offsets are free. |
| `u-r`, `u-r_err` | `float` | Color of the \(u\) measurements relative to the reference filter, and its uncertainty (mag). | Weighted mean of \(u\) minus the reference model, after the same \(5\sigma\) clip. The column name uses the reference filter, which is \(r\) when \(r\) is present. |

A color is absent when one of the two filters has no measurement. The uncertainty of an adjacent pair is the formal error of the single free offset when the other filter is the reference.

---

## 6. Shape fits on channel \(L_k\)

Written when `shape_fit_multiband` is selected. Both models are fit to the \(L_k\) magnitudes by non-linear least squares. The starting frequency is the photometric frequency of the multiband fit, and the fit may change it by up to 2%. The reduced chi-squared is \(\chi^2/(N-p)\), with \(p\) the number of free parameters.

### A. Touching-spheres model

The asteroid is modeled as two tangent spheres of radii \(r_1 = 1\) and \(r_2 = q\), rotating about the axis perpendicular to the line that joins their centers. The projected area includes mutual occultation. Five parameters are free: the separated magnitude, \(q\), the aspect, the initial phase, and the frequency.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `chi_r_sph_k` | `float` | Reduced chi-squared of the touching-spheres fit on \(L_k\). | \(\chi^2/(N-5)\). |
| `q_sph_k` | `float` | Ratio of the two lobe radii (\(r_2/r_1 \le 1\)). | \(\min(q, 1/q)\) from the least-squares fit, so the smaller sphere is reported second. |
| `q_sph_err_k` | `float` | Formal 1-sigma uncertainty of `q_sph_k`. | Uncertainty of \(q\) from the covariance. It is multiplied by \(1/q^2\) when the fitted \(q\) was greater than 1. |
| `phi_sph_deg_k` | `float` | Angle between the spin axis and the line of sight (deg). | \(\phi = \arccos\sqrt{\cos^2\phi}\). \(\phi = 90^\circ\) is equator-on and \(\phi = 0^\circ\) is pole-on. This is not the rotational phase. |
| `phi_sph_err_deg_k` | `float` | 1-sigma uncertainty of that angle (deg). | Propagated from the covariance of \(\cos^2\phi\). NaN when \(\cos^2\phi\) is consistent with zero, so the fit only gives a lower limit, usually near \(90^\circ\). |
| `amp_sph_k` | `float` | Peak-to-peak amplitude of the touching-spheres model (mag). | Faintest minus brightest magnitude of the model, \(\Delta m = -2.5\log_{10}(A_{\min}/A_{\max})\). |
| `m_mid_sph_k` | `float` | Midpoint of the model light curve (mag). | Halfway between the brightest and faintest magnitudes of the model. |
| `m_sep_k` | `float` | Magnitude when the two disks are fully separated. | Fitted magnitude of the maximum projected area, \(\pi(1+q^2)\). |
| `p_rot_sph_h_k` | `float` | Rotation period of the touching-spheres fit (h). | Twice the photometric period of the fitted frequency. |

### B. Equator-on ellipsoid

The asteroid is modeled as an ellipsoid with long axis \(a = 1\), seen equator-on and rotating about its short axis. The other two axes are equal, \(c = b\), so the only shape parameter is \(b/a\). Four parameters are free: the magnitude of the \(a\)–\(b\) face, \(b/a\), the initial phase, and the frequency.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `chi_r_ell_k` | `float` | Reduced chi-squared of the ellipsoid fit on \(L_k\). | \(\chi^2/(N-4)\). |
| `b_over_a_k` | `float` | Equatorial axis ratio (\(b/a \le 1\)). | Elongation that produces the double-peaked light curve. |
| `b_over_a_err_k` | `float` | Formal 1-sigma uncertainty of \(b/a\). | Taken from the covariance of the fit. |
| `amp_ell_k` | `float` | Peak-to-peak amplitude of the ellipsoid model (mag). | Faintest minus brightest magnitude of the model. For this equator-on biaxial ellipsoid, \(\Delta m = 2.5\log_{10}(a/b)\). |
| `m_mid_ell_k` | `float` | Midpoint of the model light curve (mag). | Halfway between the brightest and faintest magnitudes of the model. |
| `m_ab_k` | `float` | Magnitude of the \(a\)–\(b\) face. | Fitted magnitude of the maximum projected area, the brightest orientation when seen equator-on. |
| `p_rot_ell_h_k` | `float` | Rotation period of the ellipsoid fit (h). | Twice the photometric period of the fitted frequency. |

---

## 7. Reflectance spectrum

Written when `spectrum` is selected. The colors come from the joint light-curve fit, relative to \(r\). The shared rotation model removes the rotational phase, so the two filters do not have to be observed on the same night. The reflectance is normalized to the Sun and to the \(r\) band, where \(R(r) = 1\).

For a filter \(b\), the color excess is \((m_b - m_r) - (m_b - m_r)_\odot\), and

\[
R(b)/R(r) = 10^{-0.4 \times \mathrm{excess}}.
\]

The solar colors \((m_b-m_r)_\odot\) are \(g-r = 0.45\), \(i-r = -0.12\), \(z-r = -0.16\), and \(y-r = -0.18\). The uncertainty is \(R \times 0.4 \ln(10) \times \sigma(m_b-m_r)\).

For the \(g\), \(i\), \(z\), and \(y\) filters, written below as \(b\):

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `{b}-r_spec` | `float` | Color \(m_b - m_r\) (mag). | The joint-fit color relative to \(r\). For \(g\) this equals `g-r`. For \(i\) it equals `-(r-i)`. |
| `{b}-r_spec_err` | `float` | Uncertainty of that color (mag). | Uncertainty propagated from the joint fit. |
| `R_{b}` | `float` | Solar-normalized reflectance \(R(b)/R(r)\). | \(10^{-0.4 \times \mathrm{excess}}\) with the solar color above. The stored names are `R_g`, `R_i`, `R_z`, and `R_y`. |
| `R_{b}_err` | `float` | Uncertainty of the reflectance. | \(R \times 0.4 \ln(10) \times \sigma(m_b-m_r)\). |
| `n_{b}_spec` | `int` | Number of measurements in filter \(b\). | Count of measurements in that filter. It is 0 when the color is not available. |

The stored color names are `g-r_spec`, `i-r_spec`, `z-r_spec`, and `y-r_spec`, each with an `_err` column.
