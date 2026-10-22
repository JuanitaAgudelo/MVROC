# Column Descriptions for `asteroid_analysis_results.csv`

This document describes the definition, physical meaning, and calculation method for each column in `asteroid_analysis_results.csv`, as written by `catalogue_pipeline.py`.

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
| `obs_duration` | `float` | Time span of the photometric series (days). | `epoch_utc_end` minus `epoch_utc_ini`, in days. |
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
| `power` | `float` | Lomb–Scargle power of the strongest peak. | Power of the first period in `periods_h`. When `multiband_photometry_compare` is selected in the same run, `power` is the power of that comparison (section 3d). |
| `amp_lsm` | `float` | Peak-to-peak amplitude of the Lomb–Scargle model in `band_periodogram` (mag). | \(\max - \min\) of the model over one cycle at the adopted frequency. |
| `amp_photo` | `float` | Photometric amplitude in `band_periodogram` (mag). | 95th minus 5th percentile of the magnitudes, \(P_{95} - P_{5}\). |

---

## 3. Voted period

Written when `periodogram_democratic` is selected. Later fits that need a starting period use this value when this analysis is selected and neither `phase_model_compare` nor `multiband_photometry_compare` is selected, and `periodogram_democratic_weighted` is not the only vote in the run. If both votes are selected and neither comparison is, later fits keep this unweighted period. The phase analysis writes a separate folded curve for this period, titled as democratic and unweighted. The vote uses the \(N_b\) filters with the most measurements and at least five points. \(u\) and \(y\) are excluded. The default is \(N_b = 3\).

Each of those filters casts one vote with its own Lomb–Scargle maximum. Every pair of those filters casts one joint multiband vote, and, when more than two filters are used, all of them together cast one more. For \(N_b = 3\) that is seven series. Periods that agree to 3%, or whose photometric frequencies differ by \(1 \pm 0.1\) or \(2 \pm 0.1\) cycles per day, count as the same vote. A unanimous vote reports the period of the joint series that includes the most filters. Otherwise the reported period is the Lomb–Scargle period of the filter with the most measurements among the single-filter series in the winning group.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_dem_h` | `float` | Rotation period chosen by the vote (h). | \(P_{\mathrm{rot}} = 48/f\) at the Lomb–Scargle frequency of the series selected above. |
| `votes_dem` | `int` | Number of series that voted for `p_rot_dem_h`. | Size of the winning group of agreeing periods. |
| `n_series_dem` | `int` | Number of series in the vote. | The \(N_b\) filters, every pair among them, and the joint series of all of them when \(N_b > 2\). |
| `bands_dem` | `str` | Filters that entered the vote, in order of the number of measurements. | For example `Lr+Lg+Li`. This is the filter list, not the subset that won. |
| `periods_dem_h` | `str` | Period adopted by every series. | Written as `Lr=3.3296,Lg=3.3296,Lr+Lg=3.3296`, with four digits after the decimal point. |

---

## 3b. Weighted alias period

Written when `periodogram_democratic_weighted` is selected. The reference filter is the one with the most measurements among the same \(N_b\) filters, again excluding \(u\) and \(y\) and requiring five points. Its Lomb–Scargle frequency \(f_{\mathrm{ref}}\) defines the candidates \(f_k = f_{\mathrm{ref}} + k\) for \(k = -6,\ldots,6\). A candidate outside the search window or below the cadence cut of the reference filter is dropped. The series are the same ones as in the unweighted vote: each filter, every pair, and all of the filters together when more than two are used. Each series scores every remaining candidate by its Lomb–Scargle power divided by its own maximum. A joint series uses the multiband periodogram. The alias score is the sum of those relative powers weighted by the fraction of points in the series:

\[
S(k) = \sum_b \frac{N_b}{N}\,\rho_b(f_k).
\]

The adopted alias maximizes \(S(k)\). If another alias lies within 5% of that score, the multiband Lomb–Scargle power of the same filters breaks the tie. The phase analysis writes a folded curve at this period, and the title says that the period is the weighted democratic one and gives \(k\). When this is the only vote selected and neither comparison is selected, later fits start from this period.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_wdem_h` | `float` | Rotation period chosen by the weighted alias score (h). | \(P_{\mathrm{rot}} = 48/f_{k^\star}\). |
| `k_alias_wdem` | `int` | Alias index of that period relative to the reference filter. | \(k^\star\) in \(f_k = f_{\mathrm{ref}} + k\), from \(-6\) to \(6\). Zero is the reference period itself. |
| `score_wdem` | `float` | Winning alias score. | \(S(k^\star)\). Each series contributes at most its fraction of the points counted over all series, so the score is at most 1. |
| `band_wdem` | `str` | Reference filter. | The filter with the most measurements. |
| `scores_wdem` | `str` | Score of every alias that was scored. | Written as `k=-1:0.640,k=+0:0.812`, in order of \(k\). |

---

## 3c. Phase-model comparison

Written when `phase_model_compare` is selected. It needs one of the two votes. The candidates are the photometric frequency of the single-band periodogram and the frequency adopted by every series in the selected vote. Frequencies within 0.1 cycles/day are merged. Each candidate is folded on the periodogram band and compared with the Lomb–Scargle model of that band at the same frequency. The score is the reduced chi-squared of that model,

\[
\chi^2_\nu = \frac{1}{\nu}\sum_i\left(\frac{y_i - m(t_i)}{\sigma_i}\right)^2,
\]

where \(m(t)\) is the model and \(\nu\) is the number of points minus the number of model coefficients. The period with the lowest \(\chi^2_\nu\) is the one used by every later fit when `multiband_photometry_compare` is not selected: the folded model, the shape fit, the time-domain model, the literature comparison, and the multiband photometry. The comparison figure draws one panel per candidate, with the selected period marked in its title.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_cmp_h` | `float` | Rotation period with the lowest reduced chi-squared (h). | \(P_{\mathrm{rot}} = 48/f\) for the candidate that minimizes \(\chi^2_\nu\). |
| `chi2_nu_cmp` | `float` | Reduced chi-squared of that period. | \(\chi^2_\nu\) of the periodogram-band model at `p_rot_cmp_h`. |
| `periods_cmp` | `str` | Reduced chi-squared of every candidate, best first. | Written as `3.0550:1.80,4.4491:24.30`. |

---

## 3d. Multiband comparison with a two-peak model

Written when `multiband_photometry_compare` is selected. It needs `periodogram_democratic`. The candidates are the periods adopted by the series of that vote. The observations span \(T\) days. The campaign aliases are the frequencies \(1/T + k\) and \(2/T + k\). The integer \(k\) runs from the long edge of the period search, a rotation of \(48\)~h, down to a rotation of \(2\)~h. With \(P_{\mathrm{rot}} = 48/f\), those edges are the frequencies \(1\) and \(24\,\mathrm{d}^{-1}\). \(k_{\min} = \lceil 1 - 2/T \rceil\) and \(k_{\max} = \lceil 24 - 1/T \rceil\), and both families use every integer between them. A shift that falls outside the periodogram window, or below the cadence cut, is dropped. A voted frequency that lies within 0.1 cycles/day of an alias receives that label and keeps the identity of the series that voted for it. The chi-squared is computed for the voted frequencies only. An alias that no series voted for is drawn on the periodogram and is not fitted. At each voted candidate the filters \(g\), \(r\), \(i\), \(z\), and \(y\) share one harmonic model of order \(K = 2\), with a zero point per filter. Magnitudes outside the 5th to 95th percentile of their filter, extended by the same margin used in the phase curve, are left out of the fit and of the figure. Order 2 is the model with at most two brightness maxima per rotation. \(u\) is left out. The score is the reduced chi-squared of that fit on the same points for every candidate. The adopted period is the lowest \(\chi^2_\nu\) among the candidates whose model has two peaks and that are not aliases of the campaign duration. A voted period that is an alias is scored and is not adopted: a low score there means the nights follow the campaign window. If every two-peak model is an alias, the asteroid is left out of the catalogue and listed in `asteroid_processing_report.txt` under “Only campaign aliases have two peaks”. If no candidate has two peaks, it is listed under “No two-peak period in the multiband comparison”. The figure is for a visual check. The adopted period occupies a large full-width panel on top. Below it, in two columns, are the other voted periods whose model has two peaks and a higher \(\chi^2_\nu\). A model with one peak or with none is scored in the log and left off the figure. A voted period that also lands on an alias is drawn with a red title. When a period is adopted, it is the one used by the later fits, including the full multiband photometry. Each panel title gives the Lomb–Scargle power at that photometric frequency, the peak-to-peak amplitude \(A\) of the order-2 model, and the ratio \(A/s\): \(A\) is the peak-to-peak amplitude of the order-2 model, and \(s\) is the mean distance of the points to that curve. The farthest 10% of those distances are left out. Each remaining point contributes \(\sqrt{d^2+\sigma^2}\), so a large magnitude error widens \(s\) even when the point sits on the curve. `timemodel_Lk` draws that same period through the nights of the synthetic \(L_k\) light curve, with a two-harmonic model fit to \(L_k\). It needs the \(L_k\) series built from this period. The selected panel draws dotted lines at the 2.5th and 97.5th percentiles of the plotted magnitudes. When `compare` is selected, `p_rot_h` is plotted against the literature Lomb–Scargle period and against the literature Fourier period. The adopted model also yields the color indices. Each filter zero point is the magnitude of that filter minus the reference filter at the same rotational phase. The reference is \(r\) when \(r\) has data, and otherwise the filter with the most points. The filters are ordered by wavelength as \(u\), \(g\), \(r\), \(i\), \(z\), \(y\). Every pair among the filters that have a zero point is stored, with the bluer filter first. With \(g\), \(r\), and \(i\) the colors are \(g-r\), \(g-i\), and \(r-i\). With a fourth filter the three extra pairs that include it are stored as well. A color is the difference of the two zero points. When both zero points are free, the uncertainty is \(\sqrt{\mathrm{Var}(c_a)+\mathrm{Var}(c_b)-2\mathrm{Cov}(c_a,c_b)}\). When one filter is the reference, the uncertainty is the formal error of the free coefficient. \(u\) is measured afterwards, as the weighted mean of \(H_u\) minus the reference model at the same times, after a \(5\sigma\) clip. Each color \(u\) minus another filter is that offset minus the other filter's zero point, and the two uncertainties are added in quadrature. A pair is empty when either filter has no measurement.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_h` | `float` | Rotation period chosen by the two-peak comparison (h). | \(P_{\mathrm{rot}} = 48/f\) for the candidate that minimizes \(\chi^2_\nu\) among the two-peak models. |
| `chi2_nu` | `float` | Reduced chi-squared of that period. | \(\chi^2_\nu\) of the order-2 joint model. \(\nu\) is the number of fitted points minus the number of coefficients. |
| `n_peaks` | `int` | Brightness maxima of the adopted model. | Peaks of \(-\)magnitude over one rotation. The phase is treated as circular. A peak must rise at least 0.02 mag and at least 15% of the model amplitude. |
| `power` | `float` | Lomb–Scargle power of the adopted period. | Normalized multiband Lomb–Scargle power of the fitted points at the photometric frequency \(f = 48/P_{\mathrm{rot}}\). This is the value stored in `power` when this analysis is selected. |
| `amp_ratio` | `float` | Amplitude of the adopted model over the scatter about the curve. | \(A/s\). \(A\) is the peak-to-peak of the order-2 model, the value shown in the panel title. \(s\) is the mean of \(\sqrt{d^2+\sigma^2}\) after dropping the farthest 10% of the distances \(d\) from the points to the curve. |
| `amp_data` | `float` | Amplitude of the plotted magnitudes (mag). | Difference between the 97.5th and the 2.5th percentile of the points shifted onto the reference filter. The same two levels are the dotted lines on the selected panel. |
| `span_alias` | `int` | Whether the adopted period is an alias of the campaign duration. | 0 for every stored period. A voted period that lands on a campaign alias is marked and is not adopted, so a row that has a period has this flag equal to 0. |
| `u-g`, `u-g_err` | `float` | \(u-g\) and its uncertainty (mag). | \(c_u - c_g\). \(c_u\) is the post-fit mean of \(u\) minus the reference curve. The uncertainties are added in quadrature. |
| `u-r`, `u-r_err` | `float` | \(u-r\) and its uncertainty (mag). | \(c_u - c_r\). When \(r\) is the reference, \(c_r = 0\) and this is \(c_u\). |
| `u-i`, `u-i_err` | `float` | \(u-i\) and its uncertainty (mag). | \(c_u - c_i\), uncertainties in quadrature. |
| `u-z`, `u-z_err` | `float` | \(u-z\) and its uncertainty (mag). | \(c_u - c_z\), uncertainties in quadrature. |
| `u-y`, `u-y_err` | `float` | \(u-y\) and its uncertainty (mag). | \(c_u - c_y\), uncertainties in quadrature. |
| `g-r`, `g-r_err` | `float` | \(g-r\) and its uncertainty (mag). | \(c_g - c_r\). A positive value means \(g\) is fainter than \(r\) at the same phase. When \(r\) is the reference this is the fitted offset of \(g\). |
| `g-i`, `g-i_err` | `float` | \(g-i\) and its uncertainty (mag). | \(c_g - c_i\), with the covariance of the two zero points when both are free. |
| `g-z`, `g-z_err` | `float` | \(g-z\) and its uncertainty (mag). | \(c_g - c_z\). |
| `g-y`, `g-y_err` | `float` | \(g-y\) and its uncertainty (mag). | \(c_g - c_y\). |
| `r-i`, `r-i_err` | `float` | \(r-i\) and its uncertainty (mag). | \(c_r - c_i\). When \(r\) is the reference this is minus the fitted offset of \(i\). |
| `r-z`, `r-z_err` | `float` | \(r-z\) and its uncertainty (mag). | \(c_r - c_z\). |
| `r-y`, `r-y_err` | `float` | \(r-y\) and its uncertainty (mag). | \(c_r - c_y\). |
| `i-z`, `i-z_err` | `float` | \(i-z\) and its uncertainty (mag). | \(c_i - c_z\). |
| `i-y`, `i-y_err` | `float` | \(i-y\) and its uncertainty (mag). | \(c_i - c_y\). |
| `z-y`, `z-y_err` | `float` | \(z-y\) and its uncertainty (mag). | \(c_z - c_y\). |

---

## 3e. Fold at the catalogue period

Written when `justfold` is selected. It does not search for a period. The period is the one stored for the asteroid in `asteroids-with-rotation.csv`. The Lomb–Scargle value and the Fourier value are both drawn when they differ. When they agree, the figure has one panel. The model is the same joint harmonic fit of order \(K = 2\) used by the multiband comparison: filters \(g\), \(r\), \(i\), \(z\), and \(y\) share the curve, each with its own zero point, and \(u\) is left out. The same magnitude outliers are rejected. A panel title is red when that period is a daily alias of the campaign duration. `select_provid_rotfold()` runs this analysis alone. The scalar columns below refer to the first panel, which is the Lomb–Scargle period when the two values differ.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_justfold_h` | `float` | Catalogue rotation period drawn in the first panel (h). | `period_lsm_hr` when it is present. `period_fourier_hr` otherwise. |
| `chi2_nu_justfold` | `float` | Reduced chi-squared of that fold. | \(\chi^2_\nu\) of the order-2 joint model at the catalogue period. |
| `n_peaks_justfold` | `int` | Brightness maxima of that model. | Same peak rule as `n_peaks`. |
| `power_justfold` | `float` | Lomb–Scargle power at that period. | Normalized multiband power at \(f = 48/P_{\mathrm{rot}}\). |
| `amp_ratio_justfold` | `float` | Amplitude of the model over the scatter about the curve. | Same \(A/s\) as `amp_ratio`. |
| `span_alias_justfold` | `int` | Whether that period is an alias of the campaign duration. | 1 when the frequency lies within 0.1 cycles/day of \(1/T + k\) or \(2/T + k\). |
| `periods_justfold` | `str` | Every catalogue period that was drawn. | Written as `lsm:8.9000:2.10:2:0.650:1.200:0` or `lsm+fourier:...` when the two values agree. The fields are source, period, \(\chi^2_\nu\), peaks, power, \(A/s\), and the baseline-alias flag. |

---

## 4. Single-band shape fit

Written when `shape_fit` is selected. The definitions are the same as the \(L_k\) columns in section 6, without the `_k` suffix. The measurements are those in `band_periodogram`, folded at the period chosen by `phase_model_compare` when that analysis is selected, otherwise at the voted period, and at the single-band periodogram when no vote is selected.

- Touching spheres: `chi_r_sph`, `q_sph`, `q_sph_err`, `phi_sph_deg`, `phi_sph_err_deg`, `amp_sph`, `m_mid_sph`, `m_sep`, `p_rot_sph_h`.
- Equator-on ellipsoid: `chi_r_ell`, `b_over_a`, `b_over_a_err`, `amp_ell`, `m_mid_ell`, `m_ab`, `p_rot_ell_h`.

---

## 5. Multiband photometry and synthetic channel \(L_k\)

Written when `multiband_photometry` is selected. The filters \(g\), \(r\), \(i\), \(z\), and \(y\) share one harmonic model. \(u\) is left out of that model and is compared with it afterwards. The reference filter is \(r\) when \(r\) has measurements; otherwise it is the filter with the most measurements. When `multiband_photometry_compare` is selected in the same run, every color pair from that comparison stays (section 3d). The adjacent colors computed in this section are not written over those columns.

The starting photometric frequency \(f_{\mathrm{ph0}}\) is the frequency chosen by `multiband_photometry_compare` when that analysis is selected, otherwise `phase_model_compare`, otherwise the voted period, otherwise the single-band periodogram. The rotational frequency is then refined: for each harmonic order \(K\) from 2 to 10, the fit searches a grid of step \(0.0002\,\mathrm{d}^{-1}\) within \(\pm 0.02\,\mathrm{d}^{-1}\) of \(f_{\mathrm{ph0}}/2\). The adopted pair \((f_{\mathrm{rot}}, K)\) is the one with the lowest Bayesian information criterion. At each epoch the filters are shifted onto that common curve and their fluxes are added. \(L_k\) is the magnitude of the sum.

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

Written when `shape_fit_multiband` is selected. It runs after `multiband_photometry_compare` when that analysis is selected: the comparison chooses the period, \(L_k\) is built at that period, and then the shape is fit. Both models are fit to the \(L_k\) magnitudes by non-linear least squares. The starting frequency is the photometric frequency of that multiband fit, and the fit may change it by up to 2%. The reduced chi-squared is \(\chi^2/(N-p)\), with \(p\) the number of free parameters.

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
