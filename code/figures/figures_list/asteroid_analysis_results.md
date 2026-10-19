# Column Descriptions for `asteroid_analysis_results.csv`

This document describes the definition, physical meaning, and calculation method for each column in `asteroid_analysis_results.csv`.

Each row is one asteroid. The columns give the time span and number of photometric measurements, the rotation periods and amplitudes, the colors, and the parameters of the shape models.

---

## 1. Mandatory Columns (Metadata and Observation Counts)

Present at the start of every row:

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `provid` | `str` | Asteroid provisional designation from the Minor Planet Center (MPC). | Object identifier (e.g., `2002 TH52`). |
| `epoch_utc_ini` | `str` | UTC timestamp of the earliest photometric measurement. | Earliest observation time in the data, formatted as `YYYY-MM-DD HH:MM:SS`. |
| `epoch_utc_end` | `str` | UTC timestamp of the latest photometric measurement. | Latest observation time in the data, formatted as `YYYY-MM-DD HH:MM:SS`. |
| `n_obs` | `int` | Total number of photometric measurements of the asteroid, in all filters. | Number of measurements in the data for that asteroid. |
| `n_obs_u` | `int` | Number of photometric measurements in the Rubin $u$ filter. | Number of measurements in the $u$ band. |
| `n_obs_g` | `int` | Number of photometric measurements in the Rubin $g$ filter. | Number of measurements in the $g$ band. |
| `n_obs_r` | `int` | Number of photometric measurements in the Rubin $r$ filter. | Number of measurements in the $r$ band. |
| `n_obs_i` | `int` | Number of photometric measurements in the Rubin $i$ filter. | Number of measurements in the $i$ band. |
| `n_obs_z` | `int` | Number of photometric measurements in the Rubin $z$ filter. | Number of measurements in the $z$ band. |
| `n_obs_y` | `int` | Number of photometric measurements in the Rubin $y$ filter. | Number of measurements in the $y$ band. |

---

## 2. Multiband Photometry and Synthetic Channel $L_k$

The measurements in the six filters are combined into the synthetic channel $L_k$.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `p_rot_multi_used_h` | `float` | Rotation period, in hours, used to start the multiband fit. | Twice the photometric period of the frequency passed in: $P_{\text{rot}} = 48 / f_{\text{ph0}}$. This is the periodogram period, or the voted period when that vote is part of the analysis. |
| `p_rot_multi_h` | `float` | Physical rotational period in hours. | Derived from the rotational frequency: $P_{\text{rot}} = 24 / f_{\text{rot}}$. The light curve is taken to be double-peaked ($P_{\text{rot}} = 2 \times P_{\text{phot}}$). |
| `f_rot_multi` | `float` | Physical rotational frequency in cycles per day ($1/\text{day}$). | $f_{\text{rot}} = f_{\text{ph0}} / 2$, where $f_{\text{ph0}}$ is the primary photometric frequency from the periodogram. |
| `k_harm_multi` | `int` | Number of Fourier harmonics used in the joint multi-filter fit. | The value of $K$ between 2 and 10 that minimizes the Bayesian Information Criterion (BIC) at $f_{\text{rot}}$. |
| `bic` | `float` | Bayesian Information Criterion of the harmonic model adopted for channel $L_k$. | Minimum of $\chi^2 + k \ln N$ over harmonic orders $K = 2$ to $10$. Here $k$ is the number of fit parameters and $N$ is the number of measurements kept in the fit. |
| `amp_lsm_k` | `float` | Peak-to-peak amplitude of the Lomb-Scargle harmonic model on channel $L_k$ (mag). | Difference between the maximum and minimum of the model over a full phase cycle $[0, 1)$ at $f_{\text{rot}}$: $\text{amp\_lsm\_k} = \max(y_{\text{model}}) - \min(y_{\text{model}})$. |
| `amp_photo_k` | `float` | Empirical photometric amplitude on channel $L_k$, after rejecting outliers (mag). | Residuals to the harmonic model ($r_i = y_i - y_{\text{model}}(t_i)$) are clipped at $3\sigma$ ($3 \times \text{MAD}$). The amplitude is the 99th minus the 1st percentile: $\text{amp\_photo\_k} = P_{99}(y_{\text{clean}}) - P_{1}(y_{\text{clean}})$. |
| `rms_multi` | `float` | Weighted root-mean-square residual of the joint harmonic fit (mag). | Standard deviation of the weighted residuals of the multi-filter fit. |
| `n_outliers_multi` | `int` | Number of measurements rejected as outliers in the multi-filter fit. | Total number of measurements minus those that remain after iterative $3\sigma$ clipping. |
| `u-g`, `u-g_err` | `float` | $u - g$ color and its 1-sigma uncertainty (mag). | Magnitude offset between $u$ and $g$. When $u$ is left out of the periodic fit, the offset is measured against the model and the uncertainty of $g$ is included. |
| `g-r`, `g-r_err` | `float` | $g - r$ color and its 1-sigma uncertainty (mag). | Difference of the $g$ and $r$ magnitude offsets, with the covariance of the fit. |
| `r-i`, `r-i_err` | `float` | $r - i$ color and its 1-sigma uncertainty (mag). | Difference of the $r$ and $i$ magnitude offsets, with the covariance of the fit. |
| `i-z`, `i-z_err` | `float` | $i - z$ color and its 1-sigma uncertainty (mag). | Difference of the $i$ and $z$ magnitude offsets: $\sigma = \sqrt{\text{Var}(c_i) + \text{Var}(c_z) - 2\text{Cov}(c_i, c_z)}$. |
| `z-y`, `z-y_err` | `float` | $z - y$ color and its 1-sigma uncertainty (mag). | Difference of the $z$ and $y$ magnitude offsets, with the covariance of the fit. |
| `u-r`, `u-r_err` | `float` | Color of the $u$ measurements relative to the reference filter $r$, and its uncertainty (mag). | Magnitude offset between $u$ and $r$. |

---

## 3. Shape Fits on Channel $L_k$

The shape models are fit to the combined magnitudes $L_k$.

### A. Touching-spheres model

The asteroid is modeled as two tangent spheres of radii $r_1 = 1$ and $r_2 = q \le 1$, rotating about the axis perpendicular to the line joining their centers. The projected area includes mutual occultation.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `chi_r_sph_k` | `float` | Reduced chi-squared ($\chi^2/\nu$) of the touching-spheres fit on channel $L_k$. | $\chi^2 / (N - 3)$ for the radius ratio $q$, the initial phase, and the reference magnitude. |
| `q_sph_k` | `float` | Ratio of the two lobe radii ($r_2 / r_1 \le 1$). | Value from the non-linear least-squares fit. |
| `q_sph_err_k` | `float` | Formal 1-sigma uncertainty of the radius ratio $q$. | Taken from the covariance of the fit. |
| `phi_sph_deg_k` | `float` | Initial rotational phase angle, in degrees. | Phase angle at the time of the first measurement, wrapped to $[0^\circ, 360^\circ)$. |
| `phi_sph_err_deg_k` | `float` | 1-sigma uncertainty of the initial phase angle, in degrees. | Standard error of that phase angle. |
| `amp_sph_k` | `float` | Peak-to-peak amplitude of the touching-spheres model (mag). | $\Delta m = -2.5 \log_{10}(A_{\min} / A_{\max})$. |
| `m_mid_sph_k` | `float` | Mid-level reference magnitude of the touching-spheres model. | Reference magnitude of the fitted light curve. |
| `m_sep_k` | `float` | Magnitude at maximum projected separation, when the lobes do not occult each other. | Magnitude of the maximum projected area. |
| `p_rot_sph_h_k` | `float` | Rotation period of the touching-spheres fit, in hours. | Twice the photometric period of the fit. |

### B. Equator-on triaxial ellipsoid

The asteroid is modeled as a triaxial ellipsoid with semi-axes $a \ge b \ge c$, viewed in the equatorial plane.

| Column | Type | Description | Calculation Method |
| :--- | :--- | :--- | :--- |
| `chi_r_ell_k` | `float` | Reduced chi-squared ($\chi^2/\nu$) of the ellipsoid fit on channel $L_k$. | $\chi^2 / (N - 3)$ for the axis ratio $b/a$, the initial phase, and the reference magnitude. |
| `b_over_a_k` | `float` | Equatorial axis ratio ($b/a \le 1$). | Elongation that produces the double-peaked light curve. |
| `b_over_a_err_k` | `float` | Formal 1-sigma uncertainty of $b/a$. | Taken from the covariance of the fit. |
| `amp_ell_k` | `float` | Peak-to-peak amplitude of the ellipsoid model (mag). | $\Delta m = 2.5 \log_{10}(a / b)$. |
| `m_mid_ell_k` | `float` | Mean reference magnitude of the ellipsoid model. | Reference magnitude of the fitted light curve. |
| `m_ab_k` | `float` | Magnitude at maximum projected area ($\pi a b$). | Magnitude of the brightest orientation. |
| `p_rot_ell_h_k` | `float` | Rotation period of the ellipsoid fit, in hours. | Twice the photometric period of the fit. |

---

## 4. Single-Band Results

These columns use the measurements in one filter. That filter is `band_periodogram`. It is $r$ (`Lr`) unless the catalogue asks for another band. The request `Lx` means the filter with the most measurements, excluding $u$ and $y$; $r$ is kept when it ties with another filter.

### A. Periodogram

- `band_periodogram` (`str`): Filter used for the periodogram, the phase curve, the phase model, the shape fit, and the nightly time model (for example, `Lr` or `Li`).
- `periods_h` (`str`): Candidate rotation periods in hours, separated by commas, with Lomb-Scargle power of at least 80% of the strongest peak. The strongest period is listed first.
- `power` (`float`): Lomb-Scargle power of the period with the highest power, the first value in `periods_h`.
- `amp_lsm` (`float`): Peak-to-peak amplitude of the Lomb-Scargle model in `band_periodogram` ($\max - \min$).
- `amp_photo` (`float`): Photometric amplitude in `band_periodogram`, the 95th minus the 5th percentile of the magnitudes ($P_{95} - P_{5}$).

### A2. Voted period

The phase curve, the phase model, the shape fit, and the nightly time model are folded at this period.

- `p_rot_dem_h` (`float`): Rotation period, in hours, chosen by a vote. The voters are the Lomb-Scargle maxima of the filters with the most measurements, excluding $u$ and $y$, and one joint Lomb-Scargle of those same filters. The default vote uses two filters. Periods that agree to 3%, or whose photometric frequencies differ by one or two cycles per day, count as the same vote. The reported period is the joint Lomb-Scargle period when every series votes for it. Otherwise it is the Lomb-Scargle period of the filter with the most measurements among the series that cast the winning vote.
- `votes_dem` (`int`): Number of series that voted for `p_rot_dem_h`.
- `n_series_dem` (`int`): Number of series in the vote.
- `bands_dem` (`str`): Filters that voted, for example `Lr+Lg`.
- `periods_dem_h` (`str`): Period adopted by each series, written as `Lr=6.72,Lg=6.70,Lr+Lg=6.72`.

### B. Shape fit

- Touching spheres: `chi_r_sph`, `q_sph`, `q_sph_err`, `phi_sph_deg`, `phi_sph_err_deg`, `amp_sph`, `m_mid_sph`, `m_sep`, `p_rot_sph_h`.
- Triaxial ellipsoid: `chi_r_ell`, `b_over_a`, `b_over_a_err`, `amp_ell`, `m_mid_ell`, `m_ab`, `p_rot_ell_h`.

The definitions match the columns with the `_k` suffix. These values are computed from the measurements in `band_periodogram`.

---

## 5. Reflectance spectrum

Same-night colors relative to the $r$ filter. Measurements with a magnitude uncertainty above 0.08 mag are excluded. The color is the weighted mean of $m - m_r$ over the nights when both filters were observed. The reflectance is normalized to the Sun and to the $r$ band, where $R(r) = 1$.

For the $g$, $i$, $z$, and $y$ filters, written below as $b$:

- `{b}-r_spec` (`float`): Mean color $m_b - m_r$, in magnitudes.
- `{b}-r_spec_err` (`float`): Uncertainty of that color, in magnitudes.
- `R_{b}` (`float`): Solar-normalized reflectance $R(b)/R(r)$.
- `R_{b}_err` (`float`): Uncertainty of the reflectance.
- `n_{b}_spec` (`int`): Number of nights that contribute to the color.
