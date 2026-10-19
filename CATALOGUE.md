# Medellin VRO Catalogue

This document provides a comprehensive overview of the Medellin VRO Catalogue (MVROC). It begins by summarizing the overall properties of the dataset, including information about the sampling of asteroids, the time span and quantity of observations, and key statistical measures relevant to the catalogue. In addition, it contains a detailed description of each column present in the data, clarifying the meaning, measurement units, and role of every field. Users will find both a high-level summary and an in-depth column-by-column reference to facilitate the effective use and interpretation of MVROC data.

<!-- begin catalogue_props.md -->

## Sample

The catalogue contains 524 asteroids and 120 columns. The measurements run from 2025-04-22 to 2026-04-10. Each asteroid has between 50 and 765 measurements (median 370). The sample contains 203088 measurements in all.

243 asteroids are numbered. An orbit is filled for 524 asteroids and a multiband rotation period for 524. `H_r` is filled for 514 asteroids, `g-r` for 523, a JPL rotation period for 11, and a Greenstreet et al. (2026) period for 53. The vote is unanimous for 227 asteroids.

The mean in each band is the average number of measurements among asteroids with at least one measurement in that band. The last row is the mean of `n_obs`.

| Band | Asteroids with observations | Mean observations per asteroid |
|---|---:|---:|
| u | 206 | 42.7 |
| g | 523 | 105.8 |
| r | 524 | 138.0 |
| i | 524 | 113.2 |
| z | 95 | 74.9 |
| y | 35 | 6.7 |
| all | 524 | 387.6 |

| Quantity | Minimum | Median | Maximum |
|---|---:|---:|---:|
| Measurements | 50 | 370 | 765 |
| Semimajor axis [au] | 1.84 | 2.69 | 33.36 |
| Eccentricity | 0.009 | 0.130 | 0.384 |
| Inclination [deg] | 1.97 | 8.07 | 32.22 |
| $H_r$ [mag] | 6.99 | 17.77 | 21.27 |
| Rotation period [h] | 0.05 | 5.61 | 46.32 |

<!-- end catalogue_props.md -->

## Columns

### Identification

| Column | Unit | Meaning |
|---|---|---|
| `provid` | | Provisional MPC designation, `YYYY LLnn`. |
| `permid` | | Permanent number: the integer in parentheses in `mpc_designation`. NaN when the asteroid is not numbered. |
| `mpc_designation` | | MPC designation used for the geometry. Numbered objects look like `(4273) Dunhuang`; unnumbered ones keep the provisional designation. |

### Observations

The photometric time span and the number of measurements. These counts are the measurements in the Rubin data, one column per LSST filter.

| Column | Unit | Meaning |
|---|---|---|
| `epoch_utc_ini` | UTC | Time of the earliest measurement, `YYYY-MM-DD HH:MM:SS`. |
| `epoch_utc_end` | UTC | Time of the latest measurement. |
| `n_obs` | | Number of measurements in all filters. |
| `n_obs_u`, `n_obs_g`, `n_obs_r`, `n_obs_i`, `n_obs_z`, `n_obs_y` | | Number of measurements in that filter. |

### Orbit

Elements at `epoch_tdb`, taken from the geometry table at the first observation.

| Column | Unit | Meaning |
|---|---|---|
| `epoch_tdb` | TDB | Epoch of the elements, equal to the time of the first observation. |
| `a` | au | Semimajor axis. |
| `e` | | Eccentricity. |
| `q` | au | Perihelion distance. |
| `incl` | deg | Inclination to the ecliptic. |
| `node` | deg | Longitude of the ascending node. |
| `peri` | deg | Argument of perihelion. |
| `M` | deg | Mean anomaly. |

### Mean absolute magnitudes

Mean absolute magnitude in each LSST band, from `colors.parquet`. These means are not phase-corrected. The colors in the next section come from the joint light-curve fit, so `g-r` is not `H_g − H_r`.

| Column | Unit | Meaning |
|---|---|---|
| `H_u`, `H_g`, `H_r`, `H_i`, `H_z`, `H_y` | mag | Mean absolute magnitude in that LSST band. |

### Phase-corrected colors

Colors from the joint multiband fit. One rotation model is shared by the filters, so a color does not require the two filters on the same night. The reference filter is r. The u measurements are left out of the periodic model and compared with it afterwards. Each color has an `_err` column, the 1-sigma uncertainty from the fit.

| Column | Unit | Meaning |
|---|---|---|
| `u-g`, `g-r`, `r-i`, `i-z`, `z-y` | mag | Color, the difference of the two magnitude offsets. |
| `u-r` | mag | Magnitude offset of u relative to r. |

### Rotation: voted period

The period is chosen by a vote. The voters are the Lomb–Scargle maxima of the three filters with the most measurements, excluding u and y, plus a joint Lomb–Scargle of each consecutive pair and of all three together. Periods that agree to 3%, or whose photometric frequencies differ by one or two cycles per day, count as the same vote. A unanimous vote uses the joint period of all the filters in the vote. Otherwise the reported period is the Lomb–Scargle period of the filter with the most measurements among the series that cast the winning vote.

| Column | Unit | Meaning |
|---|---|---|
| `p_rot_dem_h` | h | Rotation period chosen by the vote. |
| `votes_dem` | | Number of series that voted for that period. |
| `n_series_dem` | | Number of series in the vote. |
| `bands_dem` | | Filters in the vote, for example `Lr+Lg+Li`. |
| `periods_dem_h` | h | Period of each series, written as `Lr=3.33,Lg=3.33,Lr+Lg=3.33`. |

### Rotation: multiband fit on Lk

Lk is the magnitude of the summed flux after each filter has been shifted onto one light curve by the fitted colors. Magnitudes in this section and in the shape fits are absolute magnitudes on Lk. The mean values `H_u` … `H_y` stay in the LSST bands.

The light curve is taken to be double-peaked, so the rotational period is twice the photometric period. The fit uses between 2 and 10 Fourier harmonics and keeps the order with the lowest Bayesian information criterion. `p_rot_multi_used_h` is the period that starts the fit, here the voted period. `p_rot_multi_h` is the period after the fit.

| Column | Unit | Meaning |
|---|---|---|
| `p_rot_multi_used_h` | h | Period used to start the multiband fit. |
| `p_rot_multi_h` | h | Rotational period from the fit, `24 / f_rot_multi`. |
| `f_rot_multi` | d⁻¹ | Rotational frequency. |
| `k_harm_multi` | | Number of Fourier harmonics in the adopted model. |
| `bic` | | Bayesian information criterion of that model. |
| `amp_lsm_k` | mag | Peak-to-peak amplitude of the harmonic model on Lk over one cycle. |
| `amp_photo_k` | mag | Width of the Lk magnitudes between the 1st and 99th percentiles, after clipping residuals at 3σ. |
| `rms_multi` | mag | Root-mean-square residual of the joint fit. |
| `n_outliers_multi` | | Number of measurements rejected by the clip. |

### Shape: touching spheres on Lk

Two spheres in contact, fitted to the folded Lk light curve. `q_sph_k` is reported with the smaller sphere second, so it is at most 1. The spin axis makes an angle `phi_sph_deg_k` with the line of sight: 90° is equator-on and 0° is pole-on.

| Column | Unit | Meaning |
|---|---|---|
| `q_sph_k` | | Radius ratio `r2/r1`. |
| `q_sph_err_k` | | Uncertainty of `q_sph_k`. |
| `phi_sph_deg_k` | deg | Angle between the spin axis and the line of sight. |
| `phi_sph_err_deg_k` | deg | Symmetric uncertainty of that angle. NaN when the fit only gives a lower limit. |
| `chi_r_sph_k` | | Reduced chi-squared of the fit. |
| `amp_sph_k` | mag | Peak-to-peak amplitude of the model. |
| `m_mid_sph_k` | mag | Midpoint of the model curve, halfway between its brightest and faintest magnitudes. |
| `m_sep_k` | mag | Magnitude when the two disks are fully separated. |
| `p_rot_sph_h_k` | h | Rotation period refined by the fit. |

### Shape: ellipsoid on Lk

Equator-on triaxial ellipsoid (`phi = 90°`) rotating about its short axis. The long axis is `a = 1`.

| Column | Unit | Meaning |
|---|---|---|
| `b_over_a_k` | | Axis ratio `b/a`. |
| `b_over_a_err_k` | | Uncertainty of `b/a`. |
| `amp_ell_k` | mag | Peak-to-peak amplitude of the model. |
| `m_mid_ell_k` | mag | Midpoint of the model curve. |
| `m_ab_k` | mag | Magnitude of the `a`–`b` face, the brightest aspect when seen equator-on. |
| `p_rot_ell_h_k` | h | Rotation period refined by the fit. |
| `chi_r_ell_k` | | Reduced chi-squared of the fit. |

### Reflectance

Solar-normalized reflectance from the joint light-curve colors, relative to r, where `R(r) = 1`. The solar colors subtracted from `m_b − m_r` are g−r = 0.45, i−r = −0.12, z−r = −0.16, and y−r = −0.18. The reflectance is `10^(−0.4 × excess)`.

For the g, i, z, and y filters, written below as `b`:

| Column | Unit | Meaning |
|---|---|---|
| `{b}-r_spec` | mag | Color `m_b − m_r`. For g this equals `g-r`. For i it equals `−(r−i)`. |
| `{b}-r_spec_err` | mag | Uncertainty of that color. |
| `R_b` | | Solar-normalized reflectance `R(b)/R(r)`. |
| `R_b_err` | | Uncertainty of the reflectance. |
| `n_b_spec` | | Number of measurements in filter `b`. |

The stored names are `g-r_spec`, `R_g`, `n_g_spec`, and the same pattern for i, z, and y.

### SBDB physical parameters

From the JPL Small-Body Database. The number is stored without its unit. Each quantity has an `_err` column; that uncertainty is NaN when SBDB does not give a numeric error. Taxonomic classes have no numeric error.

| Column | Unit | Meaning |
|---|---|---|
| `H_sbdb` | mag | Absolute magnitude. |
| `albedo_sbdb` | | Geometric albedo. |
| `diameter_sbdb` | km | Diameter. |
| `extent_sbdb` | km | Axis extents, kept as text `a x b` or `a x b x c`. |
| `pole_sbdb` | deg | Ecliptic pole, kept as text `longitude/latitude`. |
| `density_sbdb` | g cm⁻³ | Bulk density. |
| `rot_per_sbdb` | h | Rotation period in the literature. |
| `GM_sbdb` | km³ s⁻² | Gravitational parameter. |
| `BV_sbdb` | mag | Color B−V. |
| `UB_sbdb` | mag | Color U−B. |
| `spec_B_sbdb` | | Bus taxonomic class. |
| `spec_T_sbdb` | | Tholen taxonomic class. |

### Greenstreet et al. (2026)

Rotation periods from Greenstreet et al. (2026), taken from `data/asteroids-with-rotation.csv`. A catalogue asteroid matches a row of that table when `provid` equals `provid` or `alt_designation`. Both values are NaN when the asteroid is not in the table.

| Column | Unit | Meaning |
|---|---|---|
| `rot_per_lsm_gst` | h | Lomb–Scargle rotation period. |
| `rot_per_fr_gst` | h | Fourier rotation period. |

### Dynamics

Reserved for a mean-motion resonance. They are not computed yet, so every value is NaN.

| Column | Unit | Meaning |
|---|---|---|
| `res_body` | | Planet the asteroid is resonant with. |
| `res_flag` | | Whether a resonance is identified. |
| `res_type` | | Kind of resonance. |
| `res_kp` | | One integer of the resonance. |
| `res_k` | | The other integer of the resonance. |
| `res_orbit` | | Description of the resonant orbit. |
