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

### Photometric taxonomy

`gri_slope` is the slope of the solar-normalized reflectance between the LSST g and i bandpasses (483 nm and 755 nm), in percent per 100 nm. Reflectance is \(R(b)/R(r) = 10^{-0.4[(m_b-m_r)-(m_b-m_r)_\odot]}\), with solar colors g−r = 0.45 and i−r = −0.12. The i−r color is \(-(r-i)\).

`photo_spec_T` is the nearest of three sketch templates used on the reflectance plots: C (slope 1.25), S (14.6), and D (22.5), assigned at the midpoints 7.9 and 18.6. It is a photometric complex label, not a published Tholen class. Both columns are empty when `g-r` or `r-i` is missing.

| Column | Unit | Meaning |
|---|---|---|
| `gri_slope` | % / 100 nm | Reflectance slope from g to i. |
| `photo_spec_T` | | `C`, `S`, or `D`. |

### Rotation: multiband fit on Lk

Lk is the magnitude of the summed flux after each filter has been shifted onto one light curve by the fitted colors. Magnitudes in this section and in the shape fits are absolute magnitudes on Lk. The mean values `H_u` … `H_y` stay in the LSST bands.

| Column | Unit | Meaning |
|---|---|---|
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

Population and Tisserand parameters from the osculating elements `a`, `e`, `incl`, and `q`. They are not proper elements, and `flag_resonant` is a proximity flag, not a libration test.

The Tisserand parameter with respect to a planet of semimajor axis \(a_p\) is

\[
T = \frac{a_p}{a} + 2\cos i\sqrt{\frac{a}{a_p}(1-e^2)}.
\]

The planet distances used here are Earth 1.000 au, Mars 1.524 au, Jupiter 5.204 au, and Neptune 30.07 au. Inclination is the ecliptic inclination `incl`.

`dyn_class` is the first matching label:

| Label | Condition |
|---|---|
| `TNO` | \(a \ge 30.1\) au |
| `Trojan` | \(4.95 \le a \le 5.45\) au |
| `Centaur` | \(5.45 < a < 30.1\) au |
| `Thule` | \(4.20 \le a < 4.50\) au |
| `Hilda` | \(3.70 \le a < 4.20\) au |
| `Cybele` | \(3.28 \le a < 3.70\) au |
| `NEA` | \(q < 1.3\) au |
| `Hungaria` | \(1.78 \le a < 2.06\) au, \(i \ge 16^\circ\), \(q \ge 1.3\) au |
| `Mars-crosser` | \(q < 1.666\) au and \(a < 2.06\) au |
| `InnerMB` | \(2.06 \le a < 2.50\) au |
| `MiddleMB` | \(2.50 \le a < 2.82\) au |
| `OuterMB` | \(2.82 \le a < 3.28\) au |
| `Other` | none of the above |

`flag_resonant` is 1 when `a` falls inside the half-width of a Jupiter mean-motion resonance, and 0 otherwise. `resonant_ratio` is that resonance written as asteroid revolutions to Jupiter revolutions (`3:1`, `3:2`, `1:1`, …). If several windows overlap, the closest one in units of its half-width is kept. The nominal semimajor axis is \(a_J (n_J/n)^{2/3}\) with \(a_J = 5.204\) au.

| Ratio | Half-width |
|---|---:|
| `4:1`, `3:1`, `5:2`, `7:3` | 0.03 au |
| `2:1` | 0.04 au |
| `4:3` | 0.08 au |
| `3:2` | 0.15 au |
| `1:1` | 0.25 au |

| Column | Unit | Meaning |
|---|---|---|
| `dyn_class` | | Population label from the table above. |
| `T_earth` | | Tisserand parameter with respect to Earth. |
| `T_mars` | | Tisserand parameter with respect to Mars. |
| `T_jup` | | Tisserand parameter with respect to Jupiter. |
| `T_nep` | | Tisserand parameter with respect to Neptune. |
| `flag_resonant` | | 1 if `a` lies in a Jupiter resonance window, otherwise 0. |
| `resonant_ratio` | | Resonance ratio, empty when `flag_resonant` is 0. |

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
