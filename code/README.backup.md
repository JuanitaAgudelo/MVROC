# Medellin VRO Catalogue

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23486386.svg)](https://doi.org/10.5281/zenodo.23486386)

## Vera Rubin Observatory Catalogue of the Photometric and Dynamical properties of asteroids

The Medellin VRO Catalogue collects 524 asteroids measured in the public Rubin First Look data. Each row is one asteroid. It brings together the provisional designation, the time span and the number of measurements in each filter, the orbit at the first observation, mean absolute magnitudes, colors from the joint multiband light-curve fit, a rotation period chosen by a vote across filters, the parameters of that multiband fit on the synthetic channel Lk, two shape fits of the Lk light curve, a solar-normalized reflectance, and the physical properties reported by the JPL Small-Body Database, and the rotation periods of Greenstreet et al. (2026). A missing value is NaN.

Every public Rubin measurement of a solar-system object. The colour is the number of measurements in each cell, on a logarithmic scale. The white points are the 524 asteroids in this catalogue, one point each, at the mean position of its measurements. Most of them sit in one field near right ascension 186° and declination +5°, so they overlap at the right edge of the map. The dashed curves are the ecliptic and the galactic equator.

<img src="gallery/vro-observations-hexbin.png" alt="Density of Rubin solar-system observations" width="50%">

*Observations published by the Asteroid Institute and available publicly under [10.5281/zenodo.17047589](https://doi.org/10.5281/zenodo.17047589).*

This release is `releases/medellin-vroc-v19_10_2026.csv`. `medellin-vroc-latest.csv` is a symbolic link to that file. Git stores the table itself, not the link. The notebook `code/dev_build_catalogue.ipynb` writes the dated file in `releases/`, points the link at it, and also writes the same table to `code/medellin-vroc.csv`.

```python
import pandas as pd

catalogue = pd.read_csv("medellin-vroc-latest.csv")
catalogue.head()
```

## Scientific team

The catalogue was built by Team 5 at [Open Tools, Open Sky: Doing LSST Science with Public Rubin Data & LINCC Frameworks](https://plazasmalagon.com/lsst-science-public-rubin-data-products-hackathon-UdeA-2026/), a LATAM Hack Days workshop held in person at the Universidad de Antioquia, Medellín, Colombia, on 13–15 October 2026. 

- Jorge Iván Zuluaga Callejas, Universidad de Antioquia, Colombia
- Nicolás Pan, Universidad de la República, Uruguay
- Juanita Andrea Agudelo, Universidad de Antioquia, Colombia
- Valeria Judith Burgoa Camacho, Universidad Mayor de San Andrés, Bolivia
- Angie Melissa Gutiérrez Estupiñán, Universidad Industrial de Santander, Colombia
- Pedro Bernardinelli, mentor, Universidade de São Paulo, Brazil

<img src="gallery/Medellin-VROC-team.webp" alt="Medellin VRO Catalogue team" width="50%">

*From left to right: Jorge, Juanita, Vale, Nico, and Angie, with Pedro behind them. Behind the group are the city of Medellín, the Universidad de Antioquia station of the Medellín Metro, and the dome of the Planetario de Medellín observatory.*

## Example: Comparing MVROC Rotation Periods with Published Values

Below, we present a practical example demonstrating how to analyze and compare asteroid rotation periods using the Medellin VRO Catalogue. This example illustrates how to load the catalogue into a pandas DataFrame, extract relevant columns, filter data to cross-match MVROC periods with published values from other sources (such as the JPL Small-Body Database and Greenstreet et al. 2026), and visualize these comparisons graphically. The workflow can easily be adapted to other properties in the catalogue, providing a template for exploring and validating scientific results with the dataset.

The columns of the catalogue, and a summary of this sample, are described in [CATALOGUE.md](CATALOGUE.md).

Here, `p_rot_multi_h` is the period from the multiband fit, in hours. The black points use `rot_per_sbdb`, the period reported by the JPL Small-Body Database. The red points use `rot_per_lsm_gst`, the Lomb–Scargle period of Greenstreet et al. (2026). The legend gives the number of asteroids in each set. The dashed line is the double-period alias.

```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

catalogue_file = "medellin-vroc-latest.csv"
catalogue = pd.read_csv(catalogue_file)

# Extract periods already published in the JPL's Small-Body Database
sbdb = catalogue.dropna(subset=["p_rot_multi_h", "rot_per_sbdb"])
# Extract periods published in Greenstreet et al. (2026)
greenstreet = catalogue.dropna(subset=["p_rot_multi_h", "rot_per_lsm_gst"])

# Combine all periods
periods = np.concatenate([
    sbdb["rot_per_sbdb"], sbdb["p_rot_multi_h"],
    greenstreet["rot_per_lsm_gst"], greenstreet["p_rot_multi_h"],
])

# Plot periods
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(sbdb["rot_per_sbdb"], sbdb["p_rot_multi_h"], s=5, c="k", zorder=3,
           label=f"SBDB ({len(sbdb)})")
ax.scatter(greenstreet["rot_per_lsm_gst"], greenstreet["p_rot_multi_h"], s=5, c="r", marker="s", zorder=3,
           label=f"Greenstreet et al. (2026) ({len(greenstreet)})")

# Decoration
limits = np.array([periods.min() * 0.7, periods.max() * 1.4])
ax.plot(limits, limits, "k-", label="Equal period")
ax.plot(limits, 2 * limits, "k--", label="MVROC = 2 × published")
ax.plot(limits, 0.5 * limits, "k:", label="MVROC = 0.5 × published")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(limits)
ax.set_ylim(limits)
ax.set_aspect("equal")
ax.set_xlabel("Published period [h]")
ax.set_ylabel("MVROC period [h]")
ax.legend(loc="upper left")
ax.grid(alpha=0.3, which="both")
```

<img src="gallery/mvroc-sbdb-periods.png" alt="MVROC rotation periods compared with the SBDB and Greenstreet et al. (2026)" width="50%">

## Gallery

The figures below showcase asteroid **2023 WK24**, observed in g, r, and i bands with a total of 401 measurements to date. Its rotation period was determined by unanimous agreement across six period-finding algorithms. We highlight the best shape model obtained under stringent quality criteria: at least 250 measurements, full consensus among votes, ellipsoid amplitude between 0.25 and 1 mag, axis ratio `q` from 0.3 to 0.95, `b/a` from 0.45 to 0.85, multiband residuals under 0.05 mag, and harmonic amplitude within a factor of two of the ellipsoid amplitude. Within these constraints, the lowest reduced chi-squared for the ellipsoid model is 1.30, while the touching-spheres model yields 1.46.

### Light curve, all nights

<img src="gallery/2023_WK24_lightcurve_allnights_H.webp" alt="Light curve of 2023 WK24 on every night, in g, r, and i" width="50%">

Absolute magnitudes in g, r, and i. The top panel collects every night. Each panel below is one night, measured from its first observation.

### Multiband phase curve

<img src="gallery/2023_WK24_phase_allbands.webp" alt="Phase curve of 2023 WK24 in g, r, and i" width="50%">

Absolute magnitudes in g, r, and i, folded at the voted period of 3.33 h and shifted by the brightest point of each band.

### Multiband light-curve fit

<img src="gallery/2023_WK24_phase_model_Lk.webp" alt="Joint light-curve fit of 2023 WK24 on channel Lk" width="50%">

The same measurements, shifted onto Lk and fitted with four Fourier harmonics. The refined period is 3.328 h.

### Shape fit

<img src="gallery/2023_WK24_phase_shape_fit_Lk.webp" alt="Touching-spheres and ellipsoid fits to 2023 WK24" width="50%">

Touching spheres (`q` = 0.61, reduced chi-squared 1.46) and an equator-on ellipsoid (`b/a` = 0.71, reduced chi-squared 1.30) fitted to the Lk light curve. Lk is a band of this catalogue, not an LSST filter. It is the magnitude of the summed flux after each filter has been shifted onto one light curve by the fitted colors.

### Spectrum

<img src="gallery/2023_WK24_spectrum.webp" alt="Solar-normalized reflectance of 2023 WK24" width="50%">

Solar-normalized reflectance from the joint colors in g, r, and i, with `R(r) = 1`. The dashed curves are sketches of the C, S, and D taxonomic complexes.

The full set of light curves and phase curves is in the [SEAP-UdeA gallery](https://seap-udea.github.io/gallery/?repo=MVROC).

## Suggested papers

- [Kurlander et al. (2025)](https://doi.org/10.3847/1538-3881/add685), Predictions of the LSST Solar System Yield: Near-Earth Objects, Main Belt Asteroids, Jupiter Trojans, and Trans-Neptunian Objects
- [Murtagh et al. (2025)](https://doi.org/10.3847/1538-3881/ade1db), Predictions of the LSST Solar System Yield: Discovery Rates and Characterizations of Centaurs
- [Murtagh et al. (2026)](https://doi.org/10.3847/1538-3881/ae27bd), Predictions of the LSST Solar System Yield: Neptune Trojans
- [Merritt et al. (2025)](https://doi.org/10.3847/1538-3881/add3ec), Sorcha: A Solar System Survey Simulator for the Legacy Survey of Space and Time
- [Holman et al. (2018)](https://doi.org/10.3847/1538-3881/aad69a), HelioLinC: A Novel Approach to the Minor Planet Linking Problem
- [Greenstreet et al. (2026)](https://doi.org/10.3847/2041-8213/ae2a30), Lightcurves, Rotation Periods, and Colors for Vera C. Rubin Observatory's First Asteroid Discoveries
- [Vavilov et al. (2026)](https://arxiv.org/abs/2610.07584), Discovery of 2 km Ultra-fast Rotating Asteroid in Rubin DP2 Catalog via Varying Fourier Order Method
- [Showalter et al. (2021)](https://doi.org/10.1016/j.icarus.2020.114098), A statistical review of light curves and the prevalence of contact binaries in the Kuiper Belt
- [Chow et al. (2026)](https://doi.org/10.3847/1538-4357/ae4e23), Predictions of Imminent Earth Impactors Discovered by LSST
- [Cheng et al. (2026)](https://doi.org/10.3847/1538-3881/ae5936), Assessing the Vera C. Rubin Observatory's Ability to Discover Asteroid Impactors before They Collide with Earth
- [Dorsey et al. (2025)](https://doi.org/10.3847/PSJ/adf8ca), The Visibility of the Ōtautahi–Oxford Interstellar Object Population Model in LSST
- [Bernstein et al. (2025)](https://doi.org/10.3847/PSJ/ade3de), Determinations of Asteroid Masses Using Mutual Encounters Observed in the Legacy Survey of Space and Time
- LSST Predictions [AAS Focus issue](https://iopscience.iop.org/collections/PSJ-250528-01)

## Data access

- [B612 Foundation mirror of LSST data on the MPC](https://b612.ai/rubin-mpc-downloads/)
- Individual objects in the [Minor Planet Center](https://www.minorplanetcenter.net)

## Disclaimer

As of 19 October 2026, this repository has not been properly curated. The scripts and notebooks are not adequately documented, because the project was developed jointly during [Open Tools, Open Sky: Doing LSST Science with Public Rubin Data & LINCC Frameworks](https://plazasmalagon.com/lsst-science-public-rubin-data-products-hackathon-UdeA-2026/). As the final scientific report is prepared, the repository and the pipeline will be documented accordingly.

## Acknowledgements

The authors express their gratitude to the organizers of the 'Open Tools, Open Sky: Doing LSST Science with Public Rubin Data & LINCC Frameworks' workshop, especially Juan Carlos Muñoz-Cuartas, Andrés A. Plazas Malagón, and all participating mentors and staff. We thank the Universidad de Antioquia (Medellín, Colombia) and the Institute of Physics’ FACom research group for their hospitality. We acknowledge support from the LSST Discovery Alliance (LSST-DA), the LINCC Frameworks initiative funded by Schmidt Sciences, the Vera C. Rubin Observatory Project Office, and the Minor Planet Center (MPC).

## AI Disclosure

Portions of the code review, refactoring, inline documentation, debugging, and repository restructuring in this project were assisted by AI language models (used as coding assistants inside Cursor and related tools).

The human authors assert that all scientific ideas, the overall project conception, the design of the numerical and scientific experiments, the interpretation of the results, and the conclusions drawn from this catalogue are original contributions of the human authors. AI tools were used exclusively as coding and writing assistants — analogous to a spell-checker or a compiler — and bear no intellectual authorship over the scientific content of this work. AI models also assisted with translation from Spanish (the native language of several authors) into English, and with English spelling and grammar review.
