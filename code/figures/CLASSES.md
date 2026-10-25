## Class labels

Treat these as **tags** (multi-label), not a partition: a fold can be e.g. `clean-fold`+`phase-incomplete`, or `scattered`+`contaminated`.

| # | Glyph | Label | One-line sense |
|---|---|---|---|
| 1 | **c** | **clean-fold** | Low scatter; points hug the 2-peak model. |
| 2 | **s** | **scattered** | Large errors or high dispersion; shape still readable. |
| 3 | **x** | **contaminated** | A few points lie far from the main fold. |
| 4 | **p** | **phase-incomplete** | Only a segment of the rotation is covered. |
| 5 | **a** | **low-amplitude** | Small peak-to-peak; rotation still evident (often fast). |
| 6 | **m** | **misaligned** | Color/zero-point glitches; rotation curve still clear. |
| 7 | **u** | **undersampled** | Few points and/or large gaps in phase. |
| 8 | **e** | **excluded** | Reclassified later; drop from the science set. |
