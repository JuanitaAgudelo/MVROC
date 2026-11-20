# MVROC rotation-solution classification and score

Tags and the quality score for entries that end up in MVROC-SOL.
This note is for the pipeline / catalogue workflow; the paper only names the
labels and that a 1–9 score exists.

## Class labels

Treat these as **tags** (multi-label), not a partition: a fold can be e.g.
`clean-fold`+`phase-incomplete`, or `scattered`+`contaminated`.

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

Glyph string for a solution: sorted unique letters, e.g. `cap`, `su`, `c`.
`e` may stand alone or with other tags; once `e` is set the solution is out of
the science sample regardless of score.

## Quality score $Q$ (1 = best, 9 = worst)

$Q$ ranks how much we trust that the adopted period is a real rotation rather
than a sampling / alias / noise artefact. Lower is better.

- **1**: high confidence that the rotation is real.
- **9**: likely artefact, or needs a deep human re-review before use.
- Solutions tagged **`e`** are not used in science products; by convention
  store $Q = 9$ (or leave $Q$ empty and filter on `e`).

### Design rules

1. Start from a clean baseline and **add penalties** for every tag that weakens
   period security or fold quality.
2. Penalties are chosen so that one mild tag stays in the “usable” band
   ($Q \le 3$), two serious tags push into “review” ($Q \ge 5$), and three or
   more serious tags reach the artefact / deep-review end ($Q \ge 7$).
3. `clean-fold` (`c`) does **not** lower the score by itself; it is the default
   endorsement when no penalty tags apply. If `c` is absent, add a small
   “no clean endorsement” surcharge so untagged-looking noisy folds do not look
   as good as an explicit `c`.
4. Clip to $[1, 9]$.

### Formula

Let $T$ be the set of glyphs on the solution.

$$
Q =
\begin{cases}
9 & \text{if } \mathtt{e} \in T \\[0.4em]
\mathrm{clip}_{[1,9]}\!\left(
  1
  + \mathbf{1}[\mathtt{c}\notin T]
  + \displaystyle\sum_{t \in T \setminus \{\mathtt{c},\mathtt{e}\}} w_t
\right) & \text{otherwise}
\end{cases}
$$

| Glyph | Class | Weight $w_t$ | Role |
|---|---|---:|---|
| **a** | low-amplitude | 1 | Amplitude small; period often still fine. |
| **m** | misaligned | 1 | Zero-point / colour mess; shape still clear. |
| **s** | scattered | 2 | High scatter; shape readable but weaker. |
| **x** | contaminated | 2 | Outliers pull the harmonic / period. |
| **p** | phase-incomplete | 2 | Incomplete phase → period less locked. |
| **u** | undersampled | 3 | Undersampling / gaps → highest alias risk. |
| **c** | clean-fold | 0 | Endorsement only (see surcharge if missing). |
| **e** | excluded | — | Forces $Q = 9$ / drop from science set. |

The missing-`c` surcharge is $+1$ (the $\mathbf{1}[\mathtt{c}\notin T]$ term).

### Worked examples

| Tags | Calculation | $Q$ | Reading |
|---|---|---:|---|
| `c` | $1$ | **1** | Best: clean fold. |
| `ca` | $1 + 1$ | **2** | Clean, low amplitude. |
| `cm` | $1 + 1$ | **2** | Clean shape, colour glitch. |
| `cp` | $1 + 2$ | **3** | Clean but phase incomplete. |
| `cap` | $1 + 1 + 2$ | **4** | Usable with caveats. |
| `s` | $1 + 1_{\mathrm{no}\,c} + 2$ | **4** | Scattered, no clean tag. |
| `cs` | $1 + 2$ | **3** | Clean + scatter (human saw both). |
| `su` | $1 + 1 + 2 + 3$ | **7** | Noisy and undersampled → review. |
| `xu` | $1 + 1 + 2 + 3$ | **7** | Contaminated + undersampled. |
| `pxu` | $1 + 1 + 2 + 2 + 3$ | **9** | Deep review / artefact risk. |
| `e` or `ce` | forced | **9** | Excluded from science set. |

### Suggested use bands

| $Q$ | Band | Catalogue use |
|---|---|---|
| 1–2 | **A** — high confidence | Default science sample. |
| 3–4 | **B** — usable with caveats | Science OK if the tags are acceptable for the use case. |
| 5–6 | **C** — soft review | Prefer human check before population statistics. |
| 7–9 | **D** — deep review / exclude | Do not treat as secure rotations; `e` always here. |

When several solutions exist for one MOI object, keep the one with smallest $Q$;
ties: prefer the one that carries `c`, then higher $A/s$ / lower $\chi^2_\nu$ if
those metrics are stored.

## Optional continuous refinement (not required for $Q$)

Tags stay primary so a human can reproduce $Q$ from the glyph string alone.
If later you want a continuous prior before tagging, the same quantities already
in the analysis rows are natural ingredients (do **not** mix them into $Q$ until
the tag set is frozen):

| Symbol | Meaning | Better when |
|---|---|---|
| $\chi^2_\nu$ | Reduced chi-squared of the joint 2-peak fit | lower |
| $A/s$ | Model amplitude over residual scatter | higher |
| power | Lomb–Scargle power at the adopted frequency | higher |
| phase coverage | Fraction of rotational phase with points | higher |
| `span_alias` | Campaign-alias flag on the adopted period | 0 preferred |

A possible future hybrid is to keep the tag formula above and only use the
continuous metrics to **break ties** or to suggest tags (e.g. high scatter →
propose `s`), not to replace the human labels.

## Implementation sketch

```text
def quality_score(tags: set[str]) -> int:
    tags = {t.lower() for t in tags}
    if "e" in tags:
        return 9
    weights = {"a": 1, "m": 1, "s": 2, "x": 2, "p": 2, "u": 3}
    q = 1
    if "c" not in tags:
        q += 1
    for t, w in weights.items():
        if t in tags:
            q += w
    return int(min(9, max(1, q)))
```

Store on each SOL row: `tags` (glyph string) and `Q` (int 1–9).

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
