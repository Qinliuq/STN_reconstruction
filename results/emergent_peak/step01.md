# Emergent beta peak: Steps 0-1 (digitization and derived parameters), 2026-11-19

Source: Jeon et al. 2022, Cell Reports 38:110439. The PDF is not in the repository. Every figure in it is a single
embedded 300 ppi JPEG, extracted at native resolution with `pdfimages`, so no re-rendering was needed. The scripts are
in `tools/digitize_jeon2022/`. Digitized data and overlays of the traced points on the original panels are in
`digitized/`. The derived parameters are in `params/jeon2022_profiles.json`, committed before any Step 2 run.
Figure: `step01_profiles.png`.

## Step 0: digitization

| panel | method | resolution | reading error | result |
|---|---|---|---|---|
| 7H middle: fraction positive, *Gria1* / *Gabra1* / *Pvalb* | colour-mask line tracing per pixel column; y fitted on 6 ticks (max residual 0.0006); x from the y axis to the last data column | 0.0030 fraction/px, 0.0037 axis/px | ±0.003 (line), ±0.007 (axis) | *Gria1* 0.53 → 0.28; *Gabra1* 0.48 → 0.34; *Pvalb* 0.07 → peak 0.34 at ~0.75 → 0.26 |
| 7H SE bands | tinted-run extent around each line | as above | coarse where bands overlap | median half-width: *Gria1* 0.059, *Gabra1* 0.018, *Pvalb* 0.039 (reading error only; not used in the parameters) |
| 7I center: expression, PV+ (solid) / PV− (dotted) | colour masks; per column, upper and lower line of each colour; solid/dotted assigned at the single crossing (axis 0.43); y has no ticks, so A.U. relative to the arrow baseline | 0.0038 axis/px | ±1 px line, ±3 px baseline (2–6% of the curve height) | PV+ peaks at ~0.6–0.7; PV− peaks at ~0.1 |
| 7I right: log₂(*Gria1*:*Gabra1*) | box edges, red median and whisker ends from pixels; y fitted on 7 ticks (max residual 0.067) | 0.031 log₂/px | ±0.07 | PV− median 0.45 (IQR 0.16–0.87); PV+ median 0.02 (IQR −0.18–0.25) |
| 1E left: STN projection profiles (55 rows) | colormap inversion against the colourbar (mean RGB error 4.5); Gaussian + baseline fit per row | 184 px across the axis (0.0054 axis/px), 5.9 px per row | ±0.005 axis; LUT step 0.005 | σ median **0.219** (IQR 0.197–0.253; n = 55; median fit R² 0.96); centroids 0.40–0.94 |
| 1H right top: GPe input vs STN projection, axis 1 | ring matched filter (radius 5.25 px), with the in-panel label masked | 0.0053 box units/px | ±0.005 | 50 of 55 points (5 hidden by overlap); Spearman ρ 0.878 (paper 0.886) |
| 3C / S5 (optional, condition M) | not done | | | deferred; M is optional and doesn't count toward the verdict |

**Is 7I center density-weighted?** Partly. The PV+ curves track the 7H *Pvalb* fraction (r 0.74–0.86) and the PV−
curves track 1 − *Pvalb* (r 0.94–0.95). After dividing by that type density, 2.5–5.8-fold variation along the axis
remains. The curves are population-level, but they also carry per-cell expression and/or overall cell-density
gradients that cannot be separated from these panels. The parameters therefore use 7H for the gradients and only the
per-cell ratio from 7I right.

**Axis orientation check: consistent.** Jeon axis 1 runs ventromedial (0) to dorsolateral (1). The model axis is
(xnorm + ynorm)/2 with 1 = dorsal end. *Pvalb* rises toward the dorsolateral STN and the model's PV+ cells lie at
0.48–0.89, so both use the same 0–1 orientation.
- **Mismatch:** the model's "measured gradient" PV+ placement is much more concentrated than *Pvalb*.
  - Model: no PV+ cells below 0.45, and a local PV+ fraction of 0.75 at 0.55.
  - Data: the fraction rises smoothly from 0.07 to 0.34.
  - The configuration-D placement is a hand-shaped skew-normal, not a fit to Fig. 7H. See question 1 below.

## Step 1: derived parameters (`params/jeon2022_profiles.json`)

| quantity | value | compared with |
|---|---|---|
| *Gria1* profile | cubic fit, R² 0.976; linear 0.568 − 0.273·a (R² 0.965) | matches the cfg comment 0.55 − 0.3·a |
| Glutamate ventral (< 0.25) : dorsal (> 0.75) over model cells | *Gria1* conductance **1.56** | 0.55 − 0.3·a: 1.71; config D 0.65 − 0.5·a (a rate scale): **2.53** |
| *Gabra1* profile | cubic fit, R² 0.885; linear 0.516 − 0.143·a | config D: no gradient |
| GABA per cell, F (*Gabra1* + ratio) | ventral:dorsal 1.29; PV+ mean 1.15, PV− 0.96 | config D: PV+ 1.22 / PV− 0.94, no gradient |
| GABA, abs (*Gabra1* only) | ventral:dorsal 1.24; PV+ 0.96, PV− 1.01 | |
| GABA, noRec (flat + ratio) | ventral:dorsal 1.04; PV+ 1.19, PV− 0.95 | |
| Per-type ratio r_PV+ / r_PV− | 2^(0.45 − 0.02) = **1.35** | config D 1.13 / 0.87 = 1.30 |
| GPe→STN σ | **0.219** (IQR 0.197–0.253) | config D 0.10 |
| GPe→STN mapping | STN centroid = **0.413 + 0.568 · GPe position**; residual SD 0.072 | config D: evenly spaced homes 0.04–0.96 |

Model formula (GABA): g_i = P(a_i) · r_t / r̄(a_i), with r̄(a) = f(a) · r_PV+ + (1 − f(a)) · r_PV−. Here f is the
model's local PV+ fraction (Gaussian kernel σ 0.1), and each profile is normalized to mean 1 over the model's cells.
This is an approximation: the 7H fractions are population-level and already contain the PV mix, and dividing by r̄
removes that mix before the per-type factor is applied.

Assumptions:
- The fraction of receptor-positive cells (7H) is proportional to synaptic conductance.
- Fig. 1H has no tick labels, so the mapping rescales its STN axis to the Fig. 1E centroid range (0.40–0.94), assuming
  both panels show the same 55 injections, and the GPe axis to 0–1 over the injection range (84% of the GPe covered,
  Fig. 1B).
- The 1E widths come from published density estimates (smoothed, normalized per row), so σ is an upper bound for the
  axonal spread.

Consequences visible before any run:
- The data-derived GPe input reaches the ventral STN only through the tails of σ = 0.22, because no group's home
  position is below 0.41.
- The measured glutamate gradient is much shallower than configuration D's rate topography (1.56 vs 2.53).

The profiles are implemented as off-by-default cfg options (`cfg.ctx_gscale`, `cfg.gpe_gscale`, `gpe_params`
`topo_sigma` / `group_center_map`). A default PD run is bit-identical to the commit before them (`f871152`).

## Before Step 2

1. **PV+ placement in F.** Should "measured gradient" stay configuration D's skew-normal, or become a placement drawn
   from the digitized *Pvalb* fraction? The latter would need a new off-by-default `pvp_placement` option. The
   skew-normal is far more concentrated than the data.
2. **Condition M.** Fig. 3C / S5 not digitized yet (optional); I can add it if you want M run.
3. **Runtime.** About 7.5 min per 3 s run. 6 conditions × 2 states × 5 seeds = 60 runs ≈ 1 h 55 min at 4 in parallel
   (+20 min with M).
