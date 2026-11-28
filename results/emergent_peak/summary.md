# Emergent beta peak from measured distributions (Jeon et al. 2022): Step 2, 2026-11-22

Question: does a parkinsonian beta peak at STN axis 0.6–0.7 emerge under a **uniform** pallidal beta drive when all
spatial structure comes from measured distributions? Steps 0–1 (digitization, derived parameters): `step01.md`.
Parameters: `params/jeon2022_profiles.json`, committed before any run (`11df1f3`). Criteria and analysis code:
`analyze.py`, committed before any result existed (`fb9a83f`). Nothing was tuned after the runs.

Setup:
- Configuration D unchanged; all new features are cfg options, off by default. A default PD run is bit-identical to
  the commit before them.
- 17 condition-states × 5 seeds × 3 s = 85 runs (`sweeps/emergent_peak.py`), all succeeded.
- Uniform drive = `beta_width 1e9` (full depth in every GPe group), as in Table 3.
- PV+ placement in F / noRec / abs / sig / ref / M: drawn from the digitized *Pvalb* profile (shape only; 12 / 58 PV+
  cells kept). Local PV+ fraction vs data: `placement_local_pv_fraction.csv` (within 0.04 everywhere).
- Approximation for the per-type GABA factor: g_i = P(a_i) · r_t / r̄(a_i), r̄(a) = f(a)·r_PV+ + (1 − f(a))·r_PV−, f
  = the model's local PV+ fraction. The 7H fractions are population-level; dividing by r̄ removes the PV mix before the
  per-type factor (r_PV+ / r_PV− = 1.35, Fig. 7I right) is applied.
- Assumption: fraction of receptor-positive cells (Fig. 7H) ∝ synaptic conductance.

Files: `conditions.csv` (one row per condition and state), `profiles.csv` (axis profiles, incl. PD − healthy),
`profiles.png` (top: shared beta, PD mean ± SD over seeds and healthy; bottom: synchrony; grey band 0.55–0.75).

## Verdicts

- **Emergent peak (F): NO.** Mean-profile maximum at 0.15 (criterion 0.55–0.75), 0 / 5 seeds peak in 0.50–0.80,
  0.65 / 0.25 ratio 0.61 (criterion ≥ 1.5). Beta falls monotonically from the ventral to the dorsal STN.
- **noPV necessary: not testable** (formally the rule "fails" — ratio 0.41, peak 0.20 — but F fails too, so removing
  PV+ structure cannot show necessity).
- **noRec necessary: not testable** (ratio 0.86, peak 0.30; F fails too). noRec is the flattest profile: the measured
  receptor gradients push beta *ventrally*, not toward 0.65.
- **abs vs F:** the per-type ratio roughly doubles shared beta (0.047 vs 0.025 total; PV+ pairs 24.7 vs 14.0 ×10⁻⁶) but
  does not localize it (ratio 0.61 vs 0.40, both peak ventrally).
- **sig vs ref (data GPe topography):** ratio 0.43 vs 0.82, peaks 0.15 vs 0.05. The data topography (σ 0.22, homes
  0.41–0.98) makes beta more ventral and lowers rates (see below); with configuration D's topography it is flat-to-ventral.
- **F vs F_cfgD (PV+ placement):** only the exaggerated placement gives a hump at the right place — mean peak 0.55,
  3 / 5 seeds in 0.50–0.80, ratio 1.09 — still well short of the criterion. The measured *Pvalb* placement gives none.
- **M (assumes motor-territory origin of beta; not part of the verdict): NO peak at any locus.** Graded drive centred at
  0.572 / 0.634 / 0.697 (width 0.219): peaks 0.15 / 0.05 / 0.15, ratios 0.54 / 0.57 / 0.79, 0–1 / 5 seeds in range.
  With the measured receptor gradients and GPe topography even an imposed, correctly placed drive does not produce a
  dorsal peak (configuration D, graded drive: ratio 2.43).

## Table per condition (PD; healthy rows in `conditions.csv`)

| condition | mean peak | per-seed peaks | 0.65/0.25 | seeds in [0.50, 0.80] | beta frac | shared beta PV+ / PV− (×10⁻⁶) | rates PV+ / PV− (Hz) | rate change vs D |
|---|---|---|---|---|---|---|---|---|
| F | 0.15 | 0.40 0.15 0.05 0.10 0.25 | 0.61 | 0 | 0.29 | 24.7 / 9.0 | 28.9 / 13.3 | **−16.7%** |
| F_cfgD | 0.55 | 0.05 0.65 0.50 0.10 0.50 | 1.09 | 3 | 0.27 | 19.4 / 6.8 | 25.4 / 13.7 | **−18.8%** |
| noPV | 0.20 | 0.05 0.10 0.15 0.10 0.25 | 0.41 | 0 | 0.28 | 17.7 / 9.0 | 31.0 / 12.8 | **−16.6%** |
| noRec | 0.30 | 0.50 0.25 0.05 0.30 0.05 | 0.86 | 1 | 0.29 | 23.0 / 8.1 | 28.8 / 12.7 | **−19.0%** |
| abs | 0.15 | 0.15 0.20 0.05 0.10 0.15 | 0.40 | 0 | 0.21 | 14.0 / 3.5 | 27.9 / 13.2 | **−18.3%** |
| sig | 0.15 | 0.25 0.10 0.05 0.35 0.10 | 0.43 | 0 | 0.28 | 20.7 / 7.2 | 29.1 / 14.3 | −12.4% |
| ref | 0.05 | 0.05 0.05 0.05 0.05 0.50 | 0.82 | 1 | 0.31 | 23.0 / 12.4 | 34.3 / 16.2 | +0.7% |
| M0 (0.572) | 0.15 | 0.25 0.20 0.05 0.05 0.35 | 0.54 | 0 | 0.24 | 14.7 / 5.2 | 28.2 / 13.1 | **−18.3%** |
| M1 (0.634) | 0.05 | 0.10 0.25 0.15 0.05 0.05 | 0.57 | 0 | 0.23 | 14.4 / 5.3 | 28.4 / 13.1 | **−17.9%** |
| M2 (0.697) | 0.15 | 0.10 0.35 0.30 0.10 0.60 | 0.79 | 1 | 0.26 | 16.3 / 6.8 | 28.8 / 13.3 | **−16.7%** |

Every PD condition is parkinsonian in the usual sense (beta, synchrony and bursting above its healthy reference; PD −
healthy contrast positive along the whole axis), so the failure is about *where* beta is, not *whether* there is beta.
The PD − healthy contrast is also largest ventrally in F (0.15) and noPV; in F_cfgD it peaks at 0.55.

## Rates (> 15% vs configuration D is reported, not corrected)

Population rate change vs configuration D (healthy 37.5 / 20.3 Hz, PD 32.9 / 16.4 Hz, cell-weighted):
- **Healthy: −19.5% to −26.7% in every condition except ref (−1.2%).**
- **PD: −16.6% to −19.0% in F, F_cfgD, noPV, noRec, abs and M0–M2;** sig −12.4%; ref +0.7%.
- The drop follows the data GPe→STN topography: ref (configuration D topography) matches configuration D, sig (same
  receptors, data topography) drops. GPe rates also fall (≈ 30–31 → 27 Hz), so the STN–GPe loop gain changes; this
  was not corrected. The peak verdicts are therefore obtained at ~17–25% lower STN rates than configuration D.

## Table 3 conclusion "the PV+ gradient is necessary"

That conclusion (graded drive: ratio 2.43 with the PV+ gradient vs 1.25 with uniform PV+) was obtained with the
exaggerated configuration-D placement (no PV+ cells below axis 0.35, local PV+ fraction up to 0.43 at 0.55; data
0.07 → 0.34). With the measured *Pvalb* placement:
- Under uniform drive, the measured placement does essentially nothing for localization: F 0.61 vs noPV 0.41, both
  peaking ventrally. Only the exaggerated cluster (F_cfgD, 1.09; Table 3's uniform-drive row 1.17) leaves a hump at 0.55.
- The graded-drive comparison itself (configuration D, graded drive, measured placement) was not run here, so the
  conclusion **is not tested in its original form**. What this sweep shows is that the amplifier effect requires a
  concentrated PV+ cluster; the measured *Pvalb* gradient is too shallow to provide one. "The PV+ gradient is
  necessary" survives only as "a PV+ cluster at least as concentrated as configuration D's is necessary", and the data
  do not show such a cluster. A 10-run follow-up (configuration D, graded drive, *Pvalb* placement, PD + healthy) would
  settle the original form; not run.

## Caveats

- Seed-to-seed SD of the shared-beta profile is 30–50% of the mean; per-seed peaks are scattered in most conditions.
  F's failure is not marginal (ratio 0.61 vs 1.5; 0 / 5 seeds).
- Peaks at the ventral edge sit where n_eff is lowest (≥ 5 is required).
- Observation, not a criterion: synchrony (bottom row of `profiles.png`) peaks at 0.55–0.65 in F_cfgD and sig, but at
  0.25–0.45 in F.
- The 1E widths are smoothed density estimates, so σ 0.22 is an upper bound on axonal spread; a narrower σ would
  partly restore configuration D's topography (ref).
- Glutamate conductance ∝ *Gria1*-positive fraction ignores NMDA, synapse number and cortical rate topography.
