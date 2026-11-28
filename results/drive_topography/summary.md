# Data-derived PD drive topography — results, 2026-11-24

Spec as implemented: `spec.md` (decisions A–F). Criteria and analysis (`analyze.py`) committed before any run.
110 new runs (variant c: 40 PD + 10 healthy; variants a, b: 60 PD), 3 s each, input seeds 1–5 on one fixed network.
**All runs succeeded; every PD run has beta (beta power 0.039–0.066); no readout was degenerate; nothing was retuned.**
Files: `drive_runs.csv` (variant × condition × seed), `drive_summary.csv`, `verdict.json`, `difference_maps.png`.

Variants: **c** = both hooks follow s(x) (GPe→STN synapse scale and GPe beta depth); **a** = synapse only (depth uniform
at M); **b** = depth only (synapse uniform at 2.0). Controls 1 and 5c are the same runs in all three variants.

## Verdict: the hypothesis does not survive (in any variant)

"Maximum" = peak of the seed-mean PD − healthy shared-beta map (positions with n_eff ≥ 5), as fixed beforehand.

| criterion | c (both) | a (synapse) | b (depth) |
|---|---|---|---|
| 1. cond 2 max in [0.55, 0.80] | **fail** (0.05) | **fail** (0.05) | pass (0.60) |
| 2. cond 3 max > 0.85 | **fail** (0.05) | **fail** (0.05) | **fail** (0.05) |
| 3. cond 6 max > 0.85 | **fail** (0.70) | **fail** (0.05) | **fail** (0.60) |
| 4. cond 2 0.65/0.25 ratio > cond 1 (1.27) | pass (1.35) | **fail** (1.21) | pass (2.11) |
| **survives** | **no** | **no** | **no** |

Criteria 2 and 3 fail in every variant: with PV+ cells uniform (3) or ramped dorsally (6), the beta contrast does not
move dorsal of 0.85. Under the territory drive the PD-specific beta sits at ~0.55–0.70 whatever the PV+ placement does
(or at the ventral edge, below); it does not follow the drive × PV+ conjunction into the dorsolateral tip.

## The ventral-edge maximum (read before the tables)

Every condition, **including both uniform controls**, has a large PD − healthy value at the ventral edge (axis 0.05:
PD ≈ 13–19, healthy ≈ 2.5–5.6, ×10⁻⁶). It is present under a uniform drive, so it is a property of this network in the
PD state (one fixed wiring; decision A), not of the territory drive. Many "0.05" maxima above are this edge winning
narrowly over an interior hump. Descriptively (not a criterion), the interior local maximum of the seed-mean map:

| condition | c: edge / interior max | a: edge / interior max | b: edge / interior max |
|---|---|---|---|
| 1 uniform control | 12.3 / 0.55: 10.7 | (same) | (same) |
| 2 territory, skew-normal PV+ | 14.1 / **0.60: 13.6** | 13.3 / 0.55: 13.0 | 10.2 / **0.60: 12.8** |
| 2b territory, measured *Pvalb* | 13.0 / 0.50: 10.9 | 16.2 / 0.50: 11.1 | 13.7 / 0.50: 8.9 |
| 3 uniform PV+ | 12.7 / none | 13.0 / none | 14.6 / 0.30–0.40: 14.3 |
| 4 mouse map | 10.7 / 0.60: 10.9 | 13.7 / 0.55: 13.5 | 11.8 / **0.55: 16.7** |
| 5 early | 9.6 / 0.60: 10.7 | 12.5 / 0.60: 8.9 | 7.2 / 0.60: 11.3 |
| 5c uniform control, early | 10.6 / 0.60: 7.4 | (same) | (same) |
| 6 monotonic PV+ | 13.2 / 0.70: 14.2 | 18.0 / none | 7.6 / 0.60: 11.9 |

Even if the edge were excluded post hoc, criteria 2 and 3 would still fail (no interior maximum > 0.85 anywhere), so
the negative does not hinge on the edge. Whether the edge is wiring-specific is the question for the follow-up with
extra `seeds['loc']` / `seeds['conn']` draws.

## 2 vs 2b: the interior peak depends on the skew-normal approximation

| | cond 2 (skew-normal PV+) | cond 2b (measured *Pvalb*) |
|---|---|---|
| c: max / ratio | 0.05 (interior 0.60) / 1.35 | 0.05 (interior 0.50, barely above control) / **0.82** |
| a: max / ratio | 0.05 (interior 0.55) / 1.21 | 0.05 / 0.79 |
| b: max / ratio | **0.60 / 2.11** (criterion 1 met) | 0.05 / **0.87** |
| seeds with max in [0.55, 0.80] (c / a / b) | 1 / 2 / 2 | 0 / 0 / 1 |

**In the one variant where condition 2 meets criterion 1 (b), 2b does not, and 2b fails rule 4 in all three variants
(ratio < condition 1's 1.27).** By the rule set in decision E, the interior peak of condition 2 is an artefact of the
skew-normal PV+ placement: with the measured *Pvalb* fraction the territory drive gives no PD-specific contrast at
0.65 above the uniform control. (Contrast with FU1, where the *graded* configuration-D drive kept a peak at 0.70 with
the measured placement: the territory drive is broader and weaker at 0.65 than that hand-fitted drive.)

## Checks in passing

- **2 vs 4 (territory vs mouse map):** interior peak positions match (c 0.60 / 0.60; a 0.55 / 0.55; b 0.60 / 0.55).
  The mouse-map drive gives the stronger 0.65/0.25 contrast in a and the clearest peak in b (4 / 5 seeds in range).
- **5 vs 5c only (early weights; own uniform control at M 0.6757):** contrast gain over 5c is 1.31 (c), 0.49 (a),
  1.28 (b). The interior peak stays at 0.60 in 5 and in 5c in every variant, so where 5 sharpens contrast (c, b) it
  does so without moving the peak; in a it does not sharpen.
- **Hook decomposition:** the depth hook carries the effect. Distributing beta depth by s(x) (b) gives the largest
  interior contrast (cond 2 ratio 2.11, cond 4 interior 16.7); distributing only the synapse-scale change (a) gives
  nothing above control (cond 2 ratio 1.21 < 1.27); both together (c) sit in between (1.35). Rates in b stay within
  ~2% of configuration D.
- **PV+:PV− pairwise shared beta:** PV+ > PV− in 95 / 100 distinct PD seed-runs; ratio of seed means 1.49–3.00 per
  condition (cond 2: 2.51 / 2.24 / 2.69 in c / a / b). Consistent with the earlier result (1.85–3.96) but at the lower
  end, and below the 2.4–3.4 range in 2b, 3 and 6 (1.49–2.27). The direction is robust; the size depends on PV+
  placement.
- **Rates:** all conditions within 3% of configuration D (`rate_change_vs_cfgD_pct`).
- PV peak / sigma sweep: not run (the hypothesis already fails on criteria 2 and 3 in every variant).

## Caveats

- **Input seeds on one wiring** (decision A): every peak position is conditional on `seeds['loc'] = 1`,
  `seeds['conn'] = 4`. No wiring-robustness claim.
- **Depth scale M** (decision C): the advanced set runs at depth 0.8825, a ~12% weaker total beta drive than
  configuration D, so absolute beta is not comparable with Table 3. Condition 1 is a new run at M 0.8825, not
  `Dref_uni_pd`. Condition 5 (M 0.6757) is compared only with 5c.
- The ventral-edge maximum is outside the drive and PV+ manipulations; it is why per-seed peaks are unstable (as for
  configuration D, `Dref`: 3 / 5 seeds).
- These are post-spec facts only; no criterion, parameter or binning was changed after the data.

## One-line result

A data-derived territory drive does not reproduce the interior-PD-beta-plus-PV+-dependence pattern: dorsal PV+
placements do not move beta dorsally (criteria 2, 3 fail in all variants), condition 2's interior peak (0.60, met only
with the depth hook alone) disappears with the measured *Pvalb* placement, and the drive's effect is carried by GPe
beta depth, not by the GPe→STN synapse scale.
