# Prompt for Claude Code — data-derived drive topography, six runs

Paste everything below the line into Claude Code, running in the repo that holds the
NetPyNE STN model. Copy `stn_drive_topography.py` and `stn_drive_runs_spec.md` into
that repo first (they are in this folder).

---

I want to replace the hand-fitted graded parkinsonian drive in this STN network model
with a drive topography derived from data, and run a six-condition comparison that can
falsify the hypothesis behind it. This is Priority 1 item 3 of the project to-do.

## Background you need

The model is a NetPyNE/NEURON spiking reconstruction of the sensorimotor STN,
currently "configuration D": 58 STN cells, two cell types (PV+ / PV−), equal Kir,
parkinsonian GPe→STN potentiation scale 2.0, PV+/PV− GABA:glutamate receptor ratio.
Beta is treated as arising upstream and imposed as rhythmic pallidal inhibition — an
emergent GPe–STN loop is explicitly out of scope, so do not add one.

An earlier 2×2 (drive topography × PV+ placement) found that the parkinsonian beta
maximum sits at axis 0.60–0.75 only when a graded drive and the PV+ cluster coincide;
the PV+ gradient alone did not localise beta. The weakness of that result is that the
graded drive's peak was fitted by hand at axis 0.65 — exactly the answer it was meant
to explain. The new drive profile removes that circularity: its shape comes from a
human STN functional-territory parcellation (Plantinga et al. 2018, 7T, 34 STNs)
combined with per-channel striatal dopamine-loss weights from the post-mortem and PET
literature, passed through a saturating transfer. Nothing in it is fitted to the model.

`stn_drive_topography.py` implements this. Read it first. The functions you need are
`drive_scaling(x, **kwargs)`, `pv_density(x, ...)`, and the `CONDITIONS` dict that
defines all six runs. `stn_drive_runs_spec.md` has the full specification, the
predicted input-profile behaviour, and the acceptance criteria.

## Step 1 — orient, and report back before changing anything

Find and report:
1. Where cell positions are assigned, and how a cell's position maps to a normalised
   coordinate on the STN principal axis.
2. **Which end of that coordinate is which.** The drive module assumes x = 0 is
   ventromedial/anteromedial (limbic end) and x = 1 is dorsolateral/posterolateral
   (motor end), matching Jeon's axis1 and the axis used in the existing results draft,
   where the beta maximum was at 0.60–0.75. If the model's coordinate runs the other
   way, say so and flip it in the adapter — do not silently proceed. This is the single
   most likely way to get a confidently inverted result.
3. The two hook points: the uniform parkinsonian GPe→STN potentiation (value 2.0) and
   the uniform pallidal beta depth (value 1.0).
4. How connectivity seeds and input seeds are set, and how the existing analysis bins
   cells along the axis for the spatial beta profile.

Stop and show me this before editing.

## Step 2 — wire it in

Add a thin adapter so each cell's parkinsonian change is scaled by its position:

```python
from stn_drive_topography import drive_scaling, pv_density, CONDITIONS, SEEDS
s = drive_scaling(x_cells, **CONDITIONS[cond]["drive"])
gpe_scale_cell  = 2.0 * s
beta_depth_cell = 1.0 * s
```

Two things to get right:

- **Renormalise over the realised cell positions, not the grid.** `drive_scaling`
  returns a profile with mean 1 over whatever `x` array it is given. With only 58 cells
  the realised positions are not uniformly spaced, so divide by the mean over the actual
  cell positions in that seed. Every condition must deliver the same *total* PD change
  as the uniform control; assert `abs(s.mean() - 1) < 1e-6` per run and fail loudly if not.
- **The healthy state does not change.** Keep healthy inputs exactly as measured,
  including the existing ventral cortical bias. The claim is about the PD − healthy
  contrast, not absolute beta. Do not retune, refit or rebalance anything to make the
  result come out.

## Step 3 — run

Six conditions from `CONDITIONS` × 5 connectivity seeds. Run variant (c), both hook
points scaled, first. If that completes cleanly, also run (a) GPe→STN potentiation only
and (b) beta depth only.

| # | key | drive | PV+ placement |
|---|---|---|---|
| 1 | `1_uniform_control` | uniform, same total | gradient |
| 2 | `2_territory_adv` | territory-weighted, advanced weights | gradient |
| 3 | `3_pv_uniform` | territory-weighted, advanced | uniform |
| 4 | `4_mouse_map` | mouse/Jeon territory intervals | gradient |
| 5 | `5_territory_early` | territory-weighted, early weights | gradient |
| 6 | `6_pv_monotonic` | territory-weighted, advanced | monotonic PV+ |

Each needs a matched healthy run with the same seed and the same PV+ placement, since
the readout is a difference map.

## Step 4 — readouts

Per run, then pooled over seeds, reusing the existing analysis code's binning rather
than inventing a new one:

- PD − healthy beta difference map along the axis; peak position per seed
- ratio 0.65/0.25 on both the difference map and absolute beta
- STN beta power, synchrony, fraction of spikes in bursts
- pairwise shared beta by cell type (PV+ pairs, PV− pairs, mixed) — the earlier result
  was a PV+ : PV− ratio of 2.4–3.4; check whether it survives a data-derived drive

Write one tidy CSV (one row per condition × seed) plus a summary table, and a figure
with the six difference-map profiles overlaid.

## Step 5 — verdict

The hypothesis survives if, across 5 seeds:
1. condition 2 puts the PD − healthy beta maximum at axis 0.55–0.80, **and**
2. condition 3 moves it dorsolateral of 0.85, **and**
3. condition 6 moves it dorsolateral of 0.85, **and**
4. condition 2's difference-map ratio exceeds condition 1's.

State plainly which criteria passed and which failed. If 1 holds but 3 and 6 do not,
the peak is being set by something else in the network — geometry, boundary effects,
connectivity density — and that needs identifying before anything is reported. Do not
reframe a failure as a success; a clean negative here is a usable result.

Two predictions worth checking in passing: conditions 2 and 4 should give nearly the
same peak position (their input profiles do), so a divergence means the network is
sensitive to the species assumption in the territory map rather than the map being
decisive. And condition 5 should sharpen the contrast without moving the peak.

## Constraints

- Do not add a GPe–STN loop, retune parameters, or change the healthy state.
- Do not compare the input-profile ratios in the spec against simulated beta ratios;
  the network compresses contrast heavily (a much steeper input gave a beta ratio of
  only 2.43). Peak *position* is the claim being tested.
- `PV_PEAK = 0.68` and `PV_SIGMA = 0.22` are assumptions — Jeon reports only
  "dorsolateral and middle", with no percentages. If time allows after the main runs,
  sweep peak 0.55–0.80 and sigma 0.15–0.30 and report how far the peak position moves.
- If a run fails or a readout looks degenerate (e.g. no beta in the PD state), stop and
  report rather than adjusting parameters to rescue it.
