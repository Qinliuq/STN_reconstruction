# Step 1 diagnostics (no model changes)

## 1. Positions near the edges (fixed network, seeds[loc] = 1; 58 cells)

| placement | cells x < 0.01 | exact 0.00 (clipped) | of which PV+ | min x | max x | cells x > 0.95 |
|---|---|---|---|---|---|---|
| gradient | 1 | 1 | 0 | 0.000 | 0.918 | 0 |
| uniform | 1 | 0 | 0 | 0.005 | 0.974 | 1 |
| monotonic | 1 | 1 | 0 | 0.000 | 0.974 | 1 |
| pvalb | 1 | 1 | 0 | 0.000 | 0.918 | 0 |

Histogram, all cells (PV+ in brackets), 0.05 bins starting at 0:

| bin | gradient | uniform | monotonic | pvalb |
|---|---|---|---|---|
| 0.00–0.05 | 2 (0) | 3 (1) | 2 (0) | 2 (0) |
| 0.05–0.10 | 1 (0) | 2 (1) | 1 (0) | 1 (0) |
| 0.10–0.15 | 3 (0) | 2 (0) | 3 (0) | 5 (0) |
| 0.15–0.20 | 5 (0) | 1 (0) | 5 (0) | 3 (1) |
| 0.20–0.25 | 3 (0) | 7 (1) | 3 (0) | 2 (1) |
| 0.25–0.30 | 4 (0) | 3 (1) | 4 (0) | 5 (0) |
| 0.30–0.35 | 4 (0) | 5 (1) | 4 (0) | 3 (0) |
| 0.35–0.40 | 3 (0) | 1 (0) | 4 (2) | 3 (1) |
| 0.40–0.45 | 5 (0) | 6 (1) | 2 (0) | 6 (1) |
| 0.45–0.50 | 1 (1) | 3 (1) | 4 (0) | 3 (0) |
| 0.50–0.55 | 4 (4) | 2 (0) | 4 (1) | 2 (1) |
| 0.55–0.60 | 4 (2) | 4 (1) | 3 (1) | 4 (2) |
| 0.60–0.65 | 3 (0) | 4 (0) | 5 (2) | 3 (0) |
| 0.65–0.70 | 5 (2) | 3 (1) | 3 (1) | 5 (2) |
| 0.70–0.75 | 4 (2) | 5 (1) | 3 (1) | 4 (2) |
| 0.75–0.80 | 2 (0) | 1 (0) | 2 (1) | 2 (0) |
| 0.80–0.85 | 2 (0) | 5 (1) | 4 (2) | 2 (0) |
| 0.85–0.90 | 2 (1) | 0 (0) | 1 (0) | 2 (1) |
| 0.90–0.95 | 1 (0) | 0 (0) | 0 (0) | 1 (0) |
| 0.95–1.00 | 0 (0) | 1 (1) | 1 (1) | 0 (0) |

## 2. n_eff (sum of kernel weights, sigma 0.12) at each profile position

| x | gradient | uniform | monotonic | pvalb |
|---|---|---|---|---|
| 0.05 | 10.6 | 10.4 | 10.6 | 10.5 |
| 0.10 | 14.0 | 13.2 | 13.8 | 13.9 |
| 0.15 | 16.9 | 16.0 | 16.8 | 16.8 |
| 0.20 | 19.1 | 18.4 | 19.0 | 19.0 |
| 0.25 | 20.4 | 20.3 | 20.5 | 20.4 |
| 0.30 | 21.0 | 21.4 | 21.2 | 21.1 |
| 0.35 | 21.0 | 21.8 | 21.4 | 21.4 |
| 0.40 | 21.0 | 21.8 | 21.3 | 21.4 |
| 0.45 | 21.0 | 21.6 | 21.2 | 21.4 |
| 0.50 | 21.2 | 21.4 | 21.1 | 21.4 |
| 0.55 | 21.4 | 21.2 | 21.1 | 21.4 |
| 0.60 | 21.4 | 20.8 | 21.0 | 21.1 |
| 0.65 | 20.6 | 20.1 | 20.4 | 20.4 |
| 0.70 | 19.1 | 18.7 | 19.1 | 18.9 |
| 0.75 | 16.8 | 16.8 | 17.1 | 16.7 |
| 0.80 | 14.1 | 14.2 | 14.3 | 14.0 |
| 0.85 | 11.2 | 11.2 | 11.2 | 11.2 |
| 0.90 | 8.3 | 8.2 | 8.2 | 8.3 |
| 0.95 | 5.7 ⚠ | 5.5 ⚠ | 5.5 ⚠ | 5.7 ⚠ |

⚠ = n_eff < 8.

## 3. GPe coverage per STN cell (PD state, configuration D connectivity; mean over input seeds 1–5)

Total GPe weight = groups connected × per-arrival weight (1.0 nS × synapse scale 2.0 × type factor 1.13 PV+ / 0.87 PV−); all groups share the same baseline rate. Kernel profile (σ 0.12) over x:

| placement | groups/cell at x=0.05 | at 0.50 | at 0.95 | expected groups (Σ affinity) at axis 0 / 0.5 / 1 |
|---|---|---|---|---|
| gradient | 2.52 | 3.05 | 3.02 | 1.50 / 3.01 / 1.50 |
| uniform | 2.41 | 3.00 | 2.81 | 1.50 / 3.01 / 1.50 |
| monotonic | 2.53 | 3.02 | 2.75 | 1.50 / 3.01 / 1.50 |
| pvalb | 2.55 | 2.99 | 3.02 | 1.50 / 3.01 / 1.50 |

## 4. Delivered beta depth, conditions 1 and 2, variant b (depth only, M 0.8825, gradient placement)

Per cell: delivered depth = mean depth over its connected groups; drive amplitude = Σ over its groups of depth × per-arrival weight (= delivered depth × total weight). Kernel profiles (σ 0.12), mean over seeds; argmax over positions with n_eff ≥ 5.

| condition | delivered depth: argmax / range | depth × weight: argmax | PD − healthy map: argmax | interior (≥ 0.15) map argmax |
|---|---|---|---|---|
| 1_uniform_control | 0.30 / 0.88–0.88 | 0.50 | 0.05 | 0.55 |
| 2_territory_adv | 0.95 / 0.69–0.99 | 0.95 | 0.60 | 0.60 |

## 5. Rate and shared beta along x (seed means; ×10⁻⁶ for beta), added because pairwise cross-power scales with rate²

| x | 0.05 | 0.25 | 0.45 | 0.55 | 0.65 | 0.85 | 0.95 |
|---|---|---|---|---|---|---|---|
| cond 1 PD rate (Hz) | 20.5 | 18.8 | 21.5 | 22.7 | 21.6 | 17.2 | 15.7 |
| cond 1 PD shared beta | 15.6 | 12.9 | 15.4 | 14.9 | 13.3 | 8.7 | 7.2 |
| cond 1 PD beta / rate² (×10⁻³) | 37.0 | 36.3 | 33.6 | 29.0 | 28.6 | 29.4 | 29.0 |
| healthy (Dref) rate | 25.3 | 23.1 | 25.5 | 26.6 | 25.0 | 19.4 | 17.7 |
| healthy shared beta | 3.3 | 5.2 | 5.4 | 4.2 | 3.5 | 1.1 | 0.5 |

Observations (diagnostic, no criterion):
- Edge positions: clipping is minor (one cell at exactly 0.00 in three placements, none PV+). The real asymmetry is the
  grid: bin centres i/12 run 0–0.917, so the ventral end has a half bin centred on 0 and the dorsal end has no bin
  above 0.917 (realised max 0.918 for gradient / pvalb). n_eff ≥ 10 at 0.05 but 5.5–5.7 at 0.95 (< 8, flagged).
- GPe coverage: the two most ventral cells get 1.8 groups vs ~3 in the interior (truncated kernel; home 0.042 is the
  lowest). Kernel profile 2.5 at 0.05 vs 3.0 at 0.5. Dorsally the profile stays ~3.0 because no cell lies beyond
  0.92 (home 0.958 still covers them). Total GPe weight per cell additionally carries the PV+ type factor (1.13 vs
  0.87), so its profile peaks where the PV+ cluster sits (~0.5).
- Cond 1 (uniform depth): delivered depth is flat (0.88) by construction; depth × weight peaks at 0.50 and the
  interior map maximum is at 0.55: these line up, but the product's shape is set by GPe coverage × PV+ type factor.
  The 0.05 map maximum is NOT at the product's peak (the product is lowest there).
- Cond 2 (variant b): delivered depth rises monotonically to 0.95 (0.69–0.99); depth × weight peaks at 0.95; the map
  peaks at 0.60. The beta peak does not line up with the delivered drive.
- The 0.05 maximum of PD − healthy is driven by the healthy reference: healthy shared beta is hump-shaped (max
  ~0.35–0.40) and low at both ends, while PD shared beta under a uniform drive tracks rate (PD beta / rate² is within
  ±13% along x). PD − healthy is therefore largest where the healthy profile is lowest relative to PD, i.e. at the
  ventral edge (dorsally PD beta also falls with rate). The interior hump in PD under a uniform drive sits on the rate
  hump at 0.5–0.55 (the PV+ cluster: PV+ fire ~33 Hz vs PV− ~17 Hz).
