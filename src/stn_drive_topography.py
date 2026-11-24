"""
Data-derived drive topography for the STN network PD state.
Priority 1, item 3 of the to-do: derive the drive topography from data, not from the PV+ peak.

Replaces the hand-fitted graded drive (peak at axis 0.65) with a profile built from
  (i)  the STN functional-territory map  T_c(x)        [Plantinga 2018, 7T, 34 STNs]
  (ii) striatal dopamine loss per channel w_c          [Kish 1988/2008, Kordower 2013, Pineda-Pardo 2022]
  (iii) a saturating loss -> circuit-change transfer   g(D)

Axis convention: x in [0,1] along the STN principal axis (Jeon axis1 / Plantinga's
posterolateral-anteromedial axis).  x = 0 anteromedial/ventromedial (limbic end),
x = 1 posterolateral/dorsolateral (motor end).  Matches the axis used in
results-section-draft-v2.md, where the parkinsonian beta maximum lay at 0.60-0.75.

Hook points in the network model (both are per-cell, driven by cell position):
  * GPe->STN potentiation          currently uniform  PD_GPE_SCALE = 2.0
  * pallidal beta depth            currently uniform  BETA_DEPTH   = 1.0
Both are multiplied by drive_scaling(x), which has mean 1 across the population, so
the total PD change is unchanged relative to the uniform control.

Changes from the version supplied with the task (results/drive_topography/, decisions of 2026-11-24):
  * SEEDS: input seeds 1-5 (ctx_seed = s, gpe_seed = 100 + s) on ONE fixed network (cfg.seeds['loc'] /
    ['conn']), as in Table 3. The earlier comment "connectivity seeds, as in the 2x2" was wrong.
  * pv_density() is an analytic illustration only (the predicted input profiles); it is NOT a model input.
    PV+ placement in the model: 'gradient' = configuration D's skew-normal (tools.generate_locs), 'uniform',
    'monotonic' = generate_locs 'profile' with pv_shape_monotonic(), 'pvalb' = 'profile' with the digitized
    Jeon 2022 Pvalb fraction (params/jeon2022_profiles.json); 12 / 58 PV+ throughout.
  * The model applies the drive to the CHANGE, not the total (stn_inputs.GPeFeedback, cfg.pd_drive):
      GPe->STN synapse scale of cell i = 1 + (2.0 - 1) * s(x_i), s renormalized to mean 1 over the realised cells;
      beta depth of GPe group g = base + (M - base) * s(c_g), base = healthy depth 0.1, s renormalized so that the
      connection-weighted mean depth received by the cells is M; M is one value per comparison set (the largest
      s in the set keeps every depth <= 1), control included.
  * CONDITIONS adds 2b_territory_adv_pvalb (measured Pvalb placement) and 5c_uniform_control_early (uniform
    control on condition 5's depth scale M).
"""
import numpy as np

# ---------------------------------------------------------------- territory map
# Parcel extents calibrated to Plantinga et al. 2018 NeuroImage 168:403-411
# (17 PD patients, 34 STNs, 7T tractography; share of total STN volume):
#   motor 55.3 +- 14.0 %, associative 55.6 +- 15.9 %, limbic 20.3 +- 16.3 %
#   overlap motor&associative 17.0 +- 8.9 %, associative&limbic 15.6 +- 13.2 %
# (percentages sum to >100 % because the parcels overlap; zones "gradually
#  transition into each other", so edges are smoothed rather than hard.)
TERRITORY_INTERVALS = {
    "limbic": (0.000, 0.203),
    "assoc":  (0.061, 0.617),
    "motor":  (0.447, 1.000),
}
EDGE_WIDTH = 0.06

# mouse alternative, from Jeon et al. 2022 Cell Reports (hyperdirect):
# MOp/MOs -> CENTRAL STN, SSp/SSs -> dorsolateral.  Sensorimotor treated as one
# channel spanning centre -> dorsolateral.  Use as a species sensitivity check.
TERRITORY_INTERVALS_MOUSE = {
    "limbic": (0.000, 0.230),
    "assoc":  (0.080, 0.560),
    "motor":  (0.330, 0.950),
}

# ------------------------------------------------------- denervation weights w_c
# fraction of striatal dopamine lost in each channel
DENERVATION = {
    # advanced / DBS-candidate stage: motor putamen effectively saturated
    #   Kish 1988 (caudal putamen <1 % remaining); Kish 2008 (putamen -98 %, caudate -80 %);
    #   Kordower 2013 (dorsal putamen markers virtually gone by ~4 y post-diagnosis)
    "advanced": {"limbic": 0.45, "assoc": 0.80, "motor": 0.965},
    # early / de novo: Pineda-Pardo 2022 (hand 63.6 %, foot 58.4 % more-affected side);
    #   PPMI (contralateral putamen SBR -68 %, mean striatal -45 %)
    "early":    {"limbic": 0.20, "assoc": 0.35, "motor": 0.60},
}

# ------------------------------------------------------------- PV+ density along x
# Jeon 2022 Fig. 4C/D, Fig. 5B: glutamatergic PV+ prevalent in dorsolateral AND
# middle STN; PV- ventromedial and middle.  No percentages are given in the paper,
# so the centre and width are assumptions -- sweep them.
PV_PEAK, PV_SIGMA = 0.68, 0.22


def _smooth_box(x, lo, hi, w=EDGE_WIDTH):
    sig = lambda t: 1.0 / (1.0 + np.exp(-t / w))
    return sig(x - lo) * sig(hi - x)


def territories(x, intervals=None, w=EDGE_WIDTH):
    """T_c(x) for each channel, normalised so sum_c T_c(x) = 1."""
    intervals = intervals or TERRITORY_INTERVALS
    T = {c: _smooth_box(np.asarray(x, float), lo, hi, w) for c, (lo, hi) in intervals.items()}
    S = sum(T.values())
    return {c: T[c] / S for c in T}


def loss_index(x, stage="advanced", intervals=None):
    """D(x) = sum_c w_c T_c(x): local dopamine-loss index at STN position x."""
    T = territories(x, intervals)
    w = DENERVATION[stage]
    return sum(w[c] * T[c] for c in T)


def drive_scaling(x, stage="advanced", n=2.0, K=0.70, intervals=None, uniform=False):
    """Per-position multiplier for the PD change. Mean 1 across x, so the total
    PD change matches the uniform control.  uniform=True returns the control."""
    x = np.asarray(x, float)
    if uniform:
        return np.ones_like(x)
    D = loss_index(x, stage, intervals)
    g = D ** n / (D ** n + K ** n)
    return g / g.mean()


def pv_shape_monotonic(x):
    """Shape of the monotonic PV+ control (sigmoid ramp to the dorsolateral edge); used as a placement profile."""
    return 1 / (1 + np.exp(-(np.asarray(x, float) - 0.5) / 0.12))


def pv_density(x, peak=PV_PEAK, sigma=PV_SIGMA, monotonic=False):
    """PV+ density along x, mean-normalised. monotonic=True is the control that
    removes the dorsolateral falloff. ANALYTIC ILLUSTRATION ONLY (predicted input profiles); not a model input."""
    x = np.asarray(x, float)
    d = 1 / (1 + np.exp(-(x - 0.5) / 0.12)) if monotonic else np.exp(-0.5 * ((x - peak) / sigma) ** 2)
    return d / d.mean()


# --------------------------------------------------------------- the six runs
# Each entry gives the arguments for drive_scaling plus how PV+ is placed.
# pv_placement: "gradient" = configuration D's skew-normal PV+ placement (the model default; hand-made, more
#                            concentrated than the measured Pvalb profile),
#               "uniform"  = PV+ 20% of every grid bin,
#               "monotonic"= PV+ share ramps to the dorsolateral edge (pv_shape_monotonic, control),
#               "pvalb"    = digitized Jeon 2022 Fig. 7H Pvalb fraction.
# depth_set: comparison set sharing one depth scale M (never compare across sets).
CONDITIONS = {
    "1_uniform_control":   dict(drive=dict(uniform=True),                      pv_placement="gradient", depth_set="advanced"),
    "2_territory_adv":     dict(drive=dict(stage="advanced"),                  pv_placement="gradient", depth_set="advanced"),
    "2b_territory_adv_pvalb": dict(drive=dict(stage="advanced"),               pv_placement="pvalb",    depth_set="advanced"),
    "3_pv_uniform":        dict(drive=dict(stage="advanced"),                  pv_placement="uniform",  depth_set="advanced"),
    "4_mouse_map":         dict(drive=dict(stage="advanced",
                                           intervals=TERRITORY_INTERVALS_MOUSE), pv_placement="gradient", depth_set="advanced"),
    "5_territory_early":   dict(drive=dict(stage="early"),                     pv_placement="gradient", depth_set="early"),
    "5c_uniform_control_early": dict(drive=dict(uniform=True),                 pv_placement="gradient", depth_set="early"),
    "6_pv_monotonic":      dict(drive=dict(stage="advanced"),                  pv_placement="monotonic", depth_set="advanced"),
}
SEEDS = [1, 2, 3, 4, 5]          # INPUT seeds (ctx_seed = s, gpe_seed = 100 + s) on one fixed network, as in Table 3


def profile_table(xs=None, stage="advanced", n=2.0, K=0.70):
    xs = np.arange(0.05, 1.0, 0.10) if xs is None else np.asarray(xs, float)
    T = territories(xs)
    D = loss_index(xs, stage)
    s = drive_scaling(xs, stage, n, K)
    pv = pv_density(xs)
    rows = np.c_[xs, T["limbic"], T["assoc"], T["motor"], D, s, pv, s * pv / (s * pv).mean()]
    return rows, "x,T_limbic,T_assoc,T_motor,D,s_drive,pv_density,product"


if __name__ == "__main__":
    rows, hdr = profile_table()
    print(hdr.replace(",", "  "))
    for r in rows:
        print("  ".join(f"{v:7.3f}" for v in r))
    xs = rows[:, 0]
    f = lambda p, t: np.interp(t, xs, p)
    for name, cond in CONDITIONS.items():
        s = drive_scaling(xs, **cond["drive"])
        if cond["pv_placement"] == "pvalb":
            continue  # measured placement: no analytic illustration
        pv = pv_density(xs, monotonic=(cond["pv_placement"] == "monotonic"))
        if cond["pv_placement"] == "uniform":
            pv = np.ones_like(xs)
        prod = s * pv
        prod /= prod.mean()
        print(f"{name:20s} drive peak x={xs[s.argmax()]:.2f}  "
              f"product peak x={xs[prod.argmax()]:.2f}  ratio .65/.25={f(prod,0.65)/f(prod,0.25):.2f}")
