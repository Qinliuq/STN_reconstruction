"""Step 1: derive model parameters from the digitized Jeon et al. 2022 profiles -> params/jeon2022_profiles.json.

Committed before any network run (results/emergent_peak). Inputs: results/emergent_peak/digitized/*.csv / *.json
(written by fig7h.py, fig7i.py, fig7i_right.py, fig1e.py, fig1h.py). No simulation output is read.

- Profiles (Fig. 7H): cubic polynomial fit of the fraction-positive curve along axis 1 (and a linear fit, reported),
  tabulated at 0, 0.05, ..., 1. The model normalizes each profile to mean 1 over its own cells
  (stn_inputs.profile_scale / gpe_cell_scale), so the mean synaptic drive is unchanged.
- Glutamate: cortical conductance factor proportional to the Gria1 profile (cfg.ctx_gscale), with the cortical rate
  topography off (cfg.ctx_topography = None).
- GABA: GPe->STN conductance factor proportional to the Gabra1 profile (cfg.gpe_gscale.profile).
- Per-type ratio (Fig. 7I right): r_PV- = 1, r_PV+ = 2^(median_PV- - median_PV+) of log2(Gria1:Gabra1), applied to the
  GABA factor as g_i = P(a_i) r_t / rbar(a_i), rbar(a) = f(a) r_PV+ + (1 - f(a)) r_PV-, f = the model's local PV+
  fraction (Gaussian kernel sigma 0.1). Approximation: Fig. 7H fractions are population-level and already include the
  PV mix; dividing by rbar removes the mix before the per-type factor is applied.
- GPe->STN topography (Fig. 1E): sigma = median of the per-injection Gaussian widths (IQR reported).
- Mapping (Fig. 1H + 1E): Fig. 1H has no tick labels; its STN axis is rescaled so that the range of the data spans
  the range of projection centroids in Fig. 1E (same 55 injections), and the GPe axis to 0..1 over the data range
  (injections cover 84% of the GPe, Fig. 1B). Linear fit -> STN centroid = intercept + slope * GPe position.
Assumption: fraction of receptor-positive cells (Fig. 7H) is taken as proportional to synaptic conductance.

    python tools/digitize_jeon2022/derive_params.py
"""
import datetime
import json
import os
import subprocess

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DIG = os.path.join(ROOT, 'results', 'emergent_peak', 'digitized')
GRID = np.round(np.arange(0, 1.0001, 0.05), 2)


def fit_profile(x, y, deg=3):
    c = np.polyfit(x, y, deg)
    yf = np.polyval(c, x)
    r2 = 1 - ((y - yf) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    lin = np.polyfit(x, y, 1)
    r2l = 1 - ((y - np.polyval(lin, x)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return dict(cubic=c.tolist(), cubic_r2=float(r2), linear=dict(slope=float(lin[0]), intercept=float(lin[1]), r2=float(r2l)),
                x=GRID.tolist(), y=np.polyval(c, GRID).round(4).tolist())


def main():
    h = pd.read_csv(os.path.join(DIG, 'fig7H_fraction_positive.csv')).dropna(subset=['fraction'])
    prof = {g: fit_profile(s.axis1.values, s.fraction.values) for g, s in h.groupby('gene')}
    box = json.load(open(os.path.join(DIG, 'fig7I_right_boxes.json')))['boxes']
    r_pvp = 2 ** (box['PV-']['median'] - box['PV+']['median'])
    e = pd.read_csv(os.path.join(DIG, 'fig1E_projection_profiles.csv'))
    s1 = pd.read_csv(os.path.join(DIG, 'fig1H_gpe_vs_stn_axis1.csv'))
    mu_lo, mu_hi = float(e.mu.min()), float(e.mu.max())
    u = (s1.gpe_input_box - s1.gpe_input_box.min()) / (s1.gpe_input_box.max() - s1.gpe_input_box.min())
    c = mu_lo + (s1.stn_projection_box - s1.stn_projection_box.min()) / (s1.stn_projection_box.max() - s1.stn_projection_box.min()) * (mu_hi - mu_lo)
    slope, intercept = np.polyfit(u, c, 1)
    resid = float(np.std(c - (intercept + slope * u), ddof=2))
    gria, gabra = prof['Gria1'], prof['Gabra1']
    git = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    out = dict(
        source='Jeon H, Lee H, Kwon D-H, et al. Topographic connectivity and cellular profiling reveal detailed input '
               'pathways and functionally distinct cell types in the subthalamic nucleus. Cell Reports 38:110439 (2022). '
               'doi:10.1016/j.celrep.2022.110439',
        created=datetime.date.today().isoformat(), git_at_creation=git,
        method='tools/digitize_jeon2022/ (native 300 ppi figure images extracted with pdfimages; colour-mask line tracing, '
               'colormap inversion, ring matched filter); digitized data and overlays in results/emergent_peak/digitized/',
        axis='Jeon STN axis 1 (ventromedial 0 -> dorsolateral 1) taken as the model axis (xnorm + ynorm)/2, 1 = dorsal end; '
             'checked: Pvalb fraction rises toward the model PV+ region (both 0-1 normalized)',
        assumption_conductance='fraction of receptor-positive cells (Fig. 7H) is proportional to synaptic conductance',
        profiles={g: p for g, p in prof.items()},
        ctx_gscale={'x': gria['x'], 'y': gria['y']},
        gpe_gscale=dict(profile={'x': gabra['x'], 'y': gabra['y']}, type_ratio={'PVP': round(float(r_pvp), 4), 'PVN': 1.0},
                        pv_frac_sigma=0.1),
        gpe_gscale_flat_with_ratio=dict(profile=None, type_ratio={'PVP': round(float(r_pvp), 4), 'PVN': 1.0}, pv_frac_sigma=0.1),
        gpe_gscale_abs=dict(profile={'x': gabra['x'], 'y': gabra['y']}, type_ratio=None),
        per_type_ratio=dict(log2_gria1_gabra1_median={'PV-': box['PV-']['median'], 'PV+': box['PV+']['median']},
                            IQR={'PV-': [box['PV-']['q1'], box['PV-']['q3']], 'PV+': [box['PV+']['q1'], box['PV+']['q3']]},
                            r_PVP=float(r_pvp), note='GABA factor of PV+ relative to PV- (config D uses 1.13 / 0.87 = 1.30)'),
        gpe_topography=dict(topo_sigma=round(float(e.sigma.median()), 4), sigma_IQR=[float(e.sigma.quantile(.25)), float(e.sigma.quantile(.75))],
                            sigma_n=int(len(e)), centroid_range=[mu_lo, mu_hi],
                            group_center_map=[round(float(slope), 4), round(float(intercept), 4)], mapping_resid_sd=resid,
                            note='group home position = intercept + slope * GPe position (0..1); config D keeps sigma 0.1 and '
                                 'evenly spaced homes (0.04-0.96)'),
    )
    os.makedirs(os.path.join(ROOT, 'params'), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, 'params', 'jeon2022_profiles.json'), 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('per_type_ratio', 'gpe_topography')}, indent=1))
    for g, p in prof.items():
        print(g, 'cubic R2 %.3f; linear %.3f %+.3f*a (R2 %.3f)' % (p['cubic_r2'], p['linear']['intercept'], p['linear']['slope'], p['linear']['r2']))


if __name__ == '__main__':
    main()
