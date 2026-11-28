"""Per-cell mean cortical conductance along the STN axis for the current default network (diagnostic only, no model change).

Rebuilds the default cell positions (tools.generate_locs with cfg settings) and the cortical arrivals exactly as
stn_inputs.add_cortical_input does (generate at rate x max scale, thin per cell to scale / max scale), for both
cortical states, input seeds 1-5, 3 s runs (t >= 0.5 s, as the metrics). Mean conductance = arrivals/s x weight per
arrival x time integral of the peak-normalised Exp2Syn 'glut' kernel (tau1 1 ms, tau2 5 ms).
Output: results/edge_artefacts/ctx_input.png, ctx_input.json
"""
import dataclasses
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
from cfg import cfg  # noqa: E402
import cortical_arrivals_v2 as C  # noqa: E402
from tools import generate_locs  # noqa: E402

OUT = os.path.join(ROOT, 'results', 'edge_artefacts')
DUR, T0, SEEDS = 3.0, 0.5, [1, 2, 3, 4, 5]
TAU1, TAU2 = 1.0, 5.0  # netParams synMechParams['glut']


def kernel_integral_ms():
    tp = TAU1 * TAU2 / (TAU2 - TAU1) * math.log(TAU2 / TAU1)
    f = 1 / (math.exp(-tp / TAU2) - math.exp(-tp / TAU1))
    return f * (TAU2 - TAU1)


def main():
    pvp, pvn = generate_locs(60, 0.2, seed=cfg.seeds['loc'], plot=False, pvp_placement=cfg.pvp_placement,
                             pvp_profile=getattr(cfg, 'pvp_profile', None), grid=getattr(cfg, 'loc_grid', 'original'))
    locs = pvp + pvn
    ax = np.clip([(c['xnorm'] + c['ynorm']) / 2 for c in locs], 0, 1)
    isP = np.arange(len(locs)) < len(pvp)
    a, b = cfg.ctx_topography
    sc = a + b * ax
    sc = sc / sc.mean()
    K = kernel_integral_ms()
    nAss = max(1, math.ceil(len(locs) / cfg.ctx_cells_per_assembly))
    res = dict(axis=ax.tolist(), is_pvp=isP.tolist(), rate_scale=sc.tolist(), kernel_integral_ms=K,
               ctx_topography=list(cfg.ctx_topography), ctx_overrides=cfg.ctx_overrides, states={})
    for state in ('healthy', 'parkinsonian'):
        p = dataclasses.replace(C.STATES[state], **cfg.ctx_overrides)
        rates = []
        for s in SEEDS:
            t, nid, _ = C.generate_arrivals(n_neurons=len(locs), duration_s=DUR, dt_s=cfg.ctx_dt,
                                            params=dataclasses.replace(p, rate_hz=p.rate_hz * sc.max()), n_assemblies=nAss, seed=s)
            keep = np.random.default_rng([s, 1]).random(len(t)) < (sc / sc.max())[nid]
            t, nid = t[keep], nid[keep]
            rates.append(np.bincount(nid[t >= T0], minlength=len(locs)) / (DUR - T0))
        r = np.mean(rates, 0)
        g = r * p.conductance_nS * K / 1000  # nS, time-averaged
        res['states'][state] = dict(params=dataclasses.asdict(p), arrivals_per_s=r.tolist(), mean_g_nS=g.tolist(),
                                    expected_mean_g_nS=(p.rate_hz * sc * p.conductance_nS * K / 1000).tolist())
    json.dump(res, open(os.path.join(OUT, 'ctx_input.json'), 'w'), indent=1)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    xs = np.linspace(0, 1, 101)
    for k, state in enumerate(('healthy', 'parkinsonian')):
        A = axs[k]
        g = np.array(res['states'][state]['mean_g_nS'])
        A.scatter(ax[isP], g[isP], s=30, color='tab:orange', label='PV+ cell (realised, seeds 1–5)', zorder=3)
        A.scatter(ax[~isP], g[~isP], s=18, color='tab:green', label='PV− cell (realised, seeds 1–5)', zorder=3)
        p = res['states'][state]['params']
        ref = (a + b * xs) / (a + b * ax).mean()
        A.plot(xs, p['rate_hz'] * ref * p['conductance_nS'] * K / 1000, 'k-', lw=1, label=f'expected: {p["rate_hz"]:.0f}/s × ({a} {b:+}·x)/mean × {p["conductance_nS"]} nS × {K:.2f} ms')
        A.set_title(f'{state}: beta mod {p["beta_modulation"]} (uniform), shared {p["shared_fraction"]}, '
                    f'bursts {p["burst_rate_hz"]}/s × {p["burst_duration_s"]} s', fontsize=8)
        A.set_xlabel('STN axis (0 ventromedial, 1 dorsolateral)')
        A.grid(alpha=0.3)
    axs[0].set_ylabel('mean cortical conductance per cell (nS)')
    axs[0].legend(fontsize=7, loc='upper right')
    fig.suptitle('Cortical→STN input, current default network (configuration D, loc seed 1, 58 cells)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'ctx_input.png'), dpi=150)
    g0 = np.array(res['states']['parkinsonian']['mean_g_nS'])
    print('scale range', sc.min().round(3), sc.max().round(3), 'axis range', ax.min().round(3), ax.max().round(3))
    for st in res['states']:
        g = np.array(res['states'][st]['mean_g_nS']); r = np.array(res['states'][st]['arrivals_per_s'])
        print(st, 'arrivals/s', r.min().round(1), r.mean().round(1), r.max().round(1), '| g nS', g.min().round(3), g.mean().round(3), g.max().round(3),
              '| ventral(<0.25)/dorsal(>0.75)', (g[ax < 0.25].mean() / g[ax > 0.75].mean()).round(2),
              '| PV+ mean', g[isP].mean().round(3), 'PV- mean', g[~isP].mean().round(3))


if __name__ == '__main__':
    main()
