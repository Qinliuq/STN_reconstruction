"""Edge / coverage artefact check: readouts, table and figure (results/edge_artefacts).

Runs (variant b, input seeds 1-5, fixed network): R0 = drive_topography runs (reused) + R0_geo here; R1 = F1 (centred,
reflected position grid); R2 = F2 (GPe coverage normalization); R3 = F1 + F2; R4 (F1 + F2 + F3) held (operating point).
Conditions 1, 2, 3, 6, 2b of the drive topography + 'geo' (uniform depth, uniform PV+ placement). Healthy reference:
same fix set and same PV+ placement, seed-matched.

Per run x condition, pooled over seeds:
- pairwise cross-power: seed-mean PD - healthy map (as in drive_topography/analyze.py): argmax over positions with
  n_eff >= 5, n_eff at that bin, interior argmax (>= 0.15, descriptive), 0.65/0.25 ratio of the map and of absolute PD;
- PPC (network_metrics.beta_ppc pooled over the 5 seeds per cell): PD - healthy map, argmax over positions with PPC
  n_eff >= 5 (cells with a valid PPC), n_eff there, 0.65/0.25 ratio; fraction of cells significantly locked (PD);
- rates PV+ / PV-, burst fraction, change vs R0 of the same condition (limits 10% / 25%);
- delivered GPe input (arrivals x per-arrival weight per second, t >= t0; mean, and ventral x < 0.15 / dorsal x > 0.85
  means) where logged (R0 drive_topography runs predate the log).
The previous verdict criteria (drive_topography) are re-evaluated per run on the cross-power map.
    python results/edge_artefacts/analyze.py  ->  edge_artefacts.csv, verdict_by_run.csv, difference_maps.png
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
import network_metrics as NM  # noqa: E402

SEEDS = [1, 2, 3, 4, 5]
POS = [round(x, 2) for x in NM.PROFILE_POSITIONS]
MIN_NEFF = NM.PROFILE_MIN_CELLS
DT = os.path.join(ROOT, 'sweeps', 'out', 'drive_topography')
FU = os.path.join(ROOT, 'sweeps', 'out', 'emergent_peak_followup')
EA = os.path.join(ROOT, 'sweeps', 'out', 'edge_artefacts')
CONDS = ['1_uniform_control', '2_territory_adv', '3_pv_uniform', '6_pv_monotonic', '2b_territory_adv_pvalb', 'geo']
PLACE = {'1_uniform_control': 'gradient', '2_territory_adv': 'gradient', '3_pv_uniform': 'uniform',
         '6_pv_monotonic': 'monotonic', '2b_territory_adv_pvalb': 'pvalb', 'geo': 'uniform'}
RUNS = ['R0', 'R1', 'R2', 'R3', 'R4']
FIXES = {'R0': 'none', 'R1': 'F1', 'R2': 'F2', 'R3': 'F1+F2', 'R4': 'F1+F2+F3 (held)'}
H0 = {'gradient': (FU, 'Dref_healthy'), 'pvalb': (FU, 'FU1_healthy'), 'uniform': (DT, 'H_uniform_healthy'),
      'monotonic': (DT, 'H_monotonic_healthy')}


def pd_dir(run, cond):
    if run == 'R0':
        if cond == 'geo':
            return EA, 'R0_geo_pd'
        return DT, ('ctl_1_uniform_control_pd' if cond == '1_uniform_control' else f'b_{cond}_pd')
    return EA, f'{run}_{cond}_pd'


def h_dir(run, cond):
    return H0[PLACE[cond]] if run == 'R0' else (EA, f'{run}_H_{PLACE[cond]}_healthy')


def load(d, lab):
    out = []
    for s in SEEDS:
        f = os.path.join(d, f'{lab}_s{s}')
        if not os.path.exists(os.path.join(f, 'model_output_metrics.json')):
            return None
        out.append((json.load(open(os.path.join(f, 'model_output_metrics.json'))), dict(np.load(os.path.join(f, 'model_output_spikes.npz')))))
    return out


def argmax(prof, neff, lo=0.0):
    p = np.where((np.asarray(neff) >= MIN_NEFF) & (np.asarray(POS) >= lo), prof, -np.inf)
    i = int(np.argmax(p))
    return POS[i], float(neff[i])


def ratio(prof):
    a, b = prof[POS.index(0.65)], prof[POS.index(0.25)]
    return float(a / b) if b > 0 else np.nan


def delivered(sp):
    if 'gpe_arrivals_cell' not in sp[0][1]:
        return {}
    T = (float(sp[0][1]['t1']) - float(sp[0][1]['t0'])) / 1000
    v = np.mean([d['gpe_arrivals_cell'] * d['gpe_weight_cell_nS'] / T for _, d in sp], 0)
    ax = sp[0][1]['axis']
    return dict(gpe_in_nS_per_s=float(v.mean()), gpe_in_ventral=float(v[ax < 0.15].mean()), gpe_in_dorsal=float(v[ax > 0.85].mean()) if (ax > 0.85).any() else np.nan,
                gpe_groups_ventral=float(sp[0][1]['gpe_groups_cell'][ax < 0.15].mean()), gpe_groups_mean=float(sp[0][1]['gpe_groups_cell'].mean()))


def one(run, cond):
    pd_, h_ = pd_dir(run, cond), h_dir(run, cond)
    P, H = load(*pd_), load(*h_)
    if P is None or H is None:
        return None
    bp = np.mean([m['profile']['beta_pow'] for m, _ in P], 0) * 1e6
    bh = np.mean([m['profile']['beta_pow'] for m, _ in H], 0) * 1e6
    neff = np.mean([m['profile']['n_eff'] for m, _ in P], 0)
    diff = bp - bh
    qp, qh = NM.beta_ppc([d for _, d in P]), NM.beta_ppc([d for _, d in H])
    pp, ph = np.array(qp['profile_ppc']), np.nan_to_num(np.array(qh['profile_ppc']))
    pneff = np.array(qp['profile_ppc_n_eff'])
    pdiff = pp - ph
    mx, nx = argmax(diff, neff)
    mi, _ = argmax(diff, neff, 0.15)
    px, pn = argmax(pdiff, pneff)
    rate = lambda S, k: float(np.mean([m[k] for m, _ in S]))
    row = dict(run=run, fixes=FIXES[run], condition=cond, placement=PLACE[cond],
               xpow_diff_max=mx, xpow_neff_at_max=nx, xpow_neff_lt8=nx < 8, xpow_diff_interior_max=mi,
               xpow_diff_ratio=ratio(diff), xpow_abs_ratio=ratio(bp), xpow_diff_max_value=float(diff.max()),
               ppc_diff_max=px, ppc_neff_at_max=pn, ppc_neff_lt8=pn < 8, ppc_diff_ratio=ratio(pdiff), ppc_abs_ratio=ratio(pp),
               ppc_abs_max=argmax(pp, pneff)[0], ppc_frac_locked=qp['frac_locked'], ppc_frac_locked_PVP=qp['frac_locked_PVP'],
               ppc_frac_locked_PVN=qp['frac_locked_PVN'], ppc_PVP=qp['ppc_PVP'], ppc_PVN=qp['ppc_PVN'], ppc_healthy_frac_locked=qh['frac_locked'],
               rate_PVP=rate(P, 'rate_PVP_pop'), rate_PVN=rate(P, 'rate_PVN_pop'), bursts=rate(P, 'frac_spikes_in_bursts'),
               rate_PVP_h=rate(H, 'rate_PVP_pop'), rate_PVN_h=rate(H, 'rate_PVN_pop'), bursts_h=rate(H, 'frac_spikes_in_bursts'),
               beta_pow=rate(P, 'beta_pow'),
               pd_axis_min=float(P[0][1]['axis'].min()), pd_axis_max=float(P[0][1]['axis'].max()), **delivered(P))
    prof = dict(diff=diff, pd=bp, h=bh, neff=neff, pdiff=pdiff, ppd=pp, pneff=pneff)
    return row, prof


def main():
    rows, profs = [], {}
    for run in RUNS:
        for cond in CONDS:
            r = one(run, cond)
            if r is not None:
                rows.append(r[0])
                profs[(run, cond)] = r[1]
    df = pd.DataFrame(rows)
    ref = df[df.run == 'R0'].set_index('condition')
    for k in ('rate_PVP', 'rate_PVN', 'bursts', 'rate_PVP_h', 'rate_PVN_h', 'bursts_h'):
        df[f'{k}_change_vs_R0_pct'] = [100 * (r[k] / ref.loc[r.condition, k] - 1) if r.condition in ref.index else np.nan for _, r in df.iterrows()]
    lim = lambda r: (max(abs(r.rate_PVP_change_vs_R0_pct), abs(r.rate_PVN_change_vs_R0_pct), abs(r.rate_PVP_h_change_vs_R0_pct), abs(r.rate_PVN_h_change_vs_R0_pct)) > 10
                     or max(abs(r.bursts_change_vs_R0_pct), abs(r.bursts_h_change_vs_R0_pct)) > 25)
    df['operating_point_breach'] = [bool(lim(r)) if np.isfinite(r.rate_PVP_change_vs_R0_pct) else False for r in df.itertuples()]
    df.to_csv(os.path.join(HERE, 'edge_artefacts.csv'), index=False, float_format='%.4g')
    # previous drive-topography criteria per run (cross-power map)
    ver = []
    for run in RUNS:
        S = df[df.run == run].set_index('condition')
        if '2_territory_adv' not in S.index:
            continue
        g = lambda c, k: S.loc[c, k] if c in S.index else np.nan
        ver.append(dict(run=run, fixes=FIXES[run],
                        c1_cond2_max_in_055_080=bool(0.55 <= g('2_territory_adv', 'xpow_diff_max') <= 0.80),
                        c2_cond3_max_gt_085=bool(g('3_pv_uniform', 'xpow_diff_max') > 0.85),
                        c3_cond6_max_gt_085=bool(g('6_pv_monotonic', 'xpow_diff_max') > 0.85),
                        c4_cond2_ratio_gt_cond1=bool(g('2_territory_adv', 'xpow_diff_ratio') > g('1_uniform_control', 'xpow_diff_ratio')),
                        rule1_2b=bool(0.55 <= g('2b_territory_adv_pvalb', 'xpow_diff_max') <= 0.80),
                        rule4_2b=bool(g('2b_territory_adv_pvalb', 'xpow_diff_ratio') > g('1_uniform_control', 'xpow_diff_ratio')),
                        cond1_max=g('1_uniform_control', 'xpow_diff_max'), cond2_max=g('2_territory_adv', 'xpow_diff_max'),
                        cond3_max=g('3_pv_uniform', 'xpow_diff_max'), cond6_max=g('6_pv_monotonic', 'xpow_diff_max'),
                        cond2b_max=g('2b_territory_adv_pvalb', 'xpow_diff_max'), geo_max=g('geo', 'xpow_diff_max'),
                        geo_interior_max=g('geo', 'xpow_diff_interior_max')))
    pd.DataFrame(ver).to_csv(os.path.join(HERE, 'verdict_by_run.csv'), index=False)
    figure(profs)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(df[['run', 'condition', 'xpow_diff_max', 'xpow_neff_at_max', 'xpow_diff_interior_max', 'xpow_diff_ratio', 'xpow_abs_ratio',
                  'ppc_diff_max', 'ppc_neff_at_max', 'ppc_diff_ratio', 'ppc_abs_max', 'ppc_frac_locked', 'rate_PVP', 'rate_PVN', 'bursts',
                  'rate_PVP_change_vs_R0_pct', 'rate_PVN_change_vs_R0_pct', 'bursts_change_vs_R0_pct', 'operating_point_breach']].to_string())
        print(pd.DataFrame(ver).to_string())
    return df


def figure(profs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cols = ['1_uniform_control', '2_territory_adv', '3_pv_uniform', 'geo']
    fig, axs = plt.subplots(len(RUNS), len(cols), figsize=(15, 2.7 * len(RUNS)), sharex=True)
    for i, run in enumerate(RUNS):
        for j, cond in enumerate(cols):
            a = axs[i, j]
            if (run, cond) not in profs:
                a.text(0.5, 0.5, 'held (operating point)' if run == 'R4' else 'not run', ha='center', va='center', transform=a.transAxes, color='gray')
                a.set_yticks([])
            else:
                p = profs[(run, cond)]
                ok = p['neff'] >= 8
                a.plot(POS, p['diff'], '-', color='k', lw=1.8, label='cross-power PD − healthy (×1e-6)')
                a.plot(np.array(POS)[~ok], p['diff'][~ok], 'o', color='tab:red', ms=3, label='n_eff < 8')
                b = a.twinx()
                b.plot(POS, p['pdiff'], '--', color='tab:blue', lw=1.5, label='PPC PD − healthy')
                b.tick_params(axis='y', colors='tab:blue', labelsize=7)
                if i == 0 and j == 0:
                    a.legend(fontsize=6, loc='upper left')
                    b.legend(fontsize=6, loc='upper right')
            if i == 0:
                a.set_title(cond, fontsize=9)
            if j == 0:
                a.set_ylabel(f'{run}: {FIXES[run]}', fontsize=9)
            a.axvspan(0.55, 0.80, color='0.92', zorder=-1)
            if i == len(RUNS) - 1:
                a.set_xlabel('STN axis (0 ventromedial, 1 dorsolateral)')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'difference_maps.png'), dpi=120)


if __name__ == '__main__':
    main()
