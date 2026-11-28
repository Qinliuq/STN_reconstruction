"""
Data-derived PD drive topography: readouts and verdict (criteria fixed before any run).

    python results/drive_topography/analyze.py

Inputs: sweeps/out/drive_topography/<label>_s<seed>/model_output_metrics.json (sweeps/drive_topography_c.py, _ab.py);
healthy references per PV+ placement, same INPUT seed (fixed network; peak positions are conditional on this wiring):
gradient -> sweeps/out/emergent_peak_followup/Dref_healthy, pvalb -> .../FU1_healthy, uniform / monotonic ->
sweeps/out/drive_topography/H_<placement>_healthy.

Binning: the existing src/network_metrics.py profile (positions 0.05..0.95, Gaussian kernel 0.12 over cell positions,
pairwise beta cross-power; a position counts for a peak only where n_eff >= 5).

Per run (one row per variant x condition x seed, drive_runs.csv): PD - healthy beta difference profile and its peak,
absolute PD beta peak, 0.65/0.25 ratio of both, STN beta power, synchrony, fraction of spikes in bursts, pairwise shared
beta PV+ / PV- / mixed pairs, rates (and change vs configuration D).
Pooled (drive_summary.csv): seed-mean difference profile (peak, 0.65/0.25 ratio), per-seed peaks, seed-mean absolute
profile (peak, ratio), means of the scalars; PV+ / PV- pair-beta ratio as ratio of seed means, mean +- SD of per-seed
ratios and number of seeds with PV+ > PV-.

Verdict (variant c is the test; a and b are reported with the same rules as the hook decomposition). "Across 5 seeds"
= the seed-mean difference map; per-seed peaks are reported next to it.
  1. condition 2: difference-map maximum in [0.55, 0.80]
  2. condition 3: difference-map maximum dorsolateral of 0.85 (> 0.85, i.e. 0.90 or 0.95)
  3. condition 6: difference-map maximum > 0.85
  4. condition 2's difference-map 0.65/0.25 ratio > condition 1's
  The hypothesis survives only if 1-4 all pass. Also reported: 2b by rules 1 and 4 (measured Pvalb placement);
  2 vs 4 peak positions (species map); 5 vs 5c (own control, depth scale M 0.6757; never compared with 2/3/6 in
  magnitude): contrast gain ratio(5)/ratio(5c) next to ratio(2)/ratio(1), and 5's peak position next to 2's.
Input-profile ratios from the spec are not compared with simulated beta ratios.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(ROOT, 'sweeps', 'out', 'drive_topography')
FU = os.path.join(ROOT, 'sweeps', 'out', 'emergent_peak_followup')
sys.path.insert(0, os.path.join(ROOT, 'src'))
import stn_drive_topography as D  # noqa: E402

SEEDS = D.SEEDS
MIN_NEFF = 5
CFGD = {'healthy': (37.4867, 20.2591), 'pd': (32.8933, 16.4365)}
HEALTHY = {'gradient': (FU, 'Dref_healthy'), 'pvalb': (FU, 'FU1_healthy'),
           'uniform': (OUT, 'H_uniform_healthy'), 'monotonic': (OUT, 'H_monotonic_healthy')}
CONTROLS = ('1_uniform_control', '5c_uniform_control_early')
ADV = ['1_uniform_control', '2_territory_adv', '2b_territory_adv_pvalb', '3_pv_uniform', '4_mouse_map', '6_pv_monotonic']
EARLY = ['5c_uniform_control_early', '5_territory_early']


def label(variant, cond):
    return f'ctl_{cond}_pd' if cond in CONTROLS else f'{variant}_{cond}_pd'


def load(d, lab, seed):
    f = os.path.join(d, f'{lab}_s{seed}', 'model_output_metrics.json')
    return json.load(open(f)) if os.path.exists(f) else None


def peak(prof, neff, pos):
    p = np.where(np.asarray(neff) >= MIN_NEFF, prof, -np.inf)
    return float(pos[int(np.argmax(p))])


def ratio(prof, pos):
    a, b = prof[pos.index(0.65)], prof[pos.index(0.25)]
    return float(a / b) if b > 0 else np.nan


def runs_table():
    rows, profs = [], {}
    for variant in ('c', 'a', 'b'):
        for cond, c in D.CONDITIONS.items():
            hd, hl = HEALTHY[c['pv_placement']]
            for s in SEEDS:
                r, h = load(OUT, label(variant, cond), s), load(hd, hl, s)
                if r is None or h is None:
                    continue
                pos = [round(x, 2) for x in r['profile']['positions']]
                bp, bh = np.array(r['profile']['beta_pow']), np.array(h['profile']['beta_pow'])
                neff = r['profile']['n_eff']
                diff = bp - bh
                profs[(variant, cond, s)] = dict(pos=pos, diff=diff, pd=bp, h=bh, neff=neff, sync=r['profile']['sync'])
                ref = CFGD['pd']
                rows.append(dict(variant=variant, condition=cond, depth_set=c['depth_set'], pv_placement=c['pv_placement'], seed=s,
                                 diff_peak=peak(diff, neff, pos), diff_ratio_065_025=ratio(diff, pos),
                                 abs_peak=peak(bp, neff, pos), abs_ratio_065_025=ratio(bp, pos),
                                 beta_pow=r['beta_pow'], beta_pow_healthy=h['beta_pow'], sync=r['sync_corr10ms'],
                                 bursts=r['frac_spikes_in_bursts'],
                                 pair_PVP_e6=r['beta_pow_pair_PVP_pop'] * 1e6, pair_PVN_e6=r['beta_pow_pair_PVN_pop'] * 1e6,
                                 pair_mixed_e6=r['beta_pow_pair_cross'] * 1e6,
                                 rate_PVP=r['rate_PVP_pop'], rate_PVN=r['rate_PVN_pop'],
                                 rate_change_vs_cfgD_pct=100 * ((r['rate_PVP_pop'] * 12 + r['rate_PVN_pop'] * 46) / (ref[0] * 12 + ref[1] * 46) - 1)))
    return pd.DataFrame(rows), profs


def summary(df, profs):
    out = []
    for (variant, cond), g in df.groupby(['variant', 'condition'], sort=False):
        ks = [(variant, cond, s) for s in g.seed]
        pos = profs[ks[0]]['pos']
        neff = np.mean([profs[k]['neff'] for k in ks], axis=0)
        md = np.mean([profs[k]['diff'] for k in ks], axis=0)
        ma = np.mean([profs[k]['pd'] for k in ks], axis=0)
        pr = g.pair_PVP_e6 / g.pair_PVN_e6
        out.append(dict(variant=variant, condition=cond, depth_set=g.depth_set.iloc[0], pv_placement=g.pv_placement.iloc[0], n=len(g),
                        diff_peak_mean_profile=peak(md, neff, pos), diff_ratio_mean_profile=ratio(md, pos),
                        diff_peaks_per_seed=' '.join(f'{p:.2f}' for p in g.diff_peak),
                        seeds_diff_peak_055_080=int(((g.diff_peak >= 0.55) & (g.diff_peak <= 0.80)).sum()),
                        seeds_diff_peak_gt085=int((g.diff_peak > 0.85).sum()),
                        abs_peak_mean_profile=peak(ma, neff, pos), abs_ratio_mean_profile=ratio(ma, pos),
                        abs_peaks_per_seed=' '.join(f'{p:.2f}' for p in g.abs_peak),
                        beta_pow=g.beta_pow.mean(), beta_pow_sd=g.beta_pow.std(ddof=1), sync=g.sync.mean(), bursts=g.bursts.mean(),
                        pair_PVP_e6=g.pair_PVP_e6.mean(), pair_PVN_e6=g.pair_PVN_e6.mean(), pair_mixed_e6=g.pair_mixed_e6.mean(),
                        pvp_pvn_ratio_of_means=g.pair_PVP_e6.mean() / g.pair_PVN_e6.mean(),
                        pvp_pvn_ratio_seed_mean=pr.mean(), pvp_pvn_ratio_seed_sd=pr.std(ddof=1),
                        seeds_pvp_gt_pvn=int((g.pair_PVP_e6 > g.pair_PVN_e6).sum()),
                        rate_PVP=g.rate_PVP.mean(), rate_PVN=g.rate_PVN.mean(), rate_change_vs_cfgD_pct=g.rate_change_vs_cfgD_pct.mean()))
    return pd.DataFrame(out)


def verdict(sm):
    res = {}
    for variant in ('c', 'a', 'b'):
        S = {r.condition: r for r in sm[sm.variant == variant].itertuples()}
        ctl = {r.condition: r for r in sm[sm.condition.isin(CONTROLS)].itertuples()}
        S.update(ctl)
        if not all(k in S for k in ('1_uniform_control', '2_territory_adv', '3_pv_uniform', '6_pv_monotonic')):
            continue
        c1, c2, c3, c6 = (S[k] for k in ('1_uniform_control', '2_territory_adv', '3_pv_uniform', '6_pv_monotonic'))
        v = {'1_cond2_peak_in_055_080': bool(0.55 <= c2.diff_peak_mean_profile <= 0.80),
             '2_cond3_peak_gt_085': bool(c3.diff_peak_mean_profile > 0.85),
             '3_cond6_peak_gt_085': bool(c6.diff_peak_mean_profile > 0.85),
             '4_cond2_ratio_gt_cond1': bool(np.isfinite(c2.diff_ratio_mean_profile) and np.isfinite(c1.diff_ratio_mean_profile)
                                           and c2.diff_ratio_mean_profile > c1.diff_ratio_mean_profile)}
        v['hypothesis_survives'] = all(v.values())
        v['values'] = dict(cond1_peak=c1.diff_peak_mean_profile, cond1_ratio=c1.diff_ratio_mean_profile,
                           cond2_peak=c2.diff_peak_mean_profile, cond2_ratio=c2.diff_ratio_mean_profile,
                           cond3_peak=c3.diff_peak_mean_profile, cond6_peak=c6.diff_peak_mean_profile)
        if '2b_territory_adv_pvalb' in S:
            b = S['2b_territory_adv_pvalb']
            v['2b_rule1_peak_in_055_080'] = bool(0.55 <= b.diff_peak_mean_profile <= 0.80)
            v['2b_rule4_ratio_gt_cond1'] = bool(np.isfinite(b.diff_ratio_mean_profile) and b.diff_ratio_mean_profile > c1.diff_ratio_mean_profile)
            v['values'].update(cond2b_peak=b.diff_peak_mean_profile, cond2b_ratio=b.diff_ratio_mean_profile)
        if '4_mouse_map' in S:
            v['values']['cond4_peak'] = S['4_mouse_map'].diff_peak_mean_profile
        if '5_territory_early' in S and '5c_uniform_control_early' in S:
            c5, c5c = S['5_territory_early'], S['5c_uniform_control_early']
            v['values'].update(cond5_peak=c5.diff_peak_mean_profile, cond5_gain_over_own_control=c5.diff_ratio_mean_profile / c5c.diff_ratio_mean_profile,
                               cond2_gain_over_control=c2.diff_ratio_mean_profile / c1.diff_ratio_mean_profile)
        res[variant] = v
    return res


def figure(sm, profs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    variants = [v for v in ('c', 'a', 'b') if (sm.variant == v).any()]
    fig, axs = plt.subplots(2, len(variants), figsize=(4.6 * len(variants), 7), sharex=True, squeeze=False, constrained_layout=True)
    cols = dict(zip(ADV + EARLY, ['0.2', '#d0632b', '#8c2d04', '#2b7bba', '#5aa02c', '#9467bd', '0.2', '#e7a33e']))
    for j, v in enumerate(variants):
        for i, (group, title) in enumerate(((ADV, f'variant {v}: advanced set (M 0.8825)'), (EARLY, f'variant {v}: early set (M 0.6757)'))):
            a = axs[i, j]
            for cond in group:
                ks = [k for k in profs if k[1] == cond and (k[0] == v or cond in CONTROLS) and (cond not in CONTROLS or k[0] == 'c')]
                if not ks:
                    continue
                x = np.array(profs[ks[0]]['pos'])
                md = np.mean([profs[k]['diff'] for k in ks], axis=0) * 1e6
                a.plot(x, md, color=cols[cond], ls='--' if cond in CONTROLS else '-', lw=1.6, label=f'{cond} (n={len(ks)})')
            a.axvspan(0.55, 0.80, color='0.92', zorder=0)
            a.axvline(0.85, color='0.6', lw=0.8, ls=':')
            a.set_title(title, fontsize=9)
            a.legend(fontsize=6)
            a.set_ylabel('PD - healthy shared beta (x1e-6)')
        axs[1, j].set_xlabel('STN axis (0 ventromedial, 1 dorsolateral)')
    fig.savefig(os.path.join(HERE, 'difference_maps.png'), dpi=130)


def main():
    df, profs = runs_table()
    df.to_csv(os.path.join(HERE, 'drive_runs.csv'), index=False, float_format='%.5g')
    sm = summary(df, profs)
    sm.to_csv(os.path.join(HERE, 'drive_summary.csv'), index=False, float_format='%.4g')
    v = verdict(sm)
    json.dump(v, open(os.path.join(HERE, 'verdict.json'), 'w'), indent=1, default=float)
    figure(sm, profs)
    pd.set_option('display.width', 250)
    print(sm[['variant', 'condition', 'n', 'diff_peak_mean_profile', 'diff_peaks_per_seed', 'diff_ratio_mean_profile', 'abs_peak_mean_profile',
              'abs_ratio_mean_profile', 'beta_pow', 'sync', 'bursts', 'pvp_pvn_ratio_of_means', 'seeds_pvp_gt_pvn', 'rate_change_vs_cfgD_pct']].to_string(index=False))
    print(json.dumps(v, indent=1, default=float))


if __name__ == '__main__':
    main()
