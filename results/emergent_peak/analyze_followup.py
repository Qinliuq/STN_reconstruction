"""
POST-HOC follow-ups FU1 / FU2 and configuration-D reference reruns (sweeps/emergent_peak_followup.py); same criteria
and readouts as analyze.py (imported, unchanged), plus the PV+ / PV- pairwise shared-beta ratio across ALL conditions
of the emergent-peak test (Step 2 + follow-ups).

    python results/emergent_peak/analyze_followup.py

PV+ / PV- readout: per-seed ratio beta_pow_pair_PVP_pop / beta_pow_pair_PVN_pop (mean +- SD over seeds, as asked);
because pairwise cross-power can be ~0 or negative (healthy runs), also the ratio of seed means, the per-seed
difference PV+ - PV- (mean +- SD, x1e-6) and the number of seeds with PV+ > PV-.

Writes followup_conditions.csv, followup_profiles.csv, followup_profiles.png, pv_pair_ratio_all.csv.
"""
import os

import numpy as np
import pandas as pd

import analyze as A

FU_OUT = os.path.join(A.ROOT, 'sweeps', 'out', 'emergent_peak_followup')
# condition -> (outdir, states, healthy reference label)
FU = {'FU1': (FU_OUT, ('healthy', 'pd'), 'FU1'), 'FU2_F': (FU_OUT, ('healthy', 'pd'), 'FU2_F'),
      'FU2_M0': (FU_OUT, ('pd',), 'FU2_F'), 'FU2_M1': (FU_OUT, ('pd',), 'FU2_F'), 'FU2_M2': (FU_OUT, ('pd',), 'FU2_F'),
      'Dref': (FU_OUT, ('healthy', 'pd'), 'Dref'), 'Dref_uni': (FU_OUT, ('pd',), 'Dref')}
STEP2 = {c: (A.OUT, ('healthy', 'pd') if not c.startswith('M') else ('pd',), 'F' if c.startswith('M') else c)
         for c in A.COND}


def table(spec):
    rows, profs = [], {}
    for c, (out, states, href) in spec.items():
        for st in states:
            s = A.summarize(f'{c}_{st}', out)
            if s is None:
                continue
            profs[(c, st)] = (s, out, href)
            v = A.verdict(s) if st == 'pd' else {}
            ref = A.CFGD[st]
            pr = np.array(s['pair_ratio_seeds'])
            runs = A.load(f'{c}_{st}', out)
            pp = np.array([r['beta_pow_pair_PVP_pop'] for r in runs]) * 1e6
            pn = np.array([r['beta_pow_pair_PVN_pop'] for r in runs]) * 1e6
            rows.append(dict(condition=c, state=st, n=s['n'], mean_peak=s['mean_peak'],
                             seed_peaks=' '.join(f'{p:.2f}' for p in s['seed_peaks']), ratio_065_025=s['ratio_mean_profile'],
                             beta_pow=s['beta_pow'], beta_frac=s['beta_frac'], sync=s['sync'],
                             pair_PVP_e6=s['pair_PVP'], pair_PVN_e6=s['pair_PVN'],
                             pvp_pvn_pair_ratio_mean=pr.mean(), pvp_pvn_pair_ratio_sd=pr.std(ddof=1),
                             ratio_of_seed_means=pp.mean() / pn.mean(), diff_e6_mean=(pp - pn).mean(), diff_e6_sd=(pp - pn).std(ddof=1),
                             seeds_pvp_gt_pvn=int((pp > pn).sum()), min_pvn_e6=pn.min(),
                             rate_PVP=s['rate_PVP'], rate_PVN=s['rate_PVN'],
                             rate_change_vs_cfgD_pct=100 * ((s['rate_PVP'] * 12 + s['rate_PVN'] * 46) / (ref[0] * 12 + ref[1] * 46) - 1),
                             bursts=s['bursts'], **v))
    df = pd.DataFrame(rows)
    df['rate_flag_gt15pct'] = df.rate_change_vs_cfgD_pct.abs() > 15
    return df, profs


def main():
    fu, fprofs = table(FU)
    fu.to_csv(os.path.join(A.HERE, 'followup_conditions.csv'), index=False, float_format='%.4g')
    rows = []
    for (c, st), (s, out, href) in fprofs.items():
        h = fprofs.get((href, 'healthy'))
        for i, p in enumerate(s['positions']):
            rows.append(dict(condition=c, state=st, axis=p, beta=s['beta_profile'][i], beta_sd=s['beta_sd'][i], sync=s['sync_profile'][i],
                             pd_minus_healthy=(s['beta_profile'][i] - h[0]['beta_profile'][i]) if (st == 'pd' and h) else np.nan))
    pd.DataFrame(rows).to_csv(os.path.join(A.HERE, 'followup_profiles.csv'), index=False, float_format='%.5g')
    figure(fprofs)
    st2, _ = table(STEP2)
    allc = pd.concat([st2.assign(set='Step 2'), fu.assign(set='follow-up (post hoc)')])
    allc[['set', 'condition', 'state', 'n', 'pvp_pvn_pair_ratio_mean', 'pvp_pvn_pair_ratio_sd', 'ratio_of_seed_means',
          'diff_e6_mean', 'diff_e6_sd', 'seeds_pvp_gt_pvn', 'min_pvn_e6', 'pair_PVP_e6', 'pair_PVN_e6',
          'ratio_065_025', 'rate_change_vs_cfgD_pct']].to_csv(os.path.join(A.HERE, 'pv_pair_ratio_all.csv'), index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    print(fu[['condition', 'state', 'n', 'mean_peak', 'seed_peaks', 'ratio_065_025', 'beta_frac', 'pair_PVP_e6', 'pair_PVN_e6',
              'rate_PVP', 'rate_PVN', 'rate_change_vs_cfgD_pct', 'bursts', 'emergent', 'seeds_peak_in_050_080', 'fails_necessity_rule']].to_string(index=False))
    print(allc[['set', 'condition', 'state', 'pvp_pvn_pair_ratio_mean', 'pvp_pvn_pair_ratio_sd', 'ratio_of_seed_means',
                'diff_e6_mean', 'diff_e6_sd', 'seeds_pvp_gt_pvn', 'min_pvn_e6']].to_string(index=False))


def figure(profs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    conds = [c for c in FU if (c, 'pd') in profs]
    fig, axs = plt.subplots(2, len(conds), figsize=(2.6 * len(conds), 5.6), sharex=True, constrained_layout=True)
    axs = np.atleast_2d(axs)
    for j, c in enumerate(conds):
        s, _, href = profs[(c, 'pd')]
        h = profs.get((href, 'healthy'))
        x = np.array(s['positions'])
        a = axs[0, j]
        a.fill_between(x, s['beta_profile'] - s['beta_sd'], s['beta_profile'] + s['beta_sd'], color='#d0632b', alpha=0.2, lw=0)
        a.plot(x, s['beta_profile'], color='#d0632b', label='PD')
        if h:
            a.plot(x, h[0]['beta_profile'], color='0.4', ls='--', label='healthy')
            axs[1, j].plot(x, h[0]['sync_profile'], color='0.4', ls='--')
        a.axvspan(0.55, 0.75, color='0.9', zorder=0)
        a.set_title(f"{c} (post hoc): peak {s['mean_peak']:.2f}, 0.65/0.25 {s['ratio_mean_profile']:.2f}", fontsize=7)
        axs[1, j].plot(x, s['sync_profile'], color='#d0632b')
        axs[1, j].set_xlabel('STN axis')
    axs[0, 0].set_ylabel('shared beta (pairwise cross-power)')
    axs[1, 0].set_ylabel('synchrony')
    axs[0, 0].legend(fontsize=7)
    fig.savefig(os.path.join(A.HERE, 'followup_profiles.png'), dpi=130)


if __name__ == '__main__':
    main()
