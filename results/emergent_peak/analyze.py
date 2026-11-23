"""
Step 2 analysis for the emergent beta peak test (criteria fixed in the task before any run; see summary.md).

    python results/emergent_peak/analyze.py

Reads sweeps/out/emergent_peak/<label>_s<seed>/model_output_metrics.json. Per condition and state: seed-mean shared-beta
profile along the axis (network_metrics: weighted mean pairwise beta cross-power, positions with n_eff >= 5 count for
the peak), per-seed peak positions (beta_peak_axis), 0.65 / 0.25 ratio of the seed-mean profile (and per seed), PD -
healthy contrast profile, synchrony profile, PV+ / PV- shared beta, rates (vs configuration D, results/gaba_ratio:
healthy 37.5 / 20.3 Hz, PD 32.9 / 16.4 Hz; |change| > 15% is flagged, not corrected).

Criteria:
- Emergent peak (F): mean-profile maximum in [0.55, 0.75] AND per-seed maximum in [0.50, 0.80] in >= 4/5 seeds AND
  0.65/0.25 ratio >= 1.5.
- Each ingredient necessary: in noPV and in noRec the criterion fails, defined as ratio < 1.2 OR mean-profile maximum
  outside [0.50, 0.80].
- abs vs F: per-type ratio beyond the absolute gradients; sig vs ref: data-derived GPe topography; F vs F_cfgD:
  measured vs exaggerated PV+ placement.
- M (three loci): reported separately, "assumes motor-territory origin of beta"; not part of the verdict.
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(ROOT, 'sweeps', 'out', 'emergent_peak')
SEEDS = [1, 2, 3, 4, 5]
CFGD = {'healthy': (37.4867, 20.2591), 'pd': (32.8933, 16.4365)}
COND = ['F', 'F_cfgD', 'noPV', 'noRec', 'abs', 'sig', 'ref', 'M0', 'M1', 'M2']
MIN_NEFF = 5


def load(label, out=None):
    runs = []
    for s in SEEDS:
        f = os.path.join(out or OUT, f'{label}_s{s}', 'model_output_metrics.json')
        if os.path.exists(f):
            runs.append(json.load(open(f)))
    return runs


def ratio(prof, pos):
    i65, i25 = pos.index(0.65), pos.index(0.25)
    den = prof[i25]
    return float(prof[i65] / den) if den > 0 else np.nan


def peak(prof, neff, pos):
    ok = np.array(neff) >= MIN_NEFF
    p = np.where(ok, prof, -np.inf)
    return float(pos[int(np.argmax(p))])


def summarize(label, out=None):
    runs = load(label, out)
    if not runs:
        return None
    pos = [round(x, 2) for x in runs[0]['profile']['positions']]
    B = np.array([r['profile']['beta_pow'] for r in runs])
    S = np.array([r['profile']['sync'] for r in runs])
    N = np.array([r['profile']['n_eff'] for r in runs]).mean(0)
    mean_b = B.mean(0)
    seed_peaks = [float(r['beta_peak_axis']) for r in runs]
    seed_ratios = [ratio(b, pos) for b in B]
    return dict(label=label, n=len(runs), positions=pos, beta_profile=mean_b, beta_sd=B.std(0, ddof=1) if len(runs) > 1 else np.zeros_like(mean_b),
                sync_profile=S.mean(0), mean_peak=peak(mean_b, N, pos), seed_peaks=seed_peaks,
                ratio_mean_profile=ratio(mean_b, pos), seed_ratios=seed_ratios,
                beta_pow=np.mean([r['beta_pow'] for r in runs]), beta_frac=np.mean([r['beta_frac'] for r in runs]),
                sync=np.mean([r['sync_corr10ms'] for r in runs]),
                pair_PVP=np.mean([r['beta_pow_pair_PVP_pop'] for r in runs]) * 1e6, pair_PVN=np.mean([r['beta_pow_pair_PVN_pop'] for r in runs]) * 1e6,
                pair_cross=np.mean([r['beta_pow_pair_cross'] for r in runs]) * 1e6,
                rate_PVP=np.mean([r['rate_PVP_pop'] for r in runs]), rate_PVN=np.mean([r['rate_PVN_pop'] for r in runs]),
                bursts=np.mean([r['frac_spikes_in_bursts'] for r in runs]),
                pair_ratio_seeds=[r['beta_pow_pair_PVP_pop'] / r['beta_pow_pair_PVN_pop'] for r in runs])


def verdict(s):
    if s is None:
        return 'missing'
    in_mean = 0.55 <= s['mean_peak'] <= 0.75
    n_seed = sum(0.50 <= p <= 0.80 for p in s['seed_peaks'])
    r = s['ratio_mean_profile']
    ok = in_mean and n_seed >= 4 and np.isfinite(r) and r >= 1.5
    return dict(emergent=bool(ok), mean_peak_in_055_075=bool(in_mean), seeds_peak_in_050_080=int(n_seed), ratio=r,
                fails_necessity_rule=bool((not np.isfinite(r)) or r < 1.2 or not (0.50 <= s['mean_peak'] <= 0.80)))


def main():
    rows, profs = [], {}
    for c in COND:
        for st in ('healthy', 'pd'):
            s = summarize(f'{c}_{st}')
            if s is None:
                continue
            profs[(c, st)] = s
            v = verdict(s) if st == 'pd' else {}
            ref = CFGD[st]
            rows.append(dict(condition=c, state=st, n=s['n'], mean_peak=s['mean_peak'],
                             seed_peaks=' '.join(f'{p:.2f}' for p in s['seed_peaks']),
                             ratio_065_025=s['ratio_mean_profile'],
                             seed_ratio_mean=np.nanmean(s['seed_ratios']), seed_ratio_sd=np.nanstd(s['seed_ratios'], ddof=1),
                             beta_pow=s['beta_pow'], beta_frac=s['beta_frac'], sync=s['sync'],
                             pair_PVP_e6=s['pair_PVP'], pair_PVN_e6=s['pair_PVN'], pair_cross_e6=s['pair_cross'],
                             rate_PVP=s['rate_PVP'], rate_PVN=s['rate_PVN'],
                             rate_change_vs_cfgD_pct=100 * ((s['rate_PVP'] * 12 + s['rate_PVN'] * 46) / (ref[0] * 12 + ref[1] * 46) - 1),
                             bursts=s['bursts'], **v))
    df = pd.DataFrame(rows)
    df['rate_flag_gt15pct'] = df.rate_change_vs_cfgD_pct.abs() > 15
    df.to_csv(os.path.join(HERE, 'conditions.csv'), index=False, float_format='%.4g')
    prof_rows = []
    for (c, st), s in profs.items():
        for i, p in enumerate(s['positions']):
            h = profs.get(('F' if c.startswith('M') else c, 'healthy'))
            prof_rows.append(dict(condition=c, state=st, axis=p, beta=s['beta_profile'][i], beta_sd=s['beta_sd'][i], sync=s['sync_profile'][i],
                                  pd_minus_healthy=(s['beta_profile'][i] - h['beta_profile'][i]) if (st == 'pd' and h) else np.nan))
    pd.DataFrame(prof_rows).to_csv(os.path.join(HERE, 'profiles.csv'), index=False, float_format='%.5g')
    figure(profs)
    pd.set_option('display.width', 250)
    print(df[['condition', 'state', 'n', 'mean_peak', 'seed_peaks', 'ratio_065_025', 'beta_frac', 'pair_PVP_e6', 'pair_PVN_e6',
              'rate_PVP', 'rate_PVN', 'rate_change_vs_cfgD_pct', 'emergent', 'seeds_peak_in_050_080', 'fails_necessity_rule']].to_string(index=False))


def figure(profs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    conds = [c for c in COND if (c, 'pd') in profs]
    fig, axs = plt.subplots(2, max(len(conds), 1), figsize=(2.6 * len(conds), 5.6), sharex=True, constrained_layout=True)
    axs = np.atleast_2d(axs)
    for j, c in enumerate(conds):
        s = profs[(c, 'pd')]
        x = np.array(s['positions'])
        a = axs[0, j]
        a.fill_between(x, s['beta_profile'] - s['beta_sd'], s['beta_profile'] + s['beta_sd'], color='#d0632b', alpha=0.2, lw=0)
        a.plot(x, s['beta_profile'], color='#d0632b', label='PD')
        h = profs.get(('F' if c.startswith('M') else c, 'healthy'))
        if h:
            a.plot(x, h['beta_profile'], color='0.4', ls='--', label='healthy')
        a.axvspan(0.55, 0.75, color='0.9', zorder=0)
        a.set_title(f"{c}: peak {s['mean_peak']:.2f}, 0.65/0.25 {s['ratio_mean_profile']:.2f}", fontsize=8)
        axs[1, j].plot(x, s['sync_profile'], color='#d0632b')
        if h:
            axs[1, j].plot(x, h['sync_profile'], color='0.4', ls='--')
        axs[1, j].set_xlabel('STN axis')
    axs[0, 0].set_ylabel('shared beta (pairwise cross-power)')
    axs[1, 0].set_ylabel('synchrony')
    axs[0, 0].legend(fontsize=7)
    fig.savefig(os.path.join(HERE, 'profiles.png'), dpi=130)


if __name__ == '__main__':
    main()
