"""
H5: healthy -> parkinsonian change, configuration D vs awake primate STN (Galvan & Wichmann; Zenodo 6581800).

Refuses to run without results/heldout/H5/prediction.json (written by model/predict_H5.py before the data were
downloaded). Rules fixed here before the first run (2026-11-18); the user's adjustments of 2026-11-18 are included.

Data units (13 normal, 27 MPTP; awake, at rest) are length-matched to the model's 20 s trains: for each unit, 10 random
20 s segments (uniform start, fixed seed); every feature is computed per segment with the same function as the model
(metrics.unit_summary.h5_summary) and the per-unit value is the median over segments; a unit is oscillatory if the
majority of its segments are, and its peak is the median peak over its oscillatory segments. Units shorter than 20 s
use their full train (counted). Full-length values are reported as secondary.

Primary quantities are CHANGES healthy -> PD of population medians. Data 95% CI: 10,000 bootstrap resamples of units
within each state (independently). Model: mean +- SD over the 5 seeds of the per-seed change (the PD and healthy runs
of a seed share the network and the seed).

1. Rate change (ratio PD / healthy of the median rate): match if the model seed-mean ratio lies inside the data CI;
   partial if > 1 but outside it; mismatch if <= 1.
2. Irregularity unchanged (PD - healthy differences of the median CV, LV and pause-time fraction): match if all three
   model seed-mean differences lie inside the data CIs; partial if 1-2; mismatch if none.
   The calibration of the PD state targeted PD >= 1.5x burstier than healthy, so this sub-test also tests that
   calibration premise.
3. Oscillatory fraction, magnitude only. The direction (PD > healthy) is built into the PD definition and does not
   count. Match if the model PD fraction (seed mean) lies inside the Wilson 95% CI of the data's MPTP fraction AND the
   model healthy fraction inside the CI of the data's normal fraction; partial if one of the two; mismatch if neither.
   Primary data fractions: the length-matched counts (the comparison the length matching exists for); secondary: the
   full-length counts (reference 6/13 and 24/27).
4. Oscillation frequency ("tests the input assumption": the model's PD frequency is set by the imposed 20 Hz pallidal
   drive, not by the STN): match if >= 50% of the model's PD oscillatory cells (seed mean of the per-seed share) peak at
   4-17 Hz; else mismatch.
No overall verdict beyond the four lines. Also reported: the same changes for PV+ and PV- model cells (data have no
PV labels); within-state distributions (Mann-Whitney + KS; model cells pooled over seeds vs data units) as secondary.
Within-state levels are labelled "unclear" (Step 0): healthy rate levels were calibrated to "plausible ~30/20 Hz" with
no stated source, so the healthy PV- rate of the model (~20 Hz) next to the primate normal median (~23 Hz) is not
evidence of validity.

    python -m validation.heldout.compare.compare_H5
"""
import os

import numpy as np
import pandas as pd

from validation.heldout import data_io
from validation.heldout.compare import common as C
from validation.heldout.metrics import stats
from validation.heldout.metrics.unit_summary import h5_summary, FEATURES

SEG_S = 20.0
N_SEG = 10
N_BOOT = 10000
STEP0 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'H5_step0_table.md'), encoding='utf-8').read()


def unit_values(units):
    """Per-unit length-matched (median over 20 s segments) and full-length summaries."""
    lm, full = [], []
    for k, u in enumerate(units.itertuples()):
        st = np.concatenate([[0.0], np.cumsum(u.isi_ms)]) / 1000.0
        T = st[-1]
        rng = np.random.default_rng(10_000 + k)
        if T > SEG_S:
            starts = rng.uniform(0, T - SEG_S, N_SEG)
            segs = [h5_summary(st, a, a + SEG_S, seed=100 * k + j) for j, a in enumerate(starts)]
        else:
            segs = [h5_summary(st, 0.0, T + 1e-9, seed=100 * k)]
        segs = [s for s in segs if s is not None]
        base = dict(unit=u.unit, state=u.state, duration_s=T, n_segments=len(segs), short=T <= SEG_S)
        if segs:
            S = pd.DataFrame(segs)
            row = {**base, **{f: S[f].median() for f in FEATURES}}
            row['oscillatory'] = bool(S.oscillatory.mean() > 0.5)
            row['osc_segment_frac'] = S.oscillatory.mean()
            row['peak_hz'] = S.loc[S.oscillatory, 'peak_hz'].median() if S.oscillatory.any() else np.nan
            for c in [c for c in S.columns if c.startswith('ratio_')]:
                row[c] = S[c].median()
            lm.append(row)
        f = h5_summary(st, 0.0, T + 1e-9, seed=100 * k + 99)
        if f is not None:
            full.append({**base, **f})
    return pd.DataFrame(lm), pd.DataFrame(full)


def boot_change(h, p, kind, rng):
    h, p = np.asarray(h, float), np.asarray(p, float)
    bh = np.median(rng.choice(h, (N_BOOT, len(h))), axis=1)
    bp = np.median(rng.choice(p, (N_BOOT, len(p))), axis=1)
    b = bp / bh if kind == 'ratio' else bp - bh
    est = np.median(p) / np.median(h) if kind == 'ratio' else np.median(p) - np.median(h)
    return float(est), (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))


def inside(x, ci):
    return ci[0] <= x <= ci[1]


def main():
    pred = C.require_prediction('H5')
    M = pred['metrics']
    units = data_io.primate_stn_units()
    lm, full = unit_values(units)
    out = C.out_dir('H5')
    lm.to_csv(os.path.join(out, 'data_units_length_matched.csv'), index=False, float_format='%.6g')
    full.to_csv(os.path.join(out, 'data_units_full_length.csv'), index=False, float_format='%.6g')
    rng = np.random.default_rng(0)
    H, P = lm[lm.state == 'normal'], lm[lm.state == 'mptp']

    # ---- 1. rate change
    rate_est, rate_ci = boot_change(H.rate_hz, P.rate_hz, 'ratio', rng)
    m_rate = M['all']['rate_ratio']['mean']
    v1 = 'match' if inside(m_rate, rate_ci) else ('partial' if m_rate > 1 else 'mismatch')
    # ---- 2. irregularity
    irr = {}
    for k in ('cv', 'lv', 'pause_time_frac'):
        est, ci = boot_change(H[k], P[k], 'diff', rng)
        irr[k] = dict(data=est, ci=ci, model=M['all'][f'{k}_diff']['mean'], model_sd=M['all'][f'{k}_diff']['sd'])
        irr[k]['inside'] = inside(irr[k]['model'], ci)
    n_in = sum(v['inside'] for v in irr.values())
    v2 = 'match' if n_in == 3 else ('partial' if n_in >= 1 else 'mismatch')
    # ---- 3. oscillatory fraction (magnitude)
    kh, nh, kp, np_ = int(H.oscillatory.sum()), len(H), int(P.oscillatory.sum()), len(P)
    ci_h, ci_p = C.wilson(kh, nh), C.wilson(kp, np_)
    mh, mp = M['all']['osc_frac_healthy']['mean'], M['all']['osc_frac_pd']['mean']
    n3 = inside(mp, ci_p) + inside(mh, ci_h)
    v3 = 'match' if n3 == 2 else ('partial' if n3 == 1 else 'mismatch')
    fh_ = full[full.state == 'normal']
    fp_ = full[full.state == 'mptp']
    kh_f, kp_f = int(fh_.oscillatory.sum()), int(fp_.oscillatory.sum())
    n3f = inside(mp, C.wilson(kp_f, len(fp_))) + inside(mh, C.wilson(kh_f, len(fh_)))
    v3f = 'match' if n3f == 2 else ('partial' if n3f == 1 else 'mismatch')
    # ---- 4. oscillation frequency
    m_peak = M['all']['pd_peak_4_17_frac']['mean']
    v4 = 'match' if m_peak >= 0.5 else 'mismatch'
    dpk = P.loc[P.oscillatory, 'peak_hz']
    d_peak = float(((dpk >= 4) & (dpk < 17)).mean()) if len(dpk) else np.nan

    lines = {
        '1': f"Rate change: {v1.upper()} — model PD/healthy {C.fmt(M['all']['rate_ratio'])} vs data {rate_est:.2f} "
             f"(95% CI {rate_ci[0]:.2f}-{rate_ci[1]:.2f})",
        '2': f"Irregularity unchanged: {v2.upper()} — {n_in}/3 model differences inside the data CIs (" +
             '; '.join(f"{k} {v['model']:+.3f} vs {v['data']:+.3f} [{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}]" for k, v in irr.items()) + ')',
        '3': f"Oscillatory fraction (magnitude only): {v3.upper()} — model healthy {C.fmt(M['all']['osc_frac_healthy'])}, "
             f"PD {C.fmt(M['all']['osc_frac_pd'])} vs data (length-matched) normal {kh}/{nh} [{ci_h[0]:.2f}-{ci_h[1]:.2f}], "
             f"MPTP {kp}/{np_} [{ci_p[0]:.2f}-{ci_p[1]:.2f}]; full-length {kh_f}/{len(fh_)} and {kp_f}/{len(fp_)}: {v3f}",
        '4': f"Oscillation frequency (tests the input assumption): {v4.upper()} — model PD oscillatory cells peaking at "
             f"4-17 Hz {C.fmt(M['all']['pd_peak_4_17_frac'])} (median peak {C.fmt(M['all']['pd_peak_median_hz'], 1)} Hz); "
             f"data MPTP {d_peak:.2f} (median {dpk.median():.1f} Hz)",
    }
    verdict = (f"H5: rate change {v1}; irregularity {v2}; oscillatory fraction {v3}; oscillation frequency {v4} "
               f"(input assumption)")

    # ---- PV split (model; against the same data CIs)
    pv_rows = []
    for g in ('all', 'PV+', 'PV-'):
        r = dict(group=g, rate_ratio=C.fmt(M[g]['rate_ratio']))
        for k in ('cv', 'lv', 'pause_time_frac'):
            r[f'{k}_diff'] = C.fmt(M[g][f'{k}_diff'], 3)
        r.update(osc_healthy=C.fmt(M[g]['osc_frac_healthy']), osc_pd=C.fmt(M[g]['osc_frac_pd']),
                 pd_peak_4_17=C.fmt(M[g]['pd_peak_4_17_frac']))
        pv_rows.append(r)
    data_row = dict(group='data (length-matched)', rate_ratio=f'{rate_est:.2f} [{rate_ci[0]:.2f}, {rate_ci[1]:.2f}]')
    for k, v in irr.items():
        data_row[f'{k}_diff'] = f"{v['data']:+.3f} [{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}]"
    data_row.update(osc_healthy=f'{kh}/{nh}', osc_pd=f'{kp}/{np_}', pd_peak_4_17=f'{d_peak:.2f}')
    pv = pd.DataFrame(pv_rows + [data_row])
    pv.to_csv(os.path.join(out, 'changes.csv'), index=False)

    # ---- secondary: within-state distributions (levels: "unclear", Step 0)
    cells = pd.read_csv(os.path.join(C.ROOT, 'results', 'heldout', 'H5', 'model_cells.csv'))
    sec = []
    for state_m, state_d, D in (('healthy', 'normal', H), ('parkinsonian', 'mptp', P)):
        for k in FEATURES:
            r = stats.compare_distributions(cells.loc[cells.state == state_m, k].dropna(), D[k].dropna())
            sec.append(dict(state=state_m, feature=k, model_median=cells.loc[cells.state == state_m, k].median(),
                            data_median=D[k].median(), **{kk: r[kk] for kk in ('mw_p', 'ks_p')}))
    sec = pd.DataFrame(sec)
    sec.to_csv(os.path.join(out, 'within_state_secondary.csv'), index=False, float_format='%.4g')

    figure(cells, lm, M, rate_est, rate_ci, irr, (kh, nh, kp, np_), out)

    sections = [
        ('Sub-tests (one line each; no overall verdict)', '\n'.join(f'{i}. {l}' for i, l in lines.items()) +
         '\n\nThe calibration of the PD state targeted PD >= 1.5x burstier than healthy (results/topo_pd_state), so sub-test 2 '
         'also tests that calibration premise.\n\nSub-test 3 judges magnitude only: the direction (PD more oscillatory) is '
         'built into the PD definition. Sub-test 4 tests the input assumption (imposed 20 Hz pallidal drive), not the STN.'),
        ('Changes by PV type (model) and data (length-matched; 95% bootstrap CI)', C.md_table(pv)),
        ('Within-state levels (secondary; labelled "unclear")',
         'Within-state levels are **unclear** for validation (Step 0): healthy rates were calibrated to "plausible ~30/20 Hz" '
         'with no documented source. The model healthy PV- rate (~20 Hz) next to the primate normal median is **not** '
         'presented as validation.\n\n' + C.md_table(sec)),
        ('Data', f'{nh} normal / {np_} MPTP STN units (length-matched; {int(lm.short.sum())} shorter than 20 s use the full train). '
                 f'Full-length medians: normal rate {fh_.rate_hz.median():.1f} Hz, MPTP {fp_.rate_hz.median():.1f} Hz.'),
        ('Step 0: calibration-leak check (verbatim)', STEP0),
        ('Caveats', 'Awake primate vs a mouse/rat-parameterized model; no anatomical coordinates per unit (the axis gradient '
                    'cannot be tested); the model PD oscillation frequency is imposed by the pallidal drive; 2 animals, '
                    'units are pseudo-replicates; model seeds share one 58-cell network.'),
    ]
    C.write_summary('H5', verdict, sections, pred, table=dict(
        verdict_word=f'rate {v1} / irregularity {v2} / osc. fraction {v3} / frequency {v4}',
        model=f"rate x{M['all']['rate_ratio']['mean']:.2f}; osc {mh:.2f} -> {mp:.2f}",
        data=f"rate x{rate_est:.2f} [{rate_ci[0]:.2f}, {rate_ci[1]:.2f}]; osc {kh}/{nh} -> {kp}/{np_}",
        statistic='change of population medians, unit bootstrap CI; Wilson CI',
        subtests=lines))


def figure(cells, lm, M, rate_est, rate_ci, irr, counts, out):
    plt = C.style()
    fig, axs = plt.subplots(2, 4, figsize=(13, 6), constrained_layout=True)
    for ax, k in zip(axs[0], ['rate_hz', 'cv', 'lv', 'pause_time_frac']):
        for sm, sd, ls in (('healthy', 'normal', '--'), ('parkinsonian', 'mptp', '-')):
            for src, x, col in (('model', cells.loc[cells.state == sm, k], C.COL['model']),
                                ('data', lm.loc[lm.state == sd, k], C.COL['data'])):
                x = np.sort(x.dropna())
                ax.step(x, np.arange(1, len(x) + 1) / len(x), ls=ls, color=col, label=f'{src} {sm[:7]}')
        ax.set_title(k)
    axs[0, 0].legend(fontsize=6)
    # change plot
    ax = axs[1, 0]
    names = ['rate ratio', 'CV diff', 'LV diff', 'pause diff']
    dvals = [(rate_est, rate_ci)] + [(irr[k]['data'], irr[k]['ci']) for k in ('cv', 'lv', 'pause_time_frac')]
    mvals = [M['all']['rate_ratio']] + [M['all'][f'{k}_diff'] for k in ('cv', 'lv', 'pause_time_frac')]
    for i, ((d, ci), m) in enumerate(zip(dvals, mvals)):
        ax.errorbar(i - 0.1, d, yerr=[[d - ci[0]], [ci[1] - d]], fmt='o', color=C.COL['data'])
        ax.errorbar(i + 0.1, m['mean'], yerr=m['sd'], fmt='s', color=C.COL['model'])
    ax.axhline(0, color='k', lw=0.5)
    ax.axhline(1, color='k', lw=0.5, ls=':')
    ax.set_xticks(range(4), names)
    ax.set_title('change healthy -> PD (data CI, model seed SD)')
    # oscillatory fractions
    ax = axs[1, 1]
    kh, nh, kp, np_ = counts
    for i, (k, n, mm) in enumerate(((kh, nh, M['all']['osc_frac_healthy']), (kp, np_, M['all']['osc_frac_pd']))):
        lo, hi = C.wilson(k, n)
        ax.errorbar(i - 0.1, k / n, yerr=[[k / n - lo], [hi - k / n]], fmt='o', color=C.COL['data'])
        ax.errorbar(i + 0.1, mm['mean'], yerr=mm['sd'], fmt='s', color=C.COL['model'])
    ax.set_xticks([0, 1], ['healthy / normal', 'PD / MPTP'])
    ax.set_ylim(0, 1.05)
    ax.set_title('oscillatory fraction')
    # spectra: mean power / surrogate-median ratio per window
    rc = sorted([c for c in cells.columns if c.startswith('ratio_')], key=lambda c: float(c[6:]))
    f = [float(c[6:]) for c in rc]
    for ax, sm, sd in ((axs[1, 2], 'healthy', 'normal'), (axs[1, 3], 'parkinsonian', 'mptp')):
        ax.plot(f, cells.loc[cells.state == sm, rc].median().values, color=C.COL['model'], label='model')
        ax.plot(f, lm.loc[lm.state == sd, rc].median().values, color=C.COL['data'], label='data')
        ax.axhline(1, color='k', lw=0.5)
        ax.set_xlabel('Hz')
        ax.set_title(f'{sm}: median power / surrogate median')
        ax.legend(fontsize=6)
    fig.savefig(os.path.join(out, 'H5_summary.png'), dpi=140)


if __name__ == '__main__':
    main()
