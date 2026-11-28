"""
H3 comparison: spike phase-locking to cortex. Reads H3/prediction.json first, then computes the rat per-unit locking
(saved to H3/rat_phase_locking.csv) with the same phase functions as the model:
- reference: IpsiEEG (fallback ContraEEG), 1 kHz; units: the 45 distinct SS units in sessions with EEG
- slow (0.5-2 Hz): all of each unit's spikes
- beta (13-30 Hz): only spikes in activated-cortex windows: 4 s windows with slow-wave fraction
  (0.5-2 Hz / 0.5-45 Hz power) < 0.4 (rule fixed from the EEG alone before any phase was computed); beta over all
  epochs is reported as secondary
- units need >= 50 spikes in the analysed epochs
Verdict per band (fixed before the comparison): match if the model and rat PPC distributions do not differ
(Mann-Whitney p >= 0.05); partial if they differ but the medians are within a factor 3 and the same sign; else mismatch.
Overall: match if both bands match, mismatch if both mismatch, else partial.
Revised 2026-11-18 (on review; the original verdict is reported as `verdict_v1`): slow-band locking is set mainly by
the imposed input envelope (priors.json), so it is reported but no longer counts towards the verdict. The overall
verdict uses the beta band only; if beta is not testable, the verdict is 'uninformative'.

    python -m validation.heldout.compare.compare_H3
"""
import os

import numpy as np
import pandas as pd

from validation.heldout.compare import common as C
from validation.heldout.metrics import phase, stats

ACTIVATED_MAX_SLOW_FRAC = 0.4
WIN_S = 4.0


def rat_locking(data_io):
    trains = data_io.rat_ss_trains()
    clus = data_io.rat_dedup(data_io.rat_ss_table()).set_index(['file', 'label']).ward_cluster
    rows, epochs = [], []
    for (fname, label), st in trains.items():
        sess = fname[:-4]
        e = data_io.rat_eeg(sess)
        base = dict(file=fname, label=label, cluster='bursty' if clus.get((fname, label)) == 0 else 'regular')
        if e is None:
            rows.append(dict(base, eeg='none'))
            continue
        x, fs, t0, ch = e
        starts, frac = phase.slow_fraction_windows(x, fs, WIN_S)
        inside = (starts >= st.min()) & (starts + WIN_S <= st.max())
        act = frac < ACTIVATED_MAX_SLOW_FRAC
        mask = phase.window_mask(len(x), fs, starts, act)
        row = dict(base, eeg=ch, n_windows=int(inside.sum()), n_activated_windows=int((inside & act).sum()))
        for name, band, msk in (('slow', 'slow', None), ('beta_activated', 'beta', mask), ('beta_all', 'beta', None)):
            r = phase.locking(st, x, fs, band, t0, msk)
            row.update({f'{name}_{k}': v for k, v in r.items()})
        rows.append(row)
    return pd.DataFrame(rows)


def band_verdict(m, d):
    m, d = np.asarray(m, float), np.asarray(d, float)
    m, d = m[np.isfinite(m)], d[np.isfinite(d)]
    if len(d) < 3 or len(m) < 3:
        return 'not testable', None
    r = stats.compare_distributions(m, d)
    if r['mw_p'] >= 0.05:
        return 'match', r
    a, b = r['median_a'], r['median_b']
    if a * b > 0 and max(abs(a), abs(b)) <= 3 * min(abs(a), abs(b)):
        return 'partial', r
    return 'mismatch', r


def test_lines(items):
    out = []
    for n, v, r in items:
        txt = f"- {n}: {v}"
        if r:
            txt += f"; MW p {r['mw_p']:.3g}, KS p {r['ks_p']:.3g}, medians model {r['median_a']:.3g} / rat {r['median_b']:.3g}"
        out.append(txt)
    return '\n'.join(out)


def main():
    pred = C.require_prediction('H3')
    from validation.heldout import data_io
    out = C.out_dir('H3')
    rat = rat_locking(data_io)
    rat.to_csv(os.path.join(out, 'rat_phase_locking.csv'), index=False)
    ma = pd.read_csv(os.path.join(out, 'model_phase_slowosc.csv'))
    mb = pd.read_csv(os.path.join(out, 'model_phase_pd_beta.csv'))
    re = rat[rat.eeg != 'none']

    vs, rs = band_verdict(ma.slow_ctx_ppc, re.slow_ppc)
    vb, rb = band_verdict(mb.beta_ctx_ppc, re.beta_activated_ppc)
    vb_gpe, rb_gpe = band_verdict(mb.beta_gpe_all_ppc, re.beta_activated_ppc)
    # secondary (not in the verdict): beta over all epochs (slow-wave dominated) vs beta in the slow-oscillation runs
    vb2, rb2 = band_verdict(ma.beta_ctx_ppc, re.beta_all_ppc)
    vb2_gpe, rb2_gpe = band_verdict(ma.beta_gpe_all_ppc, re.beta_all_ppc)
    tested = [v for v in (vs, vb) if v != 'not testable']
    verdict_v1 = ('match' if tested and all(v == 'match' for v in tested) else
               'mismatch' if tested and all(v == 'mismatch' for v in tested) else 'partial')
    overall = vb if vb != 'not testable' else 'uninformative'

    rows = []
    def add(src, grp, s, ppc, pref, pval):
        s = s.copy()
        ok = np.isfinite(s[ppc])
        rows.append(dict(source=src, group=grp, n=int(ok.sum()), median_ppc=s.loc[ok, ppc].median(),
                         q25=s.loc[ok, ppc].quantile(0.25), q75=s.loc[ok, ppc].quantile(0.75),
                         frac_rayleigh_p01=float((s.loc[ok, pval] < 0.01).mean()) if ok.any() else np.nan,
                         circ_mean_pref_rad=float(np.angle(np.exp(1j * s.loc[ok, pref]).mean())) if ok.any() else np.nan))
    for grp, s in [('all', re)] + list(re.groupby('cluster')):
        add('rat slow', grp, s, 'slow_ppc', 'slow_pref_phase', 'slow_rayleigh_p')
        add('rat beta (activated)', grp, s, 'beta_activated_ppc', 'beta_activated_pref_phase', 'beta_activated_rayleigh_p')
        add('rat beta (all epochs)', grp, s, 'beta_all_ppc', 'beta_all_pref_phase', 'beta_all_rayleigh_p')
    for col, src, df in (('slow_ctx', 'model slow (PD+SO, cortical input)', ma), ('beta_ctx', 'model beta (cortical input)', mb),
                         ('beta_gpe_all', 'model beta (GPe drive)', mb), ('beta_gpe_burst', 'model beta (GPe drive, in bursts)', mb),
                         ('slow_ctx_control', 'model slow, no SO (control)', mb),
                         ('beta_ctx', 'model beta in PD+SO runs (cortical input)', ma), ('beta_gpe_all', 'model beta in PD+SO runs (GPe drive)', ma)):
        for grp, s in [('all', df)] + list(df.groupby('type')) + list(df.groupby(df.bursty_seedwise.map({True: 'bursty', False: 'regular'}))):
            add(src, grp, s, f'{col}_ppc', f'{col}_pref', f'{col}_rayleigh_p')
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(out, 'phase_locking_summary.csv'), index=False)

    plt = C.style()
    fig = plt.figure(figsize=(12, 6.5))
    ax = fig.add_subplot(2, 3, 1)
    for lab, v, c, ls in (('rat', re.slow_ppc, C.COL['data'], '-'), ('model PD+SO', ma.slow_ctx_ppc, C.COL['model'], '-'),
                          ('model PV+', ma.loc[ma.type == 'PV+', 'slow_ctx_ppc'], C.COL['PV+'], '--'),
                          ('model PV-', ma.loc[ma.type == 'PV-', 'slow_ctx_ppc'], C.COL['PV-'], '--')):
        x, y = stats.ecdf(v)
        ax.step(x, y, color=c, ls=ls, label=lab, where='post')
    ax.set_xscale('symlog', linthresh=1e-3)
    ax.set_title(f'slow 0.5-2 Hz PPC ({vs})')
    ax.legend(fontsize=6)
    ax = fig.add_subplot(2, 3, 2)
    for lab, v, c, ls in (('rat activated', re.beta_activated_ppc, C.COL['data'], '-'), ('rat all epochs', re.beta_all_ppc, C.COL['data'], ':'),
                          ('model, cortical input', mb.beta_ctx_ppc, C.COL['model'], '-'), ('model, GPe drive', mb.beta_gpe_all_ppc, C.COL['model'], '--')):
        x, y = stats.ecdf(v)
        if len(x):
            ax.step(x, y, color=c, ls=ls, label=lab, where='post')
    ax.set_xscale('symlog', linthresh=1e-3)
    x, y = stats.ecdf(ma.beta_ctx_ppc)
    ax.step(x, y, color=C.COL['model'], ls=':', label='model PD+SO, cortical input', where='post')
    ax.set_title(f'beta 13-30 Hz PPC (activated: {vb}; slow-wave epochs: {vb2})')
    ax.legend(fontsize=6)
    ax = fig.add_subplot(2, 3, 3)
    for lab, s, c in (('rat bursty', re[re.cluster == 'bursty'].slow_ppc, C.COL['bursty']), ('rat regular', re[re.cluster == 'regular'].slow_ppc, C.COL['regular'])):
        x, y = stats.ecdf(s)
        ax.step(x, y, color=c, label=lab, where='post')
    for lab, s, c in (('model bursty', ma[ma.bursty_seedwise == True].slow_ctx_ppc, C.COL['bursty']), ('model regular', ma[ma.bursty_seedwise == False].slow_ctx_ppc, C.COL['regular'])):
        x, y = stats.ecdf(s)
        ax.step(x, y, color=c, ls='--', label=lab, where='post')
    ax.set_xscale('symlog', linthresh=1e-3)
    ax.set_title('slow PPC by H1 cluster')
    ax.legend(fontsize=6)
    for k, (lab, th, c) in enumerate((('rat slow', re.slow_pref_phase, C.COL['data']), ('model slow (cortical input)', ma.slow_ctx_pref, C.COL['model']),
                                      ('model beta (GPe drive)', mb.beta_gpe_all_pref, C.COL['model']))):
        ax = fig.add_subplot(2, 3, 4 + k, projection='polar')
        th = np.asarray(th, float)
        th = th[np.isfinite(th)]
        ax.hist(th, bins=18, color=c, alpha=0.7)
        ax.set_title(f'preferred phase: {lab}', fontsize=8)
    fig.suptitle(f'H3 phase locking to cortex: {overall}', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'H3_summary.png'), dpi=150)

    n_beta = int(re.beta_activated_ppc.notna().sum())
    line = (f"H3: {overall.upper()} (original rule: {verdict_v1}; the slow-band result is set mainly by the imposed input envelope) — slow PPC {vs} (model median {ma.slow_ctx_ppc.median():.3g} vs rat {re.slow_ppc.median():.3g}); "
            f"beta PPC to cortical input {vb} (model {mb.beta_ctx_ppc.median():.3g} vs rat activated-epoch {re.beta_activated_ppc.median():.3g}, "
            f"{n_beta} rat units qualify); beta to GPe drive {vb_gpe} (model {mb.beta_gpe_all_ppc.median():.3g}). "
            f"Secondary, slow-wave epochs: beta {vb2} (model PD+SO {ma.beta_ctx_ppc.median():.3g} vs rat all-epoch {re.beta_all_ppc.median():.3g}).")
    sess = re.groupby('file').agg(units=('label', 'size'), windows=('n_windows', 'max'), activated=('n_activated_windows', 'max'))
    sections = [('Summary table', C.md_table(tab)),
                ('Activated-cortex epochs (rat)', f"Rule: 4 s windows with slow fraction < {ACTIVATED_MAX_SLOW_FRAC}. Per session "
                 f"(windows within the unit spike span):\n\n{C.md_table(sess.reset_index())}\n\n{n_beta} of {len(re)} units with EEG have >= 50 spikes in activated windows. "
                 f"{int((rat.eeg == 'none').sum())} units (pjx286c, pjx289c) have no EEG. The unit trains are 100 s epochs chosen during "
                 "slow-wave activity; the activated windows of NPR040d, NPR057 and pjx310c fall outside them, so the activated-epoch beta "
                 "test cannot be run on these data."),
                ('Tests', test_lines((('slow', vs, rs), ('beta (activated epochs), cortical input', vb, rb), ('beta (activated epochs), GPe drive', vb_gpe, rb_gpe),
                                       ('secondary: rat beta all epochs vs model beta in PD+SO runs, cortical input', vb2, rb2),
                                       ('secondary: same, GPe drive', vb2_gpe, rb2_gpe)))),
                ('Caveats', "The model has no slow-wave cortex: the slow oscillation is an imposed ~1 Hz up/down envelope on the cortical input "
                 "(priors.json: depth 0.9, duty 0.5, mean preserved), so slow locking is largely set by that choice. The model reference is the "
                 "cortical input rate (phase 0 = peak input) or the GPe beta drive (phase 0 = peak GPe rate); the rat reference is surface EEG, "
                 "whose polarity relative to cortical firing is not established here, so preferred phases are compared for concentration, "
                 "not absolute angle. Rat EEG is mostly slow-wave dominated; activated epochs exist in 3 of 6 sessions. Anesthesia: urethane (rat).")]
    C.write_summary('H3', line, sections, pred,
                    dict(verdict_word=overall, model=f"slow PPC {ma.slow_ctx_ppc.median():.3g}; beta PPC {mb.beta_ctx_ppc.median():.3g} (ctx) / {mb.beta_gpe_all_ppc.median():.3g} (GPe)",
                 data=f"slow PPC {re.slow_ppc.median():.3g}; beta (all epochs) {re.beta_all_ppc.median():.3g}; activated: {n_beta} units", statistic='Mann-Whitney on PPC (medians shown)'))


if __name__ == '__main__':
    main()
