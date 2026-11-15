"""
H4 comparison: virtual optogenetics (model, healthy state) vs mouse D2- / D1-MSN ChR2 (dopamine-intact, anesthetized,
multi-unit). Reads H4/prediction.json first.
Data targets: D2 from binned_response.csv (group D2) and d2_wc_response_xtalk_removed.csv; D1 only from
binned_response.csv (group D1; the KS4 D1 sessions have no stimulation).
Verdict per model scenario vs its pathway (fixed before the comparison):
- mismatch if the population effect has the opposite sign to the data (or no cell responds);
- match if the sign agrees, responses are one-signed (opposite-sign responders < 10% of cells) and the responder
  fraction (seed mean) is inside the 95% Wilson CI of the binned data;
- partial otherwise (sign agrees only).

    python -m validation.heldout.compare.compare_H4
"""
import os

import numpy as np
import pandas as pd

from validation.heldout.compare import common as C
from validation.heldout.metrics import response

PATHWAY = {'D2_gpe': 'D2', 'D2_ctx': 'D2', 'D1_ctx': 'D1'}
DATA_SIGN = {'D2': 'inh', 'D1': 'exc'}


def data_targets(data_io):
    b = data_io.chr2_binned()
    wc = data_io.chr2_d2_wc()
    wc['ratio'] = (1 + wc.mi) / (1 - wc.mi)
    t = {}
    for g in ('D1', 'D2'):
        s = b[b.group == g].copy()
        s['ratio'] = s.stim / s.pre
        sign = DATA_SIGN[g]
        r = s[s.resp == sign]
        t[g] = dict(source='binned', n=len(s), n_exc=int((s.resp == 'exc').sum()), n_inh=int((s.resp == 'inh').sum()),
                    frac_resp=float((s.resp == sign).mean()), ci=C.wilson(int((s.resp == sign).sum()), len(s)),
                    median_ratio_resp=float(r.ratio.median()), early=float((r.early / r.pre).median()),
                    late=float((r.late / r.pre).median()), post=float((r.post / r.pre).median()),
                    pop_ratio=float(s.stim.sum() / s.pre.sum()))
    t['D2_wc'] = dict(source='wave_clus (crosstalk removed)', n=len(wc), n_exc=int((wc.resp == 'exc').sum()),
                      n_inh=int((wc.resp == 'inh').sum()), frac_resp=float((wc.resp == 'inh').mean()),
                      median_ratio_resp=float(wc.loc[wc.resp == 'inh', 'ratio'].median()))
    return t


def main():
    pred = C.require_prediction('H4')
    from validation.heldout import data_io
    out = C.out_dir('H4')
    T = data_targets(data_io)
    M = pred['metrics']
    resp = pd.read_csv(os.path.join(out, 'model_responses.csv'))
    prof = pd.read_csv(os.path.join(out, 'model_profiles.csv'))

    rows, verdicts = [], {}
    for sc, mm in M.items():
        s = mm['summary']
        g = PATHWAY[sc]
        d = T[g]
        want = DATA_SIGN[g]
        other = 'exc' if want == 'inh' else 'inh'
        model_sign = 'exc' if s['pop_ratio']['mean'] > 1 else 'inh'
        f_want, f_other = s[f'frac_{want}']['mean'], s[f'frac_{other}']['mean']
        if f_want + f_other == 0 or model_sign != want:
            v = 'mismatch'
        elif f_other < 0.1 and d['ci'][0] <= f_want <= d['ci'][1]:
            v = 'match'
        else:
            v = 'partial'
        verdicts[sc] = v
        rt = s['responder_frac_by_type']
        rows.append(dict(scenario=sc, route=s['label'], pathway=g, verdict=v, model_pop_ratio=C.fmt(s['pop_ratio']),
                         data_pop_ratio=f"{d['pop_ratio']:.2f}", model_frac_inh=C.fmt(s['frac_inh']), model_frac_exc=C.fmt(s['frac_exc']),
                         data_inh_exc=f"{d['n_inh']}/{d['n']} inh, {d['n_exc']}/{d['n']} exc",
                         model_resp_ratio=C.fmt(s[f'median_ratio_{model_sign}']) if np.isfinite(s[f'median_ratio_{model_sign}']['mean']) else 'n/a',
                         data_resp_ratio=f"{d['median_ratio_resp']:.2f}",
                         model_onset_late_post=f"{s['onset_ratio']['mean']:.2f} / {s['late_ratio']['mean']:.2f} / {s['post_ratio']['mean']:.2f}",
                         data_early_late_post=f"{d['early']:.2f} / {d['late']:.2f} / {d['post']:.2f}",
                         resp_frac_PVp=C.fmt(rt['PV+']), resp_frac_PVm=C.fmt(rt['PV-']),
                         type_V=f"{s['type_enrichment_pooled'].get('V', np.nan):.2f}",
                         type_fisher_p=f"{s['type_enrichment_pooled'].get('fisher_p', np.nan):.3g}"))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(out, 'comparison_table.csv'), index=False)

    plt = C.style()
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    P = pred['priors']['H4_optogenetics']['protocol']
    for ax, sc in zip(axs, ['D2_gpe', 'D2_ctx', 'D1_ctx']):
        if sc not in M:
            continue
        g = PATHWAY[sc]
        ax.axvspan(0, P['on_s'], color='#f2d7a6', alpha=0.5, lw=0)
        for bw, lw in ((0.25, 1.4), (1.0, 0.8)):
            p = prof[(prof.scenario == sc) & (prof.bin_s == bw)].copy()
            p = p.merge(resp[resp.scenario == sc][['seed', 'cell', 'resp']], on=['seed', 'cell'])
            for sel, c, lab in ((p.resp.notna(), C.COL['model'], 'model, all cells'), (p.resp != 'ns', '#7b241c', 'model responders')):
                q = p[sel & (p.pre > 0)]
                if q.empty:
                    continue
                n = q.groupby(['seed', 't']).apply(lambda x: x.rate.sum() / x.pre.sum()).unstack(0)
                ax.plot(n.index, n.mean(1), color=c, lw=lw, label=f'{lab} ({bw:g} s bins)')
                if bw == 0.25:
                    ax.fill_between(n.index, n.mean(1) - n.std(1), n.mean(1) + n.std(1), color=c, alpha=0.15, lw=0)
        d = T[g]
        ax.plot([-0.5, 0.125, 1.1, P['on_s'] + 1], [1, d['early'], d['late'], d['post']], 'ks', ms=5,
                label=f'data {g} responders: pre / early / late / post (binned windows)')
        ax.axhline(1, color='grey', lw=0.5)
        ax.set_title(f"{M[sc]['summary']['label']}\n{verdicts[sc]}", fontsize=8)
        ax.set_xlabel('time from light onset (s; model protocol compressed)')
    axs[0].set_ylabel('rate / pre-stimulus rate')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=5, fontsize=6.5, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(out, 'H4_summary.png'), dpi=150)

    vs = ', '.join(f'{sc} {v}' for sc, v in verdicts.items())
    overall = 'match' if all(v == 'match' for v in verdicts.values()) else 'mismatch' if all(v == 'mismatch' for v in verdicts.values()) else 'partial'
    line = (f"H4: {overall.upper()} — {vs}. Data: D2 binned {T['D2']['n_inh']}/{T['D2']['n']} inhibited (0 excited), "
            f"responders x{T['D2']['median_ratio_resp']:.2f}; D1 {T['D1']['n_exc']}/{T['D1']['n']} excited, x{T['D1']['median_ratio_resp']:.2f}.")
    sections = [('Comparison', C.md_table(tab)),
                ('Data targets', '\n'.join(f'- {k}: ' + ', '.join(f'{a}={b if not isinstance(b, float) else round(b, 3)}' for a, b in v.items()) for k, v in T.items())),
                ('Notes', f"binned_response.csv holds {T['D1']['n']} D1 units (the task text says 60). The KS4 trial-aligned profiles are "
                 "available only as figures in the bundle, so the time course is compared with the binned windows (early / late / post; "
                 "window definitions in the ChR2 analysis) and the model protocol is compressed (1 s pre / 2 s on / 2 s post): the 20 s "
                 "recovery and 10-20 s after-effects are not tested. Route magnitudes are a priori (priors.json), not fitted."),
                ('Caveats', 'Mouse, dopamine-intact, anesthetized, multi-unit clusters (low rates ~1-4 Hz) vs a rat-based healthy model '
                 'state without anesthesia; striatal stimulation reaches the STN through routes the model represents only as input-rate changes.')]
    C.write_summary('H4', line, sections, pred)


if __name__ == '__main__':
    main()
