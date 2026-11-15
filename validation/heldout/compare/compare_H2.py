"""
H2 (pooled check, not fitting): model PD per-cell rate / CV / burst share vs the rat single units as distributions
(Mann-Whitney, KS, ECDFs). Per feature: 'similar' if both tests give p >= 0.05, else 'differs'. The test line counts
the similar features. Reads H2/prediction.json first.

    python -m validation.heldout.compare.compare_H2
"""
import os

import numpy as np
import pandas as pd

from validation.heldout.compare import common as C
from validation.heldout.metrics import stats

FEATS = [('rate_hz', 'rate (Hz)'), ('cv', 'ISI CV'), ('burst_share', 'Poisson-surprise burst share')]


def main():
    pred = C.require_prediction('H2')
    from validation.heldout import data_io
    out = C.out_dir('H2')
    pc = pred['metrics']['per_cell']
    rat = data_io.rat_ss_table()
    dedup = data_io.rat_dedup(rat)
    rows = []
    for f, _ in FEATS:
        for name, d in (('rat_73', rat), ('rat_dedup_45', dedup)):
            for grp in ('all', 'PV+', 'PV-'):
                r = stats.compare_distributions(pc[f][grp], d[f].values)
                rows.append(dict(feature=f, data=name, model_cells=grp, **r))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(out, 'distribution_tests.csv'), index=False)
    main_rows = tab[(tab.data == 'rat_73') & (tab.model_cells == 'all')]
    sim = [(r.feature, (r.mw_p >= 0.05) and (r.ks_p >= 0.05)) for r in main_rows.itertuples()]

    plt = C.style()
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
    for ax, (f, lab) in zip(axs, FEATS):
        x, y = stats.ecdf(rat[f].values)
        ax.step(x, y, color=C.COL['data'], label='rat SS (73)', where='post')
        x, y = stats.ecdf(dedup[f].values)
        ax.step(x, y, color=C.COL['data'], ls=':', label='rat dedup (45)', where='post')
        for grp, c, ls in (('all', C.COL['model'], '-'), ('PV+', C.COL['PV+'], '--'), ('PV-', C.COL['PV-'], '--')):
            x, y = stats.ecdf(pc[f][grp])
            ax.step(x, y, color=c, ls=ls, label=f'model {grp}', where='post')
        if f == 'burst_share':
            ax.set_xscale('symlog', linthresh=0.01)
        ax.set_xlabel(lab)
        ax.set_ylabel('ECDF')
    axs[0].legend(fontsize=6)
    fig.suptitle('H2 pooled check: model PD cells vs rat single units', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'H2_summary.png'), dpi=150)
    nsim = sum(s for _, s in sim)
    verdict = 'match' if nsim == 3 else ('partial' if nsim else 'mismatch')
    med = {f: (np.median(pc[f]['all']), rat[f].median()) for f, _ in FEATS}
    line = (f"H2 (pooled check): {verdict.upper()} — similar distributions for {nsim}/3 features; medians model vs rat: "
            + '; '.join(f'{f} {a:.3g} vs {b:.3g}' for f, (a, b) in med.items()) + '.')
    C.write_summary('H2', line, [('Tests', C.md_table(tab)),
                                 ('Note', 'Pooled check only: the model was not tuned to these distributions and will not be. '
                                  'Model cells are 5 seeds x 58 cells of one network (not independent); rat units are urethane-anesthetized.')], pred,
                    dict(verdict_word=verdict, model='; '.join(f'{f} {a:.3g}' for f, (a, b) in med.items()),
                 data='; '.join(f'{f} {b:.3g}' for f, (a, b) in med.items()), statistic='Mann-Whitney + KS (medians shown)'))


if __name__ == '__main__':
    main()
