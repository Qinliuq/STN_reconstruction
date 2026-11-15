"""
H1 comparison: two-population burst/pause structure (rat 6-OHDA single units vs model PD state).
Reads H1/prediction.json first (refuses to run without it), then the rat tables.
Verdict rule (fixed before the comparison): (a) bimodal = silhouette p < 0.05 vs the Gaussian null AND BIC prefers
k = 2, in >= 3 of 5 seeds; (b) bursty-cluster proportion (seed mean) inside the 95% Wilson CI of the data (33/73).
match = (a) and (b); partial = one of them; mismatch = neither.

    python -m validation.heldout.compare.compare_H1
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from validation.heldout.compare import common as C
from validation.heldout.metrics import clustering, spike_features as sf, stats


def main():
    pred = C.require_prediction('H1')
    from validation.heldout import data_io  # data are read only after the prediction gate
    out = C.out_dir('H1')
    m = pred['metrics']['primary_pd']
    model = pd.read_csv(os.path.join(C.ROOT, 'results', 'heldout', 'H1', 'model_features_pd.csv'))
    rat = data_io.rat_ss_table()
    rat['cluster'] = np.where(rat.ward_cluster == 0, 'bursty', 'regular')
    dedup = data_io.rat_dedup(rat)
    model['cluster'] = np.where(model.bursty_seedwise, 'bursty', 'regular')

    # (a) bimodality, (b) proportion
    S = m['summary']
    a = S['n_seeds_sig_p05'] >= 3 and S['n_seeds_bic_k2'] >= 3
    k, n = int((rat.cluster == 'bursty').sum()), len(rat)
    lo, hi = C.wilson(k, n)
    b = lo <= S['bursty_frac']['mean'] <= hi
    verdict = 'match' if a and b else ('partial' if a or b else 'mismatch')

    # model cells assigned to the data clusters (nearest centroid in the data's standardized space)
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(rat[sf.FEATURES].values)
    Xr = sc.transform(rat[sf.FEATURES].values)
    lab = np.where(rat.cluster == 'bursty', 0, 1)
    assigned = clustering.assign_to_reference(sc.transform(model[sf.FEATURES].values), Xr, lab)
    model['data_cluster'] = np.where(assigned == 0, 'bursty', 'regular')
    frac_data_bursty = model.groupby('seed').data_cluster.apply(lambda x: np.mean(x == 'bursty'))

    # (d) per-cluster features, model vs data (rank-sum)
    rows = []
    for cl in ('bursty', 'regular'):
        for f in sf.FEATURES:
            x, y = rat.loc[rat.cluster == cl, f], model.loc[model.cluster == cl, f]
            p = mannwhitneyu(x, y).pvalue if len(y) else np.nan
            rows.append(dict(cluster=cl, feature=f, data_median=x.median(), model_median=y.median() if len(y) else np.nan,
                             n_data=len(x), n_model=len(y), mw_p=p))
    feat = pd.DataFrame(rows)
    feat.to_csv(os.path.join(out, 'per_cluster_features.csv'), index=False)
    seeds = pd.DataFrame([dict(seed=s, n=r['sizes'][0] + r['sizes'][1], silhouette=r['silhouette'], null_p95=r['null_p95'],
                               p=r['p'], bic_k1=r['bic_k1'], bic_k2=r['bic_k2'], bursty_frac=r['bursty_frac'],
                               frac_PVp_in_bursty=r['frac_PVp_in_bursty'], cramers_V=r['vs_type'].get('V'),
                               fisher_p=r['vs_type'].get('fisher_p'))
                          for s, r in m['per_seed'].items()])
    seeds.to_csv(os.path.join(out, 'per_seed_clustering.csv'), index=False)
    ver = os.path.join(out, 'rat_verification.json')
    rv = pd.read_json(ver) if os.path.exists(ver) else None

    # figure
    plt = C.style()
    fig, axs = plt.subplots(2, 5, figsize=(13, 5.6))
    for ax, f in zip(axs.flat[:8], sf.FEATURES):
        data = [rat.loc[rat.cluster == 'bursty', f], rat.loc[rat.cluster == 'regular', f],
                model.loc[model.cluster == 'bursty', f], model.loc[model.cluster == 'regular', f]]
        pos = [0, 1, 2.5, 3.5]
        for x, d, c in zip(pos, data, ['bursty', 'regular'] * 2):
            if len(d):
                ax.boxplot([d], positions=[x], widths=0.6, showfliers=False, medianprops=dict(color=C.COL[c]))
                ax.plot(x + np.random.default_rng(0).uniform(-0.2, 0.2, len(d)), d, '.', ms=2, color=C.COL[c], alpha=0.5)
        ax.set_xticks(pos, ['rat\nbursty', 'rat\nregular', 'model\nbursty', 'model\nregular'])
        ax.set_title(f)
    ax = axs.flat[8]
    ax.bar(['rat\n73', 'rat\ndedup 45', 'model\nseed mean', 'model\nPV+ frac.'],
           [k / n, (dedup.cluster == 'bursty').mean(), S['bursty_frac']['mean'], pred['metrics']['pv_fraction']],
           yerr=[0, 0, S['bursty_frac']['sd'], 0], color=[C.COL['data'], C.COL['data'], C.COL['model'], C.COL['PV+']])
    ax.axhspan(lo, hi, color='grey', alpha=0.15)
    ax.set_title('bursty-cluster proportion')
    ax = axs.flat[9]
    ax.bar(seeds.seed.astype(str), seeds.silhouette, color=C.COL['model'], label='model silhouette')
    ax.plot(seeds.seed.astype(str), seeds.null_p95, 'k_', ms=14, label='null p95')
    ax.axhline(0.367, color=C.COL['data'], ls='--', label='rat silhouette')
    ax.set_title('silhouette per seed')
    ax.set_xlabel('seed')
    ax.legend(fontsize=6)
    fig.suptitle(f'H1 burst/pause structure in the PD state: {verdict}', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'H1_summary.png'), dpi=150)

    # contingency (pooled over seeds)
    ct = pd.crosstab(model.cluster, model.type)
    sec = pred['metrics'].get('secondary_pd_slowosc')
    sec_txt = 'not run'
    if sec:
        s2 = sec['summary']
        sec_txt = (f"PD + slow oscillation (labelled secondary): seeds with p < 0.05: {s2['n_seeds_sig_p05']}/5, BIC prefers k=2: "
                   f"{s2['n_seeds_bic_k2']}/5, silhouette {C.fmt(s2['silhouette'])}, bursty fraction {C.fmt(s2['bursty_frac'])}, "
                   f"PV+ share of the bursty cluster {C.fmt(s2['frac_PVp_in_bursty'])}, Cramer's V {C.fmt(s2['V'])}.")
    line = (f"H1: {verdict.upper()} — model PD bimodal in {S['n_seeds_sig_p05']}/5 seeds (silhouette p<0.05), BIC k=2 in "
            f"{S['n_seeds_bic_k2']}/5; bursty fraction {C.fmt(S['bursty_frac'])} vs data {k}/{n} = {k/n:.2f} "
            f"(95% CI {lo:.2f}-{hi:.2f}); bursty cluster vs PV+: Cramer's V {C.fmt(S['V'])}.")
    sections = [
        ('(a) Bimodality', f"Per seed (58 cells; seeds share the network):\n\n{C.md_table(seeds)}\n\nMean ± SD over seeds: silhouette "
         f"{C.fmt(S['silhouette'])}, null p95 {C.fmt(S['null_p95'])}, p {C.fmt(S['p'], 3)}, BIC k1 {C.fmt(S['bic_k1'], 1)} / k2 "
         f"{C.fmt(S['bic_k2'], 1)}. Pooled over seeds (290 cell-runs, not independent): silhouette {m['pooled']['silhouette']:.3f}, "
         f"p {m['pooled']['p']:.3f}, BIC k1 {m['pooled']['bic_k1']:.1f} / k2 {m['pooled']['bic_k2']:.1f}.\n\n"
         f"Data (port of cluster_test.py, labels reproduced exactly): silhouette 0.367, null p95 0.321, p 0.005 here "
         f"(reference 0.316 / 0.002; the null draws depend on the NumPy version); BIC k1 582.7 / k2 555.7 (reference 583 / 556). "
         f"With the 28 duplicated unit copies removed (45 units): silhouette 0.341, p 0.047, BIC k1 365.4 / k2 356.5, 20/25 split."),
        ('(b) Bursty proportion', f"Model {C.fmt(S['bursty_frac'])} (own clusters, seed mean); model cells assigned to the nearest *data* "
         f"cluster centroid: {frac_data_bursty.mean():.2f} ± {frac_data_bursty.std(ddof=1):.2f}. Data {k}/{n} = {k/n:.2f} "
         f"(dedup {int((dedup.cluster == 'bursty').sum())}/{len(dedup)}). Model PV+ fraction {pred['metrics']['pv_fraction']:.2f}."),
        ('(c) Bursty cluster vs PV+', f"Pooled cell-runs (rows: model cluster, columns: type):\n\n{C.md_table(ct.reset_index())}\n\n"
         f"Per seed: share of PV+ cells in the bursty cluster {C.fmt(S['frac_PVp_in_bursty'])}, share of PV+ cells that are bursty "
         f"{C.fmt(S['frac_bursty_in_PVp'])}, Cramér's V {C.fmt(S['V'])}."),
        ('(d) Per-cluster features vs data', f"Rank-sum tests, rat cluster vs the model cluster of the same name (model cells pooled over seeds):\n\n{C.md_table(feat)}"),
        ('Secondary: anesthesia-like slow oscillation', sec_txt),
        ('Caveats', "Rat: urethane-anesthetized 6-OHDA, ~100 s trains; model: 20 s per seed, unanesthetized PD regime in the primary test. "
         "Model seeds share one 58-cell network (only input noise differs). 28 of the 73 rat rows are copies of 28 other units (same title "
         "exported twice; the spike times of one copy are a subset of the other), so the effective n is 45. burst_share has length-dependent "
         "false positives (Poisson ≈ 0.10 at 300 s, 20 Hz). log_isi_bimodality follows the reference code (Pearson kurtosis in Sarle's "
         "denominator; values are about half of the usual coefficient)."),
    ]
    C.write_summary('H1', line, sections, pred)


if __name__ == '__main__':
    main()
