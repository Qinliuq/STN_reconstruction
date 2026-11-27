"""Step 1 diagnostics for the edge / interior-maximum artefact check (results/edge_artefacts). No model changes.

Rebuilds, with the same seeds as the drive-topography runs, the realised STN positions per PV+ placement
(tools.generate_locs, seeds['loc'] = 1) and the realised GPe->STN connections (gpe_feedback_v3, gpe_seed = 100 + s,
parkinsonian state, configuration D: 12 groups, homes (g + 0.5)/12, sigma 0.1, target_probability 1), and reports:
 1. positions near the edges (cells within 0.01 of 0, exact 0.00 duplicates from clipping, max, 0.05-bin histogram)
 2. n_eff at every profile position (network_metrics kernel)
 3. per-cell GPe coverage: number of GPe groups connected and total GPe->STN weight (groups x per-arrival weight,
    every group has the same baseline rate) vs axis
 4. delivered beta depth per cell (mean depth over its groups) and depth x total weight, for conditions 1 and 2 of
    variant b (depth only, M 0.8825), as kernel profiles over x, overlaid on the PD - healthy difference maps.
    python tools/edge_artefacts/step1_diagnostics.py
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
sys.path.insert(0, os.path.join(ROOT, 'tools', 'drive_topography'))
sys.path.insert(0, os.path.join(ROOT, 'results', 'drive_topography'))
import gpe_feedback_v3 as G  # noqa: E402
import network_metrics as NM  # noqa: E402
import stn_drive_topography as D  # noqa: E402
import precompute as PC  # noqa: E402

OUT = os.path.join(ROOT, 'results', 'edge_artefacts')
POS = NM.PROFILE_POSITIONS
K = NM.PROFILE_KERNEL
SEEDS = D.SEEDS
PLACEMENTS = ('gradient', 'uniform', 'monotonic', 'pvalb')
P = json.load(open(os.path.join(ROOT, 'params', 'drive_topography.json')))


def kprof(axis, v):
    w = np.exp(-(axis[None, :] - POS[:, None]) ** 2 / (2 * K ** 2))
    return (w @ v) / w.sum(1), w.sum(1)


def generator(ax, seed, shape=None, M=None):
    from cfg import cfg
    kw = dict(G.STATES['parkinsonian'])
    if shape is not None:
        kw.update(beta_depth_shape=shape, beta_depth_mean=M, beta_depth_base=0.1)
    return G.GPeInputGenerator(n_stn=len(ax), dt_s=cfg.gpe_dt, stn_axis=ax, seed=100 + seed,
                               feedback_delay_s=max(0.0, cfg.gpe_feedback_delay - cfg.gpe_dt), **kw)


def connected(gen):
    C = np.zeros((gen.n_groups, gen.n_stn), bool)
    for g, t in enumerate(gen.targets_by_gpe_group):
        C[g, t] = True
    return C


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rep = {'positions': {}, 'n_eff': {}, 'coverage': {}, 'delivered_depth': {}}
    bins = np.round(np.arange(0, 1.0001, 0.05), 2)
    lines = ['# Step 1 diagnostics (no model changes)', '']
    # 1, 2
    lines += ['## 1. Positions near the edges (fixed network, seeds[loc] = 1; 58 cells)', '',
              '| placement | cells x < 0.01 | exact 0.00 (clipped) | of which PV+ | min x | max x | cells x > 0.95 |', '|---|---|---|---|---|---|---|']
    hist_rows, neff_rows = [], []
    axes = {}
    for pl in PLACEMENTS:
        ax, npv = PC.positions(pl)
        axes[pl] = (ax, npv)
        ispv = np.arange(len(ax)) < npv
        n0 = int((ax == 0).sum())
        lines.append(f'| {pl} | {(ax < 0.01).sum()} | {n0} | {int((ax[ispv] == 0).sum())} | {ax.min():.3f} | {ax.max():.3f} | {(ax > 0.95).sum()} |')
        h, _ = np.histogram(ax, bins=bins)
        hp, _ = np.histogram(ax[ispv], bins=bins)
        hist_rows.append((pl, h, hp))
        neff = np.exp(-(ax[None, :] - POS[:, None]) ** 2 / (2 * K ** 2)).sum(1)
        neff_rows.append((pl, neff))
        rep['positions'][pl] = dict(n_lt_001=int((ax < 0.01).sum()), n_exact0=n0, max=float(ax.max()), hist_all=h.tolist(), hist_pvp=hp.tolist())
        rep['n_eff'][pl] = neff.round(2).tolist()
    lines += ['', 'Histogram, all cells (PV+ in brackets), 0.05 bins starting at 0:', '',
              '| bin | ' + ' | '.join(pl for pl, *_ in hist_rows) + ' |', '|---|' + '---|' * len(hist_rows)]
    for i in range(len(bins) - 1):
        lines.append(f'| {bins[i]:.2f}–{bins[i + 1]:.2f} | ' + ' | '.join(f'{h[i]} ({hp[i]})' for _, h, hp in hist_rows) + ' |')
    lines += ['', '## 2. n_eff (sum of kernel weights, sigma 0.12) at each profile position', '',
              '| x | ' + ' | '.join(pl for pl, _ in neff_rows) + ' |', '|---|' + '---|' * len(neff_rows)]
    for j, x in enumerate(POS):
        lines.append(f'| {x:.2f} | ' + ' | '.join(f'{n[j]:.1f}' + (' ⚠' if n[j] < 8 else '') for _, n in neff_rows) + ' |')
    lines.append('\n⚠ = n_eff < 8.')
    # 3. coverage (gradient placement = configuration D; the others for completeness)
    from cfg import cfg
    gscale = cfg.gpe_gaba_scale
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))
    lines += ['', '## 3. GPe coverage per STN cell (PD state, configuration D connectivity; mean over input seeds 1–5)', '',
              'Total GPe weight = groups connected × per-arrival weight (1.0 nS × synapse scale 2.0 × type factor 1.13 PV+ / '
              '0.87 PV−); all groups share the same baseline rate. Kernel profile (σ 0.12) over x:', '',
              '| placement | groups/cell at x=0.05 | at 0.50 | at 0.95 | expected groups (Σ affinity) at axis 0 / 0.5 / 1 |', '|---|---|---|---|---|']
    hc = G.GPeInputGenerator(n_stn=3, dt_s=0.001, stn_axis=[0, .5, 1], seed=1)
    aff = np.exp(-(np.array([0, .5, 1])[None, :] - hc.group_center[:, None]) ** 2 / (2 * 0.1 ** 2)).sum(0)
    for pl in PLACEMENTS:
        ax, npv = axes[pl]
        tw = 1.0 * 2.0 * np.where(np.arange(len(ax)) < npv, gscale['PVP'], gscale['PVN'])
        ng = np.mean([connected(generator(ax, s)).sum(0) for s in SEEDS], 0)
        pg, _ = kprof(ax, ng)
        rep['coverage'][pl] = dict(axis=ax.round(4).tolist(), groups_mean=ng.tolist(), total_weight_nS=(ng * tw).tolist(), profile_groups=pg.tolist())
        lines.append(f'| {pl} | {pg[0]:.2f} | {pg[9]:.2f} | {pg[-1]:.2f} | {aff[0]:.2f} / {aff[1]:.2f} / {aff[2]:.2f} |')
        if pl == 'gradient':
            order = np.argsort(ax)
            axs[0, 0].plot(ax, ng, 'o', ms=3, color='k', label='per cell (mean over seeds)')
            axs[0, 0].plot(POS, pg, '-', color='tab:red', label='kernel profile')
            axs[0, 0].set(title='GPe groups per STN cell (gradient placement)', xlabel='STN axis', ylabel='groups')
            axs[0, 1].plot(ax[:npv], (ng * tw)[:npv], 'o', ms=3, color='tab:green', label='PV+')
            axs[0, 1].plot(ax[npv:], (ng * tw)[npv:], 'o', ms=3, color='tab:blue', label='PV−')
            axs[0, 1].plot(POS, kprof(ax, ng * tw)[0], '-', color='tab:red', label='kernel profile')
            axs[0, 1].set(title='total GPe→STN weight per cell (nS per group-arrival)', xlabel='STN axis', ylabel='nS')
            axs[0, 0].legend(fontsize=7); axs[0, 1].legend(fontsize=7)
    for pl, n in neff_rows:
        axs[0, 2].plot(POS, n, label=pl)
    axs[0, 2].axhline(8, ls=':', c='gray'); axs[0, 2].axhline(5, ls='--', c='gray')
    axs[0, 2].set(title='n_eff per profile position', xlabel='STN axis', ylabel='n_eff'); axs[0, 2].legend(fontsize=7)
    # 4. delivered depth, conditions 1 and 2 (variant b, gradient placement, M advanced)
    import analyze as A
    _, profs = A.runs_table()
    M = P['sets']['advanced']['M']
    ax, npv = axes['gradient']
    lines += ['', '## 4. Delivered beta depth, conditions 1 and 2, variant b (depth only, M 0.8825, gradient placement)', '',
              'Per cell: delivered depth = mean depth over its connected groups; drive amplitude = Σ over its groups of depth × '
              'per-arrival weight (= delivered depth × total weight). Kernel profiles (σ 0.12), mean over seeds; argmax over '
              'positions with n_eff ≥ 5.', '',
              '| condition | delivered depth: argmax / range | depth × weight: argmax | PD − healthy map: argmax | interior (≥ 0.15) map argmax |', '|---|---|---|---|---|']
    neff = np.exp(-(ax[None, :] - POS[:, None]) ** 2 / (2 * K ** 2)).sum(1)
    ok = neff >= NM.PROFILE_MIN_CELLS
    for j, cond in enumerate(('1_uniform_control', '2_territory_adv')):
        c = D.CONDITIONS[cond]
        dd, amp = [], []
        for s in SEEDS:
            gen = generator(ax, s, shape=lambda x, kw=c['drive']: D.drive_scaling(np.asarray(x), **kw), M=M)
            C = connected(gen)
            depth_cell = (C * gen.beta_depth[:, None]).sum(0) / C.sum(0)
            tw = 1.0 * 2.0 * np.where(np.arange(len(ax)) < npv, gscale['PVP'], gscale['PVN'])
            dd.append(depth_cell); amp.append((C * gen.beta_depth[:, None]).sum(0) * tw)
        dd, amp = np.mean(dd, 0), np.mean(amp, 0)
        pdd, pamp = kprof(ax, dd)[0], kprof(ax, amp)[0]
        diff = np.mean([profs[('b', cond, s)]['diff'] for s in SEEDS], 0) * 1e6
        am = lambda p: POS[ok][np.argmax(p[ok])]
        inn = ok & (POS >= 0.15)
        rep['delivered_depth'][cond] = dict(depth_cell=dd.tolist(), amp_cell=amp.tolist(), profile_depth=pdd.tolist(), profile_amp=pamp.tolist(), diff_map=diff.tolist())
        lines.append(f'| {cond} | {am(pdd):.2f} / {dd.min():.2f}–{dd.max():.2f} | {am(pamp):.2f} | {am(diff):.2f} | {POS[inn][np.argmax(diff[inn])]:.2f} |')
        a = axs[1, j]
        a.plot(POS, diff, 'k-', lw=2, label='PD − healthy shared beta (×1e-6)')
        a2 = a.twinx()
        a2.plot(POS, pamp / pamp.mean(), '-', color='tab:red', label='depth × total weight (norm.)')
        a2.plot(POS, pdd / pdd.mean(), '--', color='tab:orange', label='delivered depth (norm.)')
        a2.plot(POS, kprof(ax, np.mean([connected(generator(ax, s)).sum(0) for s in SEEDS], 0))[0] / 3, ':', color='tab:purple', label='groups/cell ÷ 3')
        a.set(title=f'{cond}, variant b', xlabel='STN axis'); a.legend(fontsize=7, loc='upper left'); a2.legend(fontsize=7, loc='upper right')
    for i in range(len(bins) - 1):
        pass
    a = axs[1, 2]
    for pl, h, hp in hist_rows:
        a.step(bins[:-1], h, where='post', label=pl)
    a.set(title='cells per 0.05 bin (all cells)', xlabel='STN axis'); a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'step1_diagnostics.png'), dpi=120)
    json.dump(rep, open(os.path.join(OUT, 'step1_diagnostics.json'), 'w'), indent=0)
    open(os.path.join(OUT, 'step1_diagnostics.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
