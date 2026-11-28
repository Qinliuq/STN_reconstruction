"""Rate-normalised reanalysis of the stored spike data (no new simulations), 2026-11-28.

For every run x condition of the edge-artefact set (R0-R3), Table 5 (configuration D: Dref healthy / graded / uniform
drive) and the 2x2s that can be rebuilt from stored runs:
  1. firing-rate profile along the axis, per cell type (kernel 0.12; n_eff per type);
  2. rate-normalised pairwise beta: each cell's 1 ms spike-count train is divided by its own mean count, so the mean
     pairwise beta cross-power becomes C_ij / (nu_i nu_j) (dimensionless; shared beta modulation per spike, squared).
     Profiled with the same kernel and the same 'weighted mean minus auto terms' estimator as network_metrics; per-type
     pair means (PV+ pairs, PV- pairs, PV+-PV- pairs) as in Table 5. Raw cross-power is recomputed alongside as a check;
  3. PPC magnitude profile (across-episode PPC_1, pooled over seeds; magnitude, not fraction locked);
  4. edge diagnostics for the uniform-PV+ 0.05 maximum (rates, PV+ fraction and cortical drive at the ventral edge).
Outputs: rate_normalised.csv, rate_normalised_pertype.csv, rate_normalised_profiles.json, rate_normalised.png
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import welch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze as A  # noqa: E402

NM = A.NM
POS = np.array(A.POS)
K = NM.PROFILE_KERNEL
FLAG = 8.0
BETA = (13, 30)
OUT = os.path.join(A.ROOT, 'sweeps', 'out')
CTX = (0.65, -0.5)  # cfg.ctx_topography (unchanged in every run used here)


def per_seed(d):
    t0, t1 = float(d['t0']), float(d['t1'])
    T = (t1 - t0) / 1000
    gids = [g for pop in NM.STN_POPS for g in d[f'gids_{pop}']]
    nP = len(d['gids_PVP_pop'])
    spk, spkt = d['spkid'], d['spkt']
    bins = np.arange(t0, t1 + 1, 1.0)
    X = np.array([np.histogram(spkt[spk == g], bins=bins)[0] for g in gids], float)
    mu = X.mean(1)
    act = mu > 0
    Y = np.where(act[:, None], X / np.where(act, mu, 1)[:, None], 0.0)
    nper = min(512, X.shape[1])

    def bandp(x):
        f, P = welch(x - x.mean(-1, keepdims=True), fs=1000, nperseg=nper, axis=-1)
        return P[..., (f >= BETA[0]) & (f <= BETA[1])].sum(-1)
    autoY, autoX = bandp(Y), bandp(X)
    axis = np.asarray(d['axis'])
    isP = np.arange(len(gids)) < nP

    def pair_mean(Z, auto, w):  # weighted mean pairwise band cross-power, w over cells (zeros allowed)
        if (w > 0).sum() < 2:
            return np.nan
        wn = w / w.sum()
        return float((bandp(wn @ Z) - (wn ** 2) @ auto) / (1 - (wn ** 2).sum()))
    out = dict(rate=mu * 1000, axis=axis, isP=isP)
    # per type (equal weights = mean over pairs); cross from the all-cell sum
    for nm, Z, au in (('norm', Y, autoY), ('raw', X, autoX)):
        for typ, m in (('PVP', isP & act), ('PVN', ~isP & act)):
            out[f'{nm}_pair_{typ}'] = pair_mean(Z, au, m.astype(float))
        nPa, nNa = (isP & act).sum(), (~isP & act).sum()
        sP, sN = bandp(Z[isP & act].sum(0)), bandp(Z[~isP & act].sum(0))
        sA = bandp(Z[act].sum(0))
        out[f'{nm}_pair_cross'] = float((sA - sP - sN) / (2 * nPa * nNa))
        out[f'{nm}_pair_all'] = pair_mean(Z, au, act.astype(float))
        out[f'{nm}_prof'] = np.array([pair_mean(Z, au, np.exp(-(axis - x) ** 2 / (2 * K ** 2)) * act) for x in POS])
    W = np.exp(-(axis[None, :] - POS[:, None]) ** 2 / (2 * K ** 2))
    out['neff'] = (W * act).sum(1)
    for typ, m in (('all', np.ones_like(isP)), ('PVP', isP), ('PVN', ~isP)):
        Wm = W * m
        out[f'rate_prof_{typ}'] = (Wm @ out['rate']) / Wm.sum(1)
        out[f'rate_neff_{typ}'] = Wm.sum(1)
    out['rate2_prof'] = (W @ out['rate'] ** 2) / W.sum(1)  # kernel mean of nu^2 (what raw cross-power scales with)
    out['pvp_frac_prof'] = (W @ isP) / W.sum(1)
    return out


def load_set(d, lab):
    S = A.load(d, lab)
    return None if S is None else [x for _, x in S]


def summarise(runs):
    seeds = [per_seed(d) for d in runs]
    mean = lambda k: np.mean([s[k] for s in seeds], 0)
    sd = lambda k: np.std([s[k] for s in seeds], 0, ddof=1)
    res = {k: mean(k) for k in seeds[0] if k not in ('axis', 'isP', 'rate')}
    res.update({k + '_sd': sd(k) for k in ('norm_pair_PVP', 'norm_pair_PVN', 'norm_pair_cross', 'raw_pair_PVP', 'raw_pair_PVN', 'raw_pair_cross')})
    res['rate_PVP'] = float(np.mean([s['rate'][s['isP']].mean() for s in seeds]))
    res['rate_PVN'] = float(np.mean([s['rate'][~s['isP']].mean() for s in seeds]))
    if 'gpe_beta_phase' in runs[0]:
        q = NM.beta_ppc(runs, n_shuffle=0)
        res['ppc_prof'], res['ppc_neff'] = np.array(q['profile_ppc']), np.array(q['profile_ppc_n_eff'])
        res['ppc_PVP'], res['ppc_PVN'] = q['ppc_PVP'], q['ppc_PVN']
    a0 = runs[0]['axis']
    res['edge_cells'] = int((a0 < 0.15).sum())
    res['edge_pvp'] = int((a0[: len(runs[0]['gids_PVP_pop'])] < 0.15).sum())
    sc = CTX[0] + CTX[1] * a0
    res['ctx_scale_edge'] = float((sc / sc.mean())[a0 < 0.15].mean())
    return res


def argmax(prof, neff, lo=0.0):
    p = np.where((neff >= NM.PROFILE_MIN_CELLS) & (POS >= lo), prof, -np.inf)
    i = int(np.argmax(p))
    return float(POS[i]), float(neff[i])


def ratio(p):
    a, b = p[list(POS).index(0.65)], p[list(POS).index(0.25)]
    return float(a / b) if b > 0 else np.nan


SETS = []  # (group, run, condition, pd (dir, label), healthy (dir, label))
for run in ['R0', 'R1', 'R2', 'R3']:
    for cond in A.CONDS:
        SETS.append(('edge', run, cond, A.pd_dir(run, cond), A.h_dir(run, cond)))
FU, DT = A.FU, A.DT
SETS += [('table5', 'cfgD', 'graded_drive_PVgradient', (FU, 'Dref_pd'), (FU, 'Dref_healthy')),
         ('table5', 'cfgD', 'uniform_drive_PVgradient', (FU, 'Dref_uni_pd'), (FU, 'Dref_healthy'))]
R3_BREACH = {'1_uniform_control', '2_territory_adv', '6_pv_monotonic', '2b_territory_adv_pvalb'}


def main():
    cache, rows, prow, profs = {}, [], [], {}

    def get(src):
        if src not in cache:
            r = load_set(*src)
            cache[src] = None if r is None else summarise(r)
        return cache[src]
    # Table 5 healthy row
    h5 = get((FU, 'Dref_healthy'))
    prow.append(dict(group='table5', run='cfgD', condition='healthy', state='healthy',
                     **{k: h5[k] for k in ('norm_pair_PVP', 'norm_pair_PVN', 'norm_pair_cross', 'norm_pair_PVP_sd', 'norm_pair_PVN_sd',
                                           'raw_pair_PVP', 'raw_pair_PVN', 'raw_pair_cross', 'rate_PVP', 'rate_PVN')}))
    for grp, run, cond, pdsrc, hsrc in SETS:
        P, H = get(pdsrc), get(hsrc)
        if P is None or H is None:
            continue
        neff = P['neff']
        dn, dr = P['norm_prof'] - H['norm_prof'], P['raw_prof'] - H['raw_prof']
        row = dict(group=grp, run=run, condition=cond, r3_healthy_breach=(run == 'R3' and cond in R3_BREACH))
        row['raw_diff_max'], row['raw_neff_at_max'] = argmax(dr, neff)
        row['norm_diff_max'], row['norm_neff_at_max'] = argmax(dn, neff)
        row['norm_diff_interior_max'], row['norm_neff_at_interior'] = argmax(dn, neff, 0.15)
        row['norm_flag_lt8'] = row['norm_neff_at_max'] < FLAG
        row['norm_diff_ratio_065_025'] = ratio(dn)
        row['norm_pd_ratio_065_025'] = ratio(P['norm_prof'])
        row['raw_diff_ratio_065_025'] = ratio(dr)
        row['norm_diff_edge_005'] = float(dn[0])
        row['norm_diff_peak_value'] = float(dn.max())
        row['norm_diff_contrast'] = float((dn.max() - dn.min()) / abs(dn.mean())) if dn.mean() != 0 else np.nan
        row['norm_diff_mean'] = float(dn.mean())
        if 'ppc_prof' in P and 'ppc_prof' in H:
            dp = P['ppc_prof'] - np.nan_to_num(H['ppc_prof'])
            row['ppc_diff_max'], row['ppc_neff_at_max'] = argmax(dp, P['ppc_neff'])
            row['ppc_pd_max'], _ = argmax(P['ppc_prof'], P['ppc_neff'])
            row['ppc_diff_ratio_065_025'] = ratio(dp)
            row['ppc_PVP'], row['ppc_PVN'] = P['ppc_PVP'], P['ppc_PVN']
        else:
            dp = None
        # rate story
        row['rate_max_PD'], _ = argmax(P['rate_prof_all'], neff)
        row['rate2_diff_max'], _ = argmax(P['rate2_prof'] - H['rate2_prof'], neff)
        row['corr_rawdiff_vs_rate2diff'] = float(np.corrcoef(dr, P['rate2_prof'] - H['rate2_prof'])[0, 1])
        row['corr_rawpd_vs_rate2pd'] = float(np.corrcoef(P['raw_prof'], P['rate2_prof'])[0, 1])
        row['corr_normdiff_vs_rawdiff'] = float(np.corrcoef(dn, dr)[0, 1])
        row['rate_edge_PD'], row['rate_edge_H'] = float(P['rate_prof_all'][0]), float(H['rate_prof_all'][0])
        row['rate_mid_PD'] = float(P['rate_prof_all'][list(POS).index(0.5)])
        row['pvp_frac_edge'] = float(P['pvp_frac_prof'][0])
        row['edge_cells'], row['edge_pvp'] = P['edge_cells'], P['edge_pvp']
        row['ctx_scale_edge'] = P['ctx_scale_edge']
        row['rate_PVP'], row['rate_PVN'] = P['rate_PVP'], P['rate_PVN']
        rows.append(row)
        prow.append(dict(group=grp, run=run, condition=cond, state='PD',
                         **{k: P[k] for k in ('norm_pair_PVP', 'norm_pair_PVN', 'norm_pair_cross', 'norm_pair_PVP_sd', 'norm_pair_PVN_sd',
                                              'raw_pair_PVP', 'raw_pair_PVN', 'raw_pair_cross', 'rate_PVP', 'rate_PVN')},
                         norm_PVP_over_PVN=P['norm_pair_PVP'] / P['norm_pair_PVN'], raw_PVP_over_PVN=P['raw_pair_PVP'] / P['raw_pair_PVN'],
                         h_norm_pair_PVP=H['norm_pair_PVP'], h_norm_pair_PVN=H['norm_pair_PVN'],
                         norm_PVP_PD_over_H=P['norm_pair_PVP'] / H['norm_pair_PVP'], norm_PVN_PD_over_H=P['norm_pair_PVN'] / H['norm_pair_PVN']))
        profs[f'{run}|{cond}'] = dict(neff=neff, norm_pd=P['norm_prof'], norm_h=H['norm_prof'], norm_diff=dn, raw_diff=dr,
                                      rate_PVP=P['rate_prof_PVP'], rate_PVN=P['rate_prof_PVN'], rate_all=P['rate_prof_all'],
                                      rate_h=H['rate_prof_all'], rate_neff_PVP=P['rate_neff_PVP'], rate_neff_PVN=P['rate_neff_PVN'],
                                      ppc_pd=P.get('ppc_prof'), ppc_h=H.get('ppc_prof'), ppc_diff=dp, ppc_neff=P.get('ppc_neff'))
    df, pt = pd.DataFrame(rows), pd.DataFrame(prow)
    df.to_csv(os.path.join(HERE, 'rate_normalised.csv'), index=False, float_format='%.4g')
    pt.to_csv(os.path.join(HERE, 'rate_normalised_pertype.csv'), index=False, float_format='%.4g')
    ser = lambda v: None if v is None else [float(x) for x in np.asarray(v)]
    json.dump({'positions': POS.tolist(), **{k: {kk: ser(vv) for kk, vv in v.items()} for k, v in profs.items()}},
              open(os.path.join(HERE, 'rate_normalised_profiles.json'), 'w'))
    figure(profs)
    pd.set_option('display.width', 250)
    print(df.to_string())
    print(pt.to_string())


def figure(profs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cols = ['1_uniform_control', '2_territory_adv', '3_pv_uniform', 'geo']
    runs = ['R0', 'R1', 'R2', 'R3']
    fig, ax = plt.subplots(len(runs) + 1, len(cols), figsize=(18, 17), sharex=True)
    for i, run in enumerate(runs):
        for j, c in enumerate(cols):
            a = ax[i, j]
            p = profs.get(f'{run}|{c}')
            if p is None:
                a.text(0.5, 0.5, 'not run', ha='center', transform=a.transAxes, color='0.5')
                continue
            lo = p['neff'] < FLAG
            a.plot(POS, p['norm_diff'] * 1e3, 'k-', lw=1.8, label='rate-norm. PD − healthy (×1e-3)')
            a.plot(POS[lo], p['norm_diff'][lo] * 1e3, 'o', color='tab:red', ms=3, label='n_eff < 8')
            b = a.twinx()
            if p['ppc_diff'] is not None:
                b.plot(POS, p['ppc_diff'], '--', color='tab:blue', lw=1.4, label='PPC PD − healthy')
            b.tick_params(labelsize=6, colors='tab:blue')
            a.axvspan(0.55, 0.8, color='0.92', zorder=-1)
            if run == 'R3' and c in R3_BREACH:
                a.set_title('R3: healthy-burst breach', fontsize=7, color='tab:red')
            if i == 0:
                a.set_title(c if not (run == 'R3' and c in R3_BREACH) else a.get_title(), fontsize=9)
            if j == 0:
                a.set_ylabel(f'{run}: {A.FIXES[run]}')
            if i == 0 and j == 0:
                a.legend(fontsize=6, loc='upper left'); b.legend(fontsize=6, loc='upper right')
    for j, c in enumerate(cols):  # bottom row: R0 rate profiles by type
        a = ax[-1, j]
        p = profs[f'R0|{c}']
        a.plot(POS, p['rate_PVP'], '-', color='tab:orange', label='PV+ PD')
        a.plot(POS, p['rate_PVN'], '-', color='tab:green', label='PV− PD')
        a.plot(POS, p['rate_h'], ':', color='k', label='all, healthy')
        a.plot(POS[p['rate_neff_PVP'] < FLAG], p['rate_PVP'][p['rate_neff_PVP'] < FLAG], 'o', color='tab:red', ms=2)
        a.set_xlabel('STN axis (0 ventromedial, 1 dorsolateral)')
        if j == 0:
            a.set_ylabel('R0 rate (Hz)'); a.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'rate_normalised.png'), dpi=110)


if __name__ == '__main__':
    main()
