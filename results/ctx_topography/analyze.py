"""Cortical->STN topography test (sweeps/ctx_topography.py): T0-T3 x pallidal depth {graded, uniform}, F1 + F2.

Per variant x depth (PD) against the variant's healthy reference, input seeds 1-5:
  - rate profiles by cell type (kernel 0.12; per-type n_eff);
  - rate-normalised pairwise beta (C_ij / nu_i nu_j; results/edge_artefacts/rate_normalised.py) PD - healthy map;
  - PPC magnitude (across-episode PPC_1, pooled over seeds) PD - healthy map;
  - absolute modulated output nu*m: kernel-mean rate x m, m = RMS shared fractional beta-band modulation
    = sqrt(df x normalised pairwise beta power), df = Welch bin width (Hz). nu*m is in Hz (shared beta-band rate
    modulation of the local population);
  - peak positions with n_eff, seed noise (peak - profile mean vs SEM across seeds), per-seed peaks;
  - operating point: mean PV+ / PV- rate vs T0 (same state and depth); stop rule 10%;
  - depth contrast: [diff(0.65) - diff(0.25)] graded minus uniform, per seed.
Outputs: summary.csv, rates.csv, profiles.json, maps.png
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'results', 'edge_artefacts'))
import rate_normalised as RN  # noqa: E402

NM, A = RN.NM, RN.A
POS = RN.POS
OUT = os.path.join(ROOT, 'sweeps', 'out', 'ctx_topography')
VARIANTS = ['T0', 'T1', 'T2', 'T3']
DEPTHS = ['graded', 'uniform']
SEEDS = [int(s) for s in os.environ.get('CTX_TOPO_SEEDS', '1,2,3,4,5').split(',')]
DF = 1000 / 512  # Welch bin width (Hz), nperseg 512 at 1 kHz
FLAG = 8.0
I65, I25 = list(POS).index(0.65), list(POS).index(0.25)


def load(lab):
    out = []
    for s in SEEDS:
        f = os.path.join(OUT, f'{lab}_s{s}')
        if not os.path.exists(os.path.join(f, 'model_output_metrics.json')):
            return None
        out.append(dict(np.load(os.path.join(f, 'model_output_spikes.npz'))))
    return out


def argmax(prof, neff, lo=0.0):
    p = np.where((np.asarray(neff) >= NM.PROFILE_MIN_CELLS) & (POS >= lo), prof, -np.inf)
    i = int(np.argmax(p))
    return float(POS[i]), float(neff[i])


def nu_m(s):  # per seed: kernel-mean rate x RMS shared fractional beta modulation (Hz)
    return s['rate_prof_all'] * np.sqrt(np.clip(s['norm_prof'], 0, None) * DF)


def stats(D, neff):
    m = D.mean(0)
    sem = D.std(0, ddof=1) / np.sqrt(len(D)) if len(D) > 1 else np.full(len(m), np.nan)
    px, pn = argmax(m, neff)
    i = list(POS).index(px)
    seedmax = [argmax(d, neff)[0] for d in D]
    return dict(max=px, neff=pn, flag=pn < FLAG, peak_minus_mean=float(m[i] - m.mean()), sem_at_peak=float(sem[i]),
                seed_max=seedmax, ratio=float(m[I65] / m[I25]) if m[I25] > 0 else np.nan), m


def main():
    seeds, ppc = {}, {}
    for v in VARIANTS:
        for lab in [f'{v}_healthy'] + [f'{v}_{d}_pd' for d in DEPTHS]:
            r = load(lab)
            if r is None:
                continue
            seeds[lab] = [RN.per_seed(d) for d in r]
            q = NM.beta_ppc(r, n_shuffle=0)
            ppc[lab] = (np.array(q['profile_ppc']), np.array(q['profile_ppc_n_eff']), q['ppc_PVP'], q['ppc_PVN'])
    rates, rows, profs = [], [], {}
    for lab, S in seeds.items():
        v = lab.split('_')[0]
        rates.append(dict(label=lab, variant=v, state='healthy' if lab.endswith('healthy') else 'PD',
                          depth=None if lab.endswith('healthy') else lab.split('_')[1],
                          rate_PVP=float(np.mean([s['rate'][s['isP']].mean() for s in S])),
                          rate_PVN=float(np.mean([s['rate'][~s['isP']].mean() for s in S]))))
    R = pd.DataFrame(rates)
    for k in ('rate_PVP', 'rate_PVN'):
        R[f'{k}_vs_T0_pct'] = [100 * (r[k] / R[(R.variant == 'T0') & (R.state == r.state) & (R.depth.fillna('') == (r.depth or ''))][k].iloc[0] - 1)
                               if ((R.variant == 'T0') & (R.state == r.state) & (R.depth.fillna('') == (r.depth or ''))).any() else np.nan
                               for _, r in R.iterrows()]
    R['breach_10pct'] = (R[['rate_PVP_vs_T0_pct', 'rate_PVN_vs_T0_pct']].abs() > 10).any(axis=1)
    R.to_csv(os.path.join(HERE, 'rates.csv'), index=False, float_format='%.4g')
    for v in VARIANTS:
        H = seeds.get(f'{v}_healthy')
        if H is None:
            continue
        hn = np.mean([s['norm_prof'] for s in H], 0)
        hm = np.mean([nu_m(s) for s in H], 0)
        for d in DEPTHS:
            P = seeds.get(f'{v}_{d}_pd')
            if P is None:
                continue
            neff = np.mean([s['neff'] for s in P], 0)
            Dn = np.array([s['norm_prof'] - hn for s in P])
            Dm = np.array([nu_m(s) - hm for s in P])
            Am = np.array([nu_m(s) for s in P])
            sn, mn = stats(Dn, neff)
            sm, mm = stats(Dm, neff)
            sa, ma = stats(Am, neff)
            pp, pneff, pPVP, pPVN = ppc[f'{v}_{d}_pd']
            ph = np.nan_to_num(ppc[f'{v}_healthy'][0])
            px, pn = argmax(pp - ph, pneff)
            row = dict(variant=v, depth=d)
            for nm, s in (('norm_diff', sn), ('num_diff', sm), ('num_pd', sa)):
                row.update({f'{nm}_{k}': (str(val) if k == 'seed_max' else val) for k, val in s.items()})
            row.update(ppc_diff_max=px, ppc_neff=pn, ppc_flag=pn < FLAG, ppc_PVP=pPVP, ppc_PVN=pPVN,
                       contrast_norm=[float(x) for x in (Dn[:, I65] - Dn[:, I25])], contrast_num=[float(x) for x in (Dm[:, I65] - Dm[:, I25])])
            rows.append(row)
            ratepd = lambda k: np.mean([s[k] for s in P], 0)
            profs[f'{v}|{d}'] = dict(neff=neff, norm_diff=mn, num_diff=mm, num_pd=ma, num_h=hm, ppc_diff=pp - ph, ppc_neff=pneff,
                                     rate_PVP=ratepd('rate_prof_PVP'), rate_PVN=ratepd('rate_prof_PVN'),
                                     rate_neff_PVP=ratepd('rate_neff_PVP'), rate_neff_PVN=ratepd('rate_neff_PVN'),
                                     rate_PVP_h=np.mean([s['rate_prof_PVP'] for s in H], 0), rate_PVN_h=np.mean([s['rate_prof_PVN'] for s in H], 0))
    df = pd.DataFrame(rows)
    # depth contrast per variant: graded minus uniform, per seed
    for k in ('norm', 'num'):
        c = []
        for _, r in df.iterrows():
            o = df[(df.variant == r.variant) & (df.depth == 'uniform')]
            if r.depth == 'graded' and len(o) and len(SEEDS) > 1:
                g, u = np.array(r[f'contrast_{k}']), np.array(o[f'contrast_{k}'].iloc[0])
                c.append(f'{(g.mean() - u.mean()):.4g} +- {np.sqrt(g.var(ddof=1) / len(g) + u.var(ddof=1) / len(u)):.3g}')
            else:
                c.append('')
        df[f'depth_contrast_{k}_graded_minus_uniform'] = c
    df.to_csv(os.path.join(HERE, 'summary.csv'), index=False, float_format='%.4g')
    ser = lambda a: [float(x) for x in np.asarray(a)]
    json.dump({'positions': POS.tolist(), **{k: {kk: ser(vv) for kk, vv in p.items()} for k, p in profs.items()}},
              open(os.path.join(HERE, 'profiles.json'), 'w'))
    figure(profs)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 80)
    print(R.round(2).to_string())
    print(df.drop(columns=['contrast_norm', 'contrast_num']).round(4).to_string())


def figure(profs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(4, 4, figsize=(18, 15), sharex=True)
    for j, v in enumerate(VARIANTS):
        for i, (title, fn) in enumerate([('rate (Hz): PV+ / PV−, PD solid, healthy dotted', None),
                                         ('rate-normalised beta, PD − healthy (×1e-3)', 'norm_diff'),
                                         ('ν·m (Hz): PD − healthy', 'num_diff'),
                                         ('PPC magnitude, PD − healthy', 'ppc_diff')]):
            a = ax[i, j]
            for d, ls in (('graded', '-'), ('uniform', '--')):
                p = profs.get(f'{v}|{d}')
                if p is None:
                    continue
                if fn is None:
                    a.plot(POS, p['rate_PVP'], ls, color='tab:orange', label=f'PV+ PD {d}')
                    a.plot(POS, p['rate_PVN'], ls, color='tab:green', label=f'PV− PD {d}')
                    if d == 'graded':
                        a.plot(POS, p['rate_PVP_h'], ':', color='tab:orange'); a.plot(POS, p['rate_PVN_h'], ':', color='tab:green')
                        lo = p['rate_neff_PVP'] < FLAG
                        a.plot(POS[lo], p['rate_PVP'][lo], 'o', color='tab:red', ms=2)
                    continue
                y = np.array(p[fn]) * (1e3 if fn == 'norm_diff' else 1)
                ne = np.array(p['ppc_neff'] if fn == 'ppc_diff' else p['neff'])
                a.plot(POS, y, ls, color='k' if d == 'graded' else '0.5', lw=1.6, label=f'{d} depth')
                a.plot(POS[ne < FLAG], y[ne < FLAG], 'o', color='tab:red', ms=3)
            a.axvspan(0.55, 0.8, color='0.93', zorder=-1)
            if f'{v}|graded' not in profs and i == 0:
                a.text(0.5, 0.5, f'{v} held: not run\n(seed-1 smoke: PV+ rate +20% vs T0)', ha='center', va='center', transform=a.transAxes, fontsize=9)
            if i == 0:
                a.set_title(v, fontsize=10)
            if j == 0:
                a.set_ylabel(title, fontsize=8)
            if j == 0 and i < 2:
                a.legend(fontsize=6)
    for a in ax[-1]:
        a.set_xlabel('STN axis (0 ventromedial, 1 dorsolateral)')
    fig.suptitle('Cortical topography T0–T3 × pallidal depth (F1 + F2; red dots: n_eff < 8)', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'maps.png'), dpi=110)


if __name__ == '__main__':
    main()
