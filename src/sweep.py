"""
Run a parameter sweep of the network (src/init.py) in parallel and summarize it.

A sweep spec is a Python file (see sweeps/) defining:
    duration = 3000        # ms
    seeds = [1, 2, 3]      # each run sets cfg.ctx_seed = seed and cfg.gpe_seed = 100 + seed (a condition may set its own 'seeds')
    conditions = [         # 'cfg' entries override src/cfg.py; 'group' and 'x' are used for plotting
        {'label': 'healthy', 'group': 'healthy', 'x': 1.0, 'cfg': {'ctx_state': 'healthy'}},
        {'label': 'pd_kir0.5', 'group': 'parkinsonian', 'x': 0.5, 'cfg': {'ctx_state': 'parkinsonian', 'channel_block': {'gKir': 0.5}}},
    ]
    x_label = 'Kir factor'             # optional, x axis of the summary plot
    reference_group = 'healthy'        # optional: plotted as a band; used for "gap closed" (needs baseline_label)
    band_groups = ['healthy', 'untreated']  # optional: groups plotted as horizontal bands (default: [reference_group])
    baseline_label = 'pd_kir1'         # optional: untreated condition for "gap closed"

Usage (from the repo root, after `nrnivmodl mod`):
    python src/sweep.py run sweeps/kir_block.py --jobs 4
    python src/sweep.py summary sweeps/kir_block.py
Each run goes to sweeps/out/<spec name>/<label>_s<seed>/ (config, log, metrics); runs that already have
metrics are skipped unless --force. LFP and NetPyNE plots are turned off in sweep runs for speed.
"""
import argparse
import csv
import json
import os
import runpy
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

METRICS = [('beta_pow', 'STN beta power (13-30 Hz)'), ('sync_corr10ms', 'synchrony (corr, 10 ms)'),
           ('rate_PVP_pop', 'PV+ rate (Hz)'), ('rate_PVN_pop', 'PV- rate (Hz)'),
           ('frac_spikes_in_bursts', 'spikes in bursts'), ('low_pow_3_12', '3-12 Hz power'),
           ('peak_hz', 'PSD peak (Hz)'), ('gpe_rate', 'GPe rate (Hz)')]


def load_spec(path):
    spec = runpy.run_path(path)
    spec['name'] = os.path.splitext(os.path.basename(path))[0]
    spec['outdir'] = os.path.join('sweeps', 'out', spec['name'])
    return spec


def run_dir(spec, label, seed):
    return os.path.join(spec['outdir'], f'{label}_s{seed}')


def write_cfg(spec, cond, seed, d):
    lines = ["exec(open('src/cfg.py').read())", f"cfg.duration = {spec['duration']}"]
    lines += [f"cfg.{k} = {v!r}" for k, v in cond['cfg'].items()]
    lines += [f"cfg.ctx_seed = {seed}", f"cfg.gpe_seed = {100 + seed}", f"cfg.filename = {os.path.join(d, 'model_output')!r}",
              "cfg.recordLFP = None", "cfg.analysis = {}", "cfg.savePickle = False"]
    path = os.path.join(d, 'cfg_run.py')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return path


def run(spec, jobs, force):
    todo = []
    for cond in spec['conditions']:
        for seed in cond.get('seeds', spec['seeds']):
            d = run_dir(spec, cond['label'], seed)
            if not force and os.path.exists(os.path.join(d, 'model_output_metrics.json')):
                continue
            os.makedirs(d, exist_ok=True)
            todo.append((d, write_cfg(spec, cond, seed, d)))
    print(f"{len(todo)} runs to do in {spec['outdir']} ({jobs} in parallel)")

    def one(item):
        d, cfgPath = item
        with open(os.path.join(d, 'run.log'), 'w') as log:
            ret = subprocess.run([sys.executable, 'src/init.py', f'simConfig={cfgPath}', 'netParams=src/netParams.py'],
                                 stdout=log, stderr=subprocess.STDOUT, env={**os.environ, 'MPLBACKEND': 'Agg'})
        ok = os.path.exists(os.path.join(d, 'model_output_metrics.json'))
        print(f"  {os.path.basename(d)}: {'done' if ok else f'FAILED (exit {ret.returncode}, see run.log)'}", flush=True)
        return ok

    with ThreadPoolExecutor(jobs) as ex:
        results = list(ex.map(one, todo))
    print(f"{sum(results)}/{len(results)} runs succeeded")


def summary(spec):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    rows = []
    for cond in spec['conditions']:
        for seed in cond.get('seeds', spec['seeds']):
            d = run_dir(spec, cond['label'], seed)
            mf = os.path.join(d, 'model_output_metrics.json')
            if os.path.exists(mf):
                psd = np.load(os.path.join(d, 'model_output_psd.npz'))
                rows.append(dict(label=cond['label'], group=cond.get('group', cond['label']), x=cond.get('x', 0.0),
                                 seed=seed, f=psd['f'], P=psd['P'], **json.load(open(mf))))
    if not rows:
        print('no results yet')
        return
    keys = [(k, lab) for k, lab in METRICS if k in rows[0]]
    labels = [c['label'] for c in spec['conditions'] if any(r['label'] == c['label'] for r in rows)]
    byLabel = {lab: [r for r in rows if r['label'] == lab] for lab in labels}
    mean = lambda rs, k: float(np.mean([r[k] for r in rs]))
    sd = lambda rs, k: float(np.std([r[k] for r in rs]))

    print(f"{'condition':18s} {'n':>2s} " + ' '.join(f'{lab[:15]:>17s}' for _, lab in keys))
    for lab in labels:
        rs = byLabel[lab]
        print(f"{lab:18s} {len(rs):2d} " + ' '.join(f"{mean(rs, k):9.4g}±{sd(rs, k):<7.2g}" for k, _ in keys))
    with open(os.path.join(spec['outdir'], 'summary.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['condition', 'n'] + [f'{k}_{s}' for k, _ in keys for s in ('mean', 'sd')])
        for lab in labels:
            w.writerow([lab, len(byLabel[lab])] + [f'{v:.6g}' for k, _ in keys for v in (mean(byLabel[lab], k), sd(byLabel[lab], k))])

    ref, base = spec.get('reference_group'), spec.get('baseline_label')
    refRows = [r for r in rows if r['group'] == ref] if ref else []
    if refRows and base in byLabel:
        print(f"\nfraction of the gap between '{base}' and '{ref}' closed (1 = reaches {ref}, < 0 = moves away)")
        for k, lab in keys[:2]:
            gap = mean(byLabel[base], k) - mean(refRows, k)
            print(f"  {lab:28s}" + ''.join(f"  {l}: {(mean(byLabel[base], k) - mean(byLabel[l], k)) / gap:5.2f}"
                                          for l in labels if l != base and byLabel[l][0]['group'] != ref))

    groups = []
    for r in rows:
        if r['group'] not in groups:
            groups.append(r['group'])
    colors = dict(zip(groups, plt.cm.tab10.colors))
    fig, axs = plt.subplots(2, 4, figsize=(15, 6.5))
    for ax, (k, lab) in zip(axs.flat, keys):
        for g in groups:
            gl = [l for l in labels if byLabel[l][0]['group'] == g]
            xs = [byLabel[l][0]['x'] for l in gl]
            if g in spec.get('band_groups', [ref]):
                mu, s_ = mean(byLabel[gl[0]], k), sd(byLabel[gl[0]], k)
                ax.axhline(mu, color=colors[g], label=g)
                ax.axhspan(mu - s_, mu + s_, color=colors[g], alpha=0.15)
            else:
                order = np.argsort(xs)
                for l in gl:
                    ax.plot([byLabel[l][0]['x']] * len(byLabel[l]), [r[k] for r in byLabel[l]], 'o', color=colors[g], alpha=0.3, ms=4)
                ax.plot(np.array(xs)[order], [mean(byLabel[gl[i]], k) for i in order], '-o', color=colors[g], label=g)
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel(spec.get('x_label', 'x'))
    axs[0, 0].legend(fontsize=8)
    fig.suptitle(f"{spec['name']} ({len(spec['seeds'])} seeds, {spec['duration']} ms)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(spec['outdir'], 'summary.png'), dpi=110)

    fig, ax = plt.subplots(figsize=(7, 4))
    for lab in labels:
        rs = byLabel[lab]
        ax.semilogy(rs[0]['f'], np.mean([r['P'] for r in rs], 0), label=lab, lw=2 if rs[0]['group'] == ref else 1)
    ax.axvspan(13, 30, color='grey', alpha=0.15)
    ax.set_xlim(0, 60)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD of STN population rate')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(spec['outdir'], 'psd.png'), dpi=110)
    if 'profile' in rows[0]:
        pkeys = [('beta_pow', 'beta power (13-30 Hz)'), ('beta_frac', 'beta fraction'), ('sync', 'synchrony (corr, 10 ms)'),
                 ('bursts', 'spikes in bursts'), ('rate', 'rate (Hz)'), ('n_eff', 'cells per position (n_eff)')]
        fig, axs = plt.subplots(2, 3, figsize=(14, 7))
        print(f"\n{'condition':18s} beta peak (axis)   synchrony peak (axis)")
        for lab in labels:
            rs = byLabel[lab]
            pos = rs[0]['profile']['positions']
            for ax, (k, title) in zip(axs.flat, pkeys):
                ax.plot(pos, np.mean([r['profile'][k] for r in rs], 0), label=lab, lw=2 if rs[0]['group'] == ref else 1)
                ax.set_title(title, fontsize=10)
                ax.set_xlabel('position along STN axis')
            peaks = lambda k: ', '.join(f'{r[k]:.2f}' for r in rs)
            print(f"{lab:18s} {mean(rs, 'beta_peak_axis'):5.2f} ({peaks('beta_peak_axis')})"
                  f"   {mean(rs, 'sync_peak_axis'):5.2f} ({peaks('sync_peak_axis')})")
        axs[0, 0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(spec['outdir'], 'profile.png'), dpi=110)
    print(f"\nsaved {spec['outdir']}/summary.png, psd.png, summary.csv (and profile.png if available)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run or summarize a network parameter sweep')
    parser.add_argument('action', choices=['run', 'summary'])
    parser.add_argument('spec', help='sweep spec file, e.g. sweeps/kir_block.py')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--force', action='store_true', help='rerun runs that already have metrics')
    args = parser.parse_args()
    spec = load_spec(args.spec)
    if args.action == 'run':
        run(spec, args.jobs, args.force)
    summary(spec)
