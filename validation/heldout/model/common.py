"""Model-side helpers for the held-out predictions. Never imports the experimental data (tests/test_isolation.py)."""
import datetime
import glob
import hashlib
import json
import os
import runpy
import subprocess

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
PRIORS = json.load(open(os.path.join(ROOT, 'validation', 'heldout', 'priors.json'), encoding='utf-8'))
T0_S = PRIORS['transient_ms'] / 1000


def out_dir(test):
    d = os.path.join(ROOT, 'results', 'heldout', test)
    os.makedirs(d, exist_ok=True)
    return d


def git_info():
    run = lambda *a: subprocess.run(['git', *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return dict(hash=run('rev-parse', 'HEAD'), dirty=bool(run('status', '--porcelain', '--', 'src', 'validation')))


def load_runs(spec_name):
    """[(label, seed, npz dict, cfg text, path)] for all finished runs of a held-out spec."""
    spec = runpy.run_path(os.path.join(ROOT, 'validation', 'heldout', 'model', 'specs', f'{spec_name}.py'))
    runs = []
    for cond in spec['conditions']:
        for seed in spec['seeds']:
            d = os.path.join(ROOT, spec['outdir'], f"{cond['label']}_s{seed}")
            f = os.path.join(d, 'model_output_spikes.npz')
            if not os.path.exists(f):
                raise SystemExit(f'missing run {f}: run python -m validation.heldout.model.run_all first')
            runs.append((cond['label'], seed, dict(np.load(f)), open(os.path.join(d, 'cfg_run.py'), encoding='utf-8').read(), f))
    return runs


def cell_trains(d, t0_s=T0_S):
    """[(cell index, 'PV+'/'PV-', axis position, spike times in s (t >= t0_s))] in STN index order (PV+ first)."""
    out = []
    i = 0
    for pop, typ in (('PVP_pop', 'PV+'), ('PVN_pop', 'PV-')):
        for g in d[f'gids_{pop}']:
            st = np.sort(d['spkt'][d['spkid'] == g]) / 1000
            out.append((i, typ, float(d['axis'][i]), st[st >= t0_s]))
            i += 1
    return out


def sha256(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def write_prediction(test, metrics, runs, files=()):
    """results/heldout/<test>/prediction.json: model metrics + git hash + config, written before any data is read."""
    d = out_dir(test)
    pred = dict(test=test, created=datetime.datetime.now().isoformat(timespec='seconds'), git=git_info(),
                priors=PRIORS, runs=[dict(label=l, seed=s, spikes=os.path.relpath(p, ROOT), spikes_sha256=sha256(p), cfg=c)
                                     for l, s, _, c, p in runs],
                files={os.path.relpath(f, ROOT): sha256(f) for f in files}, metrics=metrics)
    with open(os.path.join(d, 'prediction.json'), 'w', encoding='utf-8') as fh:
        json.dump(pred, fh, indent=1, default=_json)
    print(f'wrote {os.path.join(d, "prediction.json")}')


def _json(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def mean_sd(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    return dict(mean=float(v.mean()) if len(v) else np.nan, sd=float(v.std(ddof=1)) if len(v) > 1 else np.nan, n=len(v),
                values=v.tolist())
