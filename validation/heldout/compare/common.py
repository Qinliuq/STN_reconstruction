"""Comparison-side helpers: the prediction gate, verdict writing and plotting defaults."""
import datetime
import hashlib
import json
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))


def out_dir(test):
    return os.path.join(ROOT, 'results', 'heldout', test)


def require_prediction(test):
    """Load results/heldout/<test>/prediction.json; refuse to run if it is missing or its files have changed since."""
    p = os.path.join(out_dir(test), 'prediction.json')
    if not os.path.exists(p):
        raise SystemExit(f'{p} is missing: run the model prediction (validation/heldout/model/predict_*.py) first. '
                         'The comparison data are not read before the prediction exists.')
    pred = json.load(open(p))
    for f, h in pred['files'].items():
        if hashlib.sha256(open(os.path.join(ROOT, f), 'rb').read()).hexdigest() != h:
            raise SystemExit(f'{f} changed after prediction.json was written; rerun the prediction')
    print(f"prediction for {test}: written {pred['created']}, git {pred['git']['hash'][:8]}{' (dirty)' if pred['git']['dirty'] else ''}")
    pred['_compared'] = datetime.datetime.now().isoformat(timespec='seconds')
    return pred


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    c = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    m = (p + z * z / (2 * n)) / (1 + z * z / n)
    return (m - c, m + c)


def write_summary(test, verdict_line, sections, pred, table=None):
    with open(os.path.join(out_dir(test), 'summary.md'), 'w') as fh:
        fh.write(f'# {test}\n\n**{verdict_line}**\n\n')
        fh.write(f"Prediction written {pred['created']} (git {pred['git']['hash'][:8]}"
                 f"{', uncommitted changes' if pred['git']['dirty'] else ''}); compared {pred['_compared']}.\n\n")
        for title, body in sections:
            fh.write(f'## {title}\n\n{body}\n\n')
    with open(os.path.join(out_dir(test), 'comparison.json'), 'w') as fh:
        json.dump(dict(verdict=verdict_line, **(table or {})), fh, indent=1)
    print(verdict_line)


def fmt(ms, nd=2):
    return f"{ms['mean']:.{nd}f} ± {ms['sd']:.{nd}f}" if np.isfinite(ms.get('sd', np.nan)) else f"{ms['mean']:.{nd}f}"


def md_table(df, floatfmt='.3g'):
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        lines.append('| ' + ' | '.join(f'{v:{floatfmt}}' if isinstance(v, (float, np.floating)) else str(v) for v in r.values) + ' |')
    return '\n'.join(lines)


def style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    return plt


COL = {'data': '#4a4a4a', 'model': '#d0632b', 'bursty': '#c0392b', 'regular': '#2e86c1', 'PV+': '#8e44ad', 'PV-': '#16a085'}
