"""Build results/heldout/README.md from the per-test prediction.json / comparison outputs.

    python -m validation.heldout.compare.make_readme
"""
import json
import os

import pandas as pd

from validation.heldout.compare import common as C

ROWS = [
    ('H1', 'Rat 6-OHDA single units (unit_cluster_SSonly.csv, 73 SS rows / 45 distinct units)',
     'bimodality (silhouette vs Gaussian null, BIC), bursty fraction, cluster vs PV+'),
    ('H2', 'Rat 6-OHDA single units (73 SS rows)', 'pooled rate / CV / burst share (pooled check)'),
    ('H3', 'Rat 6-OHDA units + IpsiEEG (37 units with EEG)', 'PPC + preferred phase, slow 0.5-2 Hz and beta 13-30 Hz'),
    ('H4', 'Mouse D1/D2-MSN ChR2 (binned_response.csv; D2 wave_clus)', 'responder sign / fraction / magnitude, time course, PV enrichment'),
]
CAVEATS = {
    'H1': 'urethane anesthesia; 28/73 rat rows are copies of other units (n_eff 45); model seeds share one 58-cell network; 20 s model trains vs ~100 s rat trains',
    'H2': 'pooled check only, never tuned; anesthesia; species/state of rates differ',
    'H3': 'model slow oscillation is an imposed input envelope (priors.json); EEG polarity vs cortical firing unknown; no rat activated epochs inside unit epochs (beta test on activated cortex not possible)',
    'H4': 'mouse, dopamine-intact, anesthetized, multi-unit (~1-4 Hz); compressed 1/2/2 s protocol (slow recovery / after-effects untested); route magnitudes a priori',
}


def main():
    lines = ['# Held-out validation of configuration D', '',
             'Experimental targets that were never used for calibration. Each model prediction (`H<n>/prediction.json`, with git hash, '
             'configuration and spike-file hashes) was written before the comparison script read the data; the comparison scripts '
             'refuse to run without it. Input-regime parameters were fixed in `validation/heldout/priors.json` before any run. '
             'Configuration-D parameters were not changed. Code: `validation/heldout/` (see its README).', '',
             '| test | data source | model metric | model | data | statistic | verdict | caveats |', '|---|---|---|---|---|---|---|---|']
    for test, src, metric in ROWS:
        f = os.path.join(C.out_dir(test), 'comparison.json')
        if not os.path.exists(f):
            lines.append(f'| {test} | {src} | {metric} | pending | | | pending | {CAVEATS[test]} |')
            continue
        v = json.load(open(f))
        lines.append(f"| {test} | {src} | {metric} | {v.get('model', '')} | {v.get('data', '')} | {v.get('statistic', '')} | "
                     f"**{v.get('verdict_word', '')}** | {CAVEATS[test]} |")
    lines += ['', 'One-line verdicts:', '']
    for test, _, _ in ROWS:
        f = os.path.join(C.out_dir(test), 'comparison.json')
        if os.path.exists(f):
            lines.append(f"- {json.load(open(f))['verdict']}")
    lines += ['', 'Details: `H<n>/summary.md`, figures `H<n>/H<n>_summary.png`, tables `H<n>/*.csv`.']
    with open(os.path.join(C.ROOT, 'results', 'heldout', 'README.md'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
