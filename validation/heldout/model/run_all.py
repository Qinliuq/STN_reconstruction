"""
Run all held-out model simulations in one pool (order: H1, H3, H4), reusing src/sweep.py (same cfg writing, seed
convention and output layout). Runs with metrics already present are skipped.

    python -m validation.heldout.model.run_all --jobs 4
"""
import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, 'src')
import sweep  # noqa: E402

SPECS = ['validation/heldout/model/specs/H1_pd.py', 'validation/heldout/model/specs/H3_slowosc.py',
         'validation/heldout/model/specs/H4_opto.py']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--only', nargs='*', help='spec names to run, e.g. H1_pd')
    args = ap.parse_args()
    todo = []
    for path in SPECS:
        spec = sweep.load_spec(path)
        if args.only and spec['name'] not in args.only:
            continue
        for cond in spec['conditions']:
            for seed in cond.get('seeds', spec['seeds']):
                d = sweep.run_dir(spec, cond['label'], seed)
                if os.path.exists(os.path.join(d, 'model_output_spikes.npz')):
                    continue
                os.makedirs(d, exist_ok=True)
                todo.append((d, sweep.write_cfg(spec, cond, seed, d)))
    print(f'{len(todo)} runs ({args.jobs} in parallel)', flush=True)

    def one(item):
        d, cfgPath = item
        with open(os.path.join(d, 'run.log'), 'w') as log:
            subprocess.run([sys.executable, 'src/init.py', f'simConfig={cfgPath}', 'netParams=src/netParams.py'],
                           stdout=log, stderr=subprocess.STDOUT, env={**os.environ, 'MPLBACKEND': 'Agg'})
        ok = os.path.exists(os.path.join(d, 'model_output_spikes.npz'))
        print(f"  {d}: {'done' if ok else 'FAILED (see run.log)'}", flush=True)

    with ThreadPoolExecutor(args.jobs) as ex:
        list(ex.map(one, todo))


if __name__ == '__main__':
    main()
