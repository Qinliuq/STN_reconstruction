"""
Time-varying input-rate modulation (held-out validation hooks; off by default).

cfg.ctx_rate_mod / cfg.gpe_rate_mod are None by default, which leaves the model unchanged. Otherwise a dict:

- {'type': 'slow_osc', 'freq_hz': 1.0, 'depth': 0.9, 'duty': 0.5, 'edge_s': 0.05, 'period_cv': 0.1, 'seed': 0}
  anesthesia-like up/down states (cortex only): a train of cycles with period 1/freq_hz (jittered with CV
  period_cv), an up state for the first `duty` of each cycle and a down state at (1 - depth) x the up-state
  level, with raised-cosine edges of edge_s. The envelope is normalized to mean 1, so the mean cortical
  drive stays at its calibrated value. Phase (0..2pi per cycle, 0 = up-state onset) is returned too.
- {'type': 'blocks', 'pre_s': 1.0, 'on_s': 2.0, 'post_s': 2.0, 'factor': 0.5, 't_start_s': 0.5}
  a repeating stimulation protocol (pre / on / post) from t_start_s: the rate is multiplied by `factor` during
  the on windows and by 1 elsewhere. For the GPe, `factor` may be a list (one factor per GPe group).
"""
import numpy as np


def slow_osc_envelope(t_s, freq_hz=1.0, depth=0.9, duty=0.5, edge_s=0.05, period_cv=0.1, seed=0):
    t_s = np.asarray(t_s, dtype=float)
    rng = np.random.default_rng(seed)
    T = t_s[-1] + 1.0 / freq_hz if len(t_s) else 0.0
    periods = []
    while sum(periods) < T + 2.0 / freq_hz:
        periods.append(max(0.3 / freq_hz, rng.normal(1.0 / freq_hz, period_cv / freq_hz)))
    starts = np.concatenate([[0.0], np.cumsum(periods)]) - rng.uniform(0, periods[0]) # random start phase
    k = np.searchsorted(starts, t_s, side='right') - 1
    P = starts[k + 1] - starts[k]
    x = (t_s - starts[k]) # time into the cycle
    phase = 2 * np.pi * x / P
    up = duty * P
    # up-state indicator with raised-cosine edges
    def ramp(u):
        return 0.5 - 0.5 * np.cos(np.pi * np.clip(u, 0, 1))
    e = edge_s if edge_s > 0 else 1e-9
    level = np.where(x < up, ramp(x / e), 1 - ramp((x - up) / e))
    env = (1 - depth) + depth * level
    return env / env.mean(), phase


def block_factor(t_s, pre_s, on_s, post_s, factor, t_start_s=0.0):
    """Multiplicative factor at times t_s (scalar factor -> array over t; list factor -> t x groups)."""
    t_s = np.atleast_1d(np.asarray(t_s, dtype=float))
    cyc = pre_s + on_s + post_s
    x = np.mod(t_s - t_start_s, cyc)
    on = (t_s >= t_start_s) & (x >= pre_s) & (x < pre_s + on_s)
    f = np.asarray(factor, dtype=float)
    if f.ndim == 0:
        return np.where(on, float(f), 1.0)
    return np.where(on[:, None], f[None, :], 1.0)


def stim_onsets(duration_s, pre_s, on_s, post_s, t_start_s=0.0, **_):
    cyc = pre_s + on_s + post_s
    on = t_start_s + pre_s + np.arange(0, 10000) * cyc
    # only complete cycles
    return on[on + on_s + post_s <= duration_s + 1e-9]


def ctx_envelope(spec, n_steps, dt_s):
    """Envelope (n_steps,) for the cortical generator and a dict of signals to save."""
    t = np.arange(n_steps) * dt_s
    kw = {k: v for k, v in spec.items() if k != 'type'}
    if spec['type'] == 'slow_osc':
        env, phase = slow_osc_envelope(t, **kw)
        return env, {'ctx_env': env.astype(np.float32), 'ctx_so_phase': phase.astype(np.float32)}
    if spec['type'] == 'blocks':
        env = block_factor(t, **kw)
        return env, {'ctx_env': env.astype(np.float32)}
    raise ValueError(f"unknown ctx_rate_mod type {spec['type']}")


def gpe_factor(spec, t_s, n_groups):
    if spec['type'] != 'blocks':
        raise ValueError(f"gpe_rate_mod supports only 'blocks', got {spec['type']}")
    kw = {k: v for k, v in spec.items() if k != 'type'}
    f = block_factor([t_s], **kw)[0]
    return np.broadcast_to(f, (n_groups,)) if np.ndim(f) == 0 else f
