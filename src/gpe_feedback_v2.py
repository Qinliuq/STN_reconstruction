"""Low-dimensional GPe input generator with delayed STN-to-GPe feedback (v2: with a pallidal beta drive).

v2 changes from gpe_feedback.py: in addition to the STN feedback, the GPe output can carry a
beta-band (default 20 Hz) modulation during beta bursts, shared by all GPe groups (common phase,
small per-group jitter). It models parkinsonian pallidal beta and makes GPe the primary beta source
for the STN; the modulation leaves the mean GPe rate unchanged. STATES holds healthy and
parkinsonian settings (pass them as keyword arguments).

Each GPe source group has one rate variable. Its rate relaxes toward a target
that increases with delayed mean STN firing in that group's assigned STN
subpopulation. GPe events are then generated independently for connected
STN cells and returned as inhibitory synaptic arrivals.

This is deliberately reduced: it omits GPe neuron dynamics, striatal input,
GPe-GPe connections, and detailed axonal/synaptic structure.

Run a standalone demonstration (NumPy required):
    python src/gpe_feedback_v2.py --neurons 2000 --duration 2 --out work/gpe_events.npz

Use in an STN simulator:
    gpe = GPeInputGenerator(n_stn=2000, dt_s=0.001, seed=4)
    for step in range(n_steps):
        stn_spiked_ids = ...  # IDs that fired in this STN step
        inhibitory_arrival_ids = gpe.step(stn_spiked_ids)
        # Apply one inhibitory conductance jump to each returned STN ID.
        # An ID can occur more than once if multiple arrivals occurred.

Returned arrival IDs are for the current step. STN spikes supplied to step()
modify GPe firing after the configured feedback delay.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

STATES = {
    # healthy: weak, brief beta episodes; parkinsonian: long, deep beta bursts (timing as in the cortical model).
    # Every GPe group projects to every STN cell (target_probability 1.0): with the default 0.35, each STN cell
    # received too little GPe input for the pallidal beta to entrain it (calibrated in sweeps/gpe_beta_tuning*.py).
    "healthy": dict(beta_modulation=0.10, burst_rate_hz=0.5, burst_duration_s=0.15, target_probability=1.0),
    "parkinsonian": dict(beta_modulation=1.00, burst_rate_hz=1.0, burst_duration_s=0.50, target_probability=1.0),
}


class GPeInputGenerator:
    """A few rate-based GPe groups that emit inhibitory arrivals to STN cells."""

    def __init__(
        self,
        n_stn: int,
        dt_s: float,
        n_gpe_groups: int = 12,
        target_probability: float = 0.35,
        gpe_baseline_hz: float = 25.0,
        stn_reference_hz: float = 10.0,
        stn_to_gpe_gain: float = 0.5,
        gpe_rate_tau_s: float = 0.020,
        feedback_delay_s: float = 0.005,
        gpe_rate_max_hz: float = 100.0,
        beta_hz: float = 20.0,
        beta_modulation: float = 0.0,
        burst_rate_hz: float = 0.0,
        burst_duration_s: float = 0.0,
        beta_phase_jitter: float = 0.3,
        seed: int = 1,
    ) -> None:
        if n_stn <= 0 or dt_s <= 0 or n_gpe_groups <= 0:
            raise ValueError("n_stn, dt_s, and n_gpe_groups must be positive")
        if not 0.0 <= target_probability <= 1.0:
            raise ValueError("target_probability must be in [0, 1]")
        if gpe_rate_tau_s <= 0 or feedback_delay_s < 0:
            raise ValueError("rate time constant must be positive; delay nonnegative")

        self.n_stn = n_stn
        self.dt_s = dt_s
        self.n_groups = min(n_gpe_groups, n_stn)
        self.rng = np.random.default_rng(seed)
        self.gpe_baseline_hz = gpe_baseline_hz
        self.stn_reference_hz = stn_reference_hz
        self.stn_to_gpe_gain = stn_to_gpe_gain
        self.gpe_rate_tau_s = gpe_rate_tau_s
        self.delay_steps = int(round(feedback_delay_s / dt_s))
        self.gpe_rate_max_hz = gpe_rate_max_hz
        self.step_index = 0

        # Divide STN cells into balanced groups; each GPe group receives
        # feedback from one such subgroup. Randomize IDs to avoid index artifacts.
        shuffled = self.rng.permutation(n_stn)
        self.stn_to_gpe_group = np.empty(n_stn, dtype=np.int32)
        self.stn_group_sizes = np.zeros(self.n_groups, dtype=np.int32)
        for group_id, members in enumerate(np.array_split(shuffled, self.n_groups)):
            self.stn_to_gpe_group[members] = group_id
            self.stn_group_sizes[group_id] = len(members)

        # Each GPe group inhibits a random subset of the STN. This makes the
        # source groups distinct while keeping connectivity sparse.
        self.targets_by_gpe_group = []
        for _ in range(self.n_groups):
            mask = self.rng.random(n_stn) < target_probability
            targets = np.flatnonzero(mask).astype(np.int32)
            # Ensure each source has at least one target.
            if len(targets) == 0:
                targets = np.array([self.rng.integers(n_stn)], dtype=np.int32)
            self.targets_by_gpe_group.append(targets)

        # pallidal beta drive: emitted rate = rate_hz * (1 + beta_modulation * burst * cos(phase)), clipped at 0;
        # bursts start as a Poisson process and last burst_duration_s; all groups share the phase up to a jitter (rad)
        if not 0.0 <= beta_modulation <= 1.0:
            raise ValueError("beta_modulation must be in [0, 1]")
        self.beta_hz = beta_hz
        self.beta_modulation = beta_modulation
        self.burst_start_prob = burst_rate_hz * dt_s
        self.burst_steps = int(round(burst_duration_s / dt_s))
        self.burst_steps_left = 0
        self.beta_phase = self.rng.uniform(0.0, 2.0 * np.pi) + self.rng.normal(0.0, beta_phase_jitter, self.n_groups)
        self.in_burst = False

        self.rate_hz = np.full(self.n_groups, gpe_baseline_hz, dtype=np.float64)
        self.emitted_rate_hz = self.rate_hz.copy()
        self.history = np.full(
            (self.delay_steps + 1, self.n_groups), stn_reference_hz, dtype=np.float64
        )

    def step(self, stn_spiked_ids: np.ndarray | list[int]) -> np.ndarray:
        """Advance one time step and return STN IDs receiving GPe arrivals.

        Parameters
        ----------
        stn_spiked_ids:
            STN neuron IDs that fired during this time bin. An empty sequence
            means no STN spikes. IDs are treated as distinct spikes; repeated
            IDs represent multiple spikes from the same neuron in the bin.
        """
        ids = np.asarray(stn_spiked_ids, dtype=np.int64)
        if ids.size and (ids.min() < 0 or ids.max() >= self.n_stn):
            raise ValueError("STN spike ID outside [0, n_stn)")

        # Use the delayed rate history to update the GPe population rates.
        if self.step_index >= self.delay_steps:
            delayed_slot = (self.step_index - self.delay_steps) % len(self.history)
            delayed_stn_rate = self.history[delayed_slot]
        else:
            delayed_stn_rate = np.full(self.n_groups, self.stn_reference_hz)

        target_rate = self.gpe_baseline_hz + self.stn_to_gpe_gain * (
            delayed_stn_rate - self.stn_reference_hz
        )
        target_rate = np.clip(target_rate, 0.0, self.gpe_rate_max_hz)
        alpha = 1.0 - np.exp(-self.dt_s / self.gpe_rate_tau_s)
        self.rate_hz += alpha * (target_rate - self.rate_hz)

        # Beta burst envelope and modulation of the emitted rate.
        if self.burst_steps_left > 0:
            self.burst_steps_left -= 1
        elif self.burst_steps > 0 and self.rng.random() < self.burst_start_prob:
            self.burst_steps_left = self.burst_steps
        self.in_burst = self.burst_steps_left > 0
        modulation = 1.0
        if self.in_burst and self.beta_modulation > 0:
            t = self.step_index * self.dt_s
            modulation = np.maximum(0.0, 1.0 + self.beta_modulation * np.cos(2.0 * np.pi * self.beta_hz * t + self.beta_phase))
        self.emitted_rate_hz = self.rate_hz * modulation

        # Draw independent Poisson arrivals on each connected GPe->STN edge.
        arrival_chunks = []
        for group_id, targets in enumerate(self.targets_by_gpe_group):
            counts = self.rng.poisson(self.emitted_rate_hz[group_id] * self.dt_s, size=len(targets))
            if np.any(counts):
                arrival_chunks.append(np.repeat(targets, counts))
        arrivals = (
            np.concatenate(arrival_chunks).astype(np.int32, copy=False)
            if arrival_chunks
            else np.empty(0, dtype=np.int32)
        )

        # Save current-bin STN group rates for future delayed feedback.
        if ids.size:
            spike_counts = np.bincount(
                self.stn_to_gpe_group[ids], minlength=self.n_groups
            )
        else:
            spike_counts = np.zeros(self.n_groups, dtype=np.int64)
        stn_rates = spike_counts / (self.stn_group_sizes * self.dt_s)
        write_slot = self.step_index % len(self.history)
        self.history[write_slot] = stn_rates
        self.step_index += 1
        return arrivals


def run_demo(
    n_stn: int,
    duration_s: float,
    dt_s: float,
    n_gpe_groups: int,
    stn_demo_rate_hz: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run with synthetic STN spikes so the script works standalone."""
    rng = np.random.default_rng(seed)
    gpe = GPeInputGenerator(
        n_stn=n_stn,
        dt_s=dt_s,
        n_gpe_groups=n_gpe_groups,
        seed=seed + 1,
    )
    n_steps = int(np.ceil(duration_s / dt_s))
    time_chunks: list[np.ndarray] = []
    id_chunks: list[np.ndarray] = []
    rate_history = np.empty((n_steps, gpe.n_groups), dtype=np.float32)

    for step_idx in range(n_steps):
        # Demo only: independent baseline STN spiking. Replace this with the
        # actual neuron IDs that fired in your STN model at this time step.
        stn_spikes = np.flatnonzero(rng.random(n_stn) < stn_demo_rate_hz * dt_s)
        arrivals = gpe.step(stn_spikes)
        rate_history[step_idx] = gpe.rate_hz
        if arrivals.size:
            time_chunks.append(np.full(arrivals.size, step_idx * dt_s, dtype=np.float64))
            id_chunks.append(arrivals)

    if time_chunks:
        time_s = np.concatenate(time_chunks)
        neuron_id = np.concatenate(id_chunks)
        order = np.argsort(time_s, kind="stable")
        time_s, neuron_id = time_s[order], neuron_id[order]
    else:
        time_s = np.empty(0, dtype=np.float64)
        neuron_id = np.empty(0, dtype=np.int32)
    return time_s, neuron_id, rate_history, gpe.stn_to_gpe_group


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neurons", type=int, default=2000)
    parser.add_argument("--duration", type=float, default=2.0, help="seconds")
    parser.add_argument("--dt", type=float, default=0.001, help="seconds")
    parser.add_argument("--gpe-groups", type=int, default=12)
    parser.add_argument("--stn-demo-rate", type=float, default=10.0, help="Hz")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--conductance-ns", type=float, default=0.5)
    parser.add_argument("--out", type=Path, default=Path("work/gpe_events.npz"))
    args = parser.parse_args()

    time_s, neuron_id, rate_history, stn_group = run_demo(
        n_stn=args.neurons,
        duration_s=args.duration,
        dt_s=args.dt,
        n_gpe_groups=args.gpe_groups,
        stn_demo_rate_hz=args.stn_demo_rate,
        seed=args.seed,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        time_s=time_s,
        neuron_id=neuron_id,
        conductance_nS=np.full(time_s.shape, args.conductance_ns, dtype=np.float32),
        gpe_rate_hz=rate_history,
        stn_to_gpe_group=stn_group,
        dt_s=np.float64(args.dt),
    )
    print(f"Saved {len(time_s):,} inhibitory arrivals to {args.out}")
    print("The standalone demo uses synthetic 10 Hz STN spikes; replace them with your STN model's spike IDs.")


if __name__ == "__main__":
    main()
