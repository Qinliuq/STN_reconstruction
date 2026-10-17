"""Generate synthetic cortical synaptic arrival events for an STN network (v2).

v2 changes from cortical_arrivals.py (STATES only): GPe (gpe_feedback_v2.py) is the primary beta
source, so the parkinsonian cortical beta modulation is reduced from 0.70 to 0.20 (a supplement);
the shared (assembly) fraction and burst timing are unchanged. The conductance per arrival is
1.0 nS in both states (calibrated in the network; the v1 default was 0.2 nS).

The output is an event list: (time_s, neuron_id, conductance_nS). It is a
reduced Poisson input model, not a reconstruction of individual cortical cells.

Example:
    python src/cortical_arrivals_v2.py --state parkinsonian --neurons 2000 --duration 2 --out work/cortical_events.npz

Requires NumPy. Load output with:
    data = numpy.load("work/cortical_events.npz")
    times_s = data["time_s"]
    neuron_ids = data["neuron_id"]
    conductance_nS = data["conductance_nS"]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class InputParams:
    # Mean cortical event arrivals per STN neuron, before beta modulation.
    rate_hz: float
    # Fraction of a neuron's arrivals supplied by a shared assembly process.
    shared_fraction: float
    # Sinusoidal rate modulation depth, 0..1. Applied during beta bursts.
    beta_modulation: float
    beta_hz: float
    # Burst envelope: average burst starts per second and fixed duration.
    burst_rate_hz: float
    burst_duration_s: float
    # Effective excitatory conductance jump per arrival (nS).
    conductance_nS: float


STATES = {
    # Illustrative defaults only; calibrate to the firing rates and synaptic
    # conductances of your STN model and experimental condition.
    "healthy": InputParams(
        rate_hz=80.0,
        shared_fraction=0.05,
        beta_modulation=0.10,
        beta_hz=20.0,
        burst_rate_hz=0.5,
        burst_duration_s=0.15,
        conductance_nS=1.0,
    ),
    "parkinsonian": InputParams(
        rate_hz=80.0,
        shared_fraction=0.25,
        beta_modulation=0.20,  # v1: 0.70; cortical beta is a supplement to the GPe beta
        beta_hz=20.0,
        burst_rate_hz=1.0,
        burst_duration_s=0.50,
        conductance_nS=1.0,
    ),
}


def make_beta_envelope(
    n_steps: int,
    dt_s: float,
    burst_rate_hz: float,
    burst_duration_s: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return a 0/1 burst envelope with Poisson-distributed burst onsets."""
    envelope = np.zeros(n_steps, dtype=np.float32)
    if burst_rate_hz <= 0 or burst_duration_s <= 0:
        return envelope

    duration_s = n_steps * dt_s
    n_bursts = rng.poisson(burst_rate_hz * duration_s)
    burst_steps = max(1, int(round(burst_duration_s / dt_s)))
    starts = rng.integers(0, n_steps, size=n_bursts)
    for start in starts:
        end = min(n_steps, start + burst_steps)
        envelope[start:end] = 1.0
    return envelope


def generate_arrivals(
    n_neurons: int,
    duration_s: float,
    dt_s: float,
    params: InputParams,
    n_assemblies: int = 64,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate cortical arrivals using a vectorized discrete-time Poisson draw.

    Neurons are divided among cortical assemblies. Within each time bin, an
    assembly gets one shared Poisson event count, delivered to every member of
    that assembly. Each STN neuron also gets independent private events. Thus
    shared_fraction controls the fraction of its mean arrivals that are common
    to the other neurons in its assembly.

    The returned time is the bin's left edge. Choose dt_s small relative to the
    STN integration step / synaptic time constants. For very long simulations,
    generate and save in chunks to limit memory use.
    """
    if n_neurons <= 0 or duration_s <= 0 or dt_s <= 0:
        raise ValueError("neurons, duration, and dt must be positive")
    if n_assemblies <= 0:
        raise ValueError("n_assemblies must be positive")
    if params.rate_hz < 0 or not 0 <= params.shared_fraction <= 1:
        raise ValueError("rate must be nonnegative and shared_fraction in [0, 1]")
    if not 0 <= params.beta_modulation <= 1:
        raise ValueError("beta_modulation must be in [0, 1]")

    rng = np.random.default_rng(seed)
    n_steps = int(np.ceil(duration_s / dt_s))
    actual_duration_s = n_steps * dt_s
    n_groups = min(n_assemblies, n_neurons)
    # Balanced contiguous membership. Shuffle so assembly membership is not
    # tied to neuron ID ordering.
    shuffled_ids = rng.permutation(n_neurons)
    groups = np.array_split(shuffled_ids, n_groups)

    burst = make_beta_envelope(
        n_steps, dt_s, params.burst_rate_hz, params.burst_duration_s, rng
    )
    phase = rng.uniform(0.0, 2.0 * np.pi)
    t = np.arange(n_steps, dtype=np.float64) * dt_s
    beta = np.cos(2.0 * np.pi * params.beta_hz * t + phase).astype(np.float32)
    # Keep instantaneous rate nonnegative. Beta modulation is active only
    # inside bursts; outside bursts, the rate is the unmodulated baseline.
    modulation = np.maximum(0.0, 1.0 + params.beta_modulation * burst * beta)

    shared_rate = params.rate_hz * params.shared_fraction
    private_rate = params.rate_hz * (1.0 - params.shared_fraction)

    time_chunks: list[np.ndarray] = []
    id_chunks: list[np.ndarray] = []

    # Generate private arrivals independently for each STN neuron.
    private_means = (private_rate * modulation * dt_s).astype(np.float64)
    for start in range(0, n_steps, 256):
        stop = min(n_steps, start + 256)
        counts = rng.poisson(private_means[start:stop, None], size=(stop - start, n_neurons))
        step_idx, neuron_idx = np.nonzero(counts)
        if len(step_idx):
            multiplicity = counts[step_idx, neuron_idx]
            time_chunks.append((start + step_idx).repeat(multiplicity).astype(np.int64))
            id_chunks.append(neuron_idx.repeat(multiplicity).astype(np.int32))

    # Each assembly has a shared event process. A shared event reaches every
    # STN cell in that assembly, modeling common cortical drive.
    shared_means = (shared_rate * modulation * dt_s).astype(np.float64)
    for member_ids in groups:
        for start in range(0, n_steps, 256):
            stop = min(n_steps, start + 256)
            counts = rng.poisson(shared_means[start:stop])
            steps = np.flatnonzero(counts)
            if len(steps):
                multiplicity = counts[steps]
                event_steps = (start + steps).repeat(multiplicity).astype(np.int64)
                event_ids = np.tile(member_ids.astype(np.int32), len(event_steps))
                event_steps = np.repeat(event_steps, len(member_ids))
                time_chunks.append(event_steps)
                id_chunks.append(event_ids)

    if time_chunks:
        step_idx = np.concatenate(time_chunks)
        neuron_id = np.concatenate(id_chunks)
        order = np.argsort(step_idx, kind="stable")
        step_idx = step_idx[order]
        neuron_id = neuron_id[order]
        time_s = step_idx.astype(np.float64) * dt_s
        conductance = np.full(time_s.shape, params.conductance_nS, dtype=np.float32)
    else:
        time_s = np.empty(0, dtype=np.float64)
        neuron_id = np.empty(0, dtype=np.int32)
        conductance = np.empty(0, dtype=np.float32)

    return time_s, neuron_id, conductance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", choices=STATES, default="parkinsonian")
    parser.add_argument("--neurons", type=int, default=2000)
    parser.add_argument("--duration", type=float, default=2.0, help="seconds")
    parser.add_argument("--dt", type=float, default=0.001, help="seconds")
    parser.add_argument("--assemblies", type=int, default=64)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", type=Path, default=Path("work/cortical_events.npz"))
    args = parser.parse_args()

    params = STATES[args.state]
    time_s, neuron_id, conductance = generate_arrivals(
        n_neurons=args.neurons,
        duration_s=args.duration,
        dt_s=args.dt,
        params=params,
        n_assemblies=args.assemblies,
        seed=args.seed,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        time_s=time_s,
        neuron_id=neuron_id,
        conductance_nS=conductance,
        n_neurons=np.int32(args.neurons),
        duration_s=np.float64(args.duration),
        dt_s=np.float64(args.dt),
        state=np.array(args.state),
    )
    print(f"Saved {len(time_s):,} events for {args.neurons:,} STN neurons to {args.out}")
    print(f"Mean arrivals per neuron: {len(time_s) / args.neurons / args.duration:.2f} / s")
    print("To create a rate-matched control, use --state healthy with the same seed and sizes.")


if __name__ == "__main__":
    main()
