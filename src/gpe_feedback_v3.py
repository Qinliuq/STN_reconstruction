"""Low-dimensional GPe input generator with delayed STN-to-GPe feedback (v3: topographic).

v3 changes from gpe_feedback_v2.py:
- Topography: each GPe group has a home position along the STN axis (evenly spaced, 0..1). A GPe
  group inhibits STN cell i with probability target_probability * exp(-(axis_i - c_g)^2 / (2 topo_sigma^2)),
  so each STN cell receives from a few nearby groups (sparse and local, not all-to-all), and it
  receives STN feedback weighted by the same Gaussian (reciprocal local STN-GPe loops).
- The parkinsonian beta drive is graded along the axis: group g's modulation depth is
  beta_floor + (beta_modulation - beta_floor) * exp(-(c_g - beta_center)^2 / (2 beta_width^2)),
  i.e. strongest in groups projecting to the upper-middle STN.
- STATES also carry synapse_scale: the GPe->STN conductance multiplier of the state (applied by
  stn_inputs.GPeFeedback, not by this generator); parkinsonian > 1 models the strengthened GPe-STN
  synapses after dopamine depletion (Fan et al. 2012).

The STN cell positions along the axis (0..1) are passed as stn_axis.
"""

from __future__ import annotations

import numpy as np

STATES = {
    # healthy: weak, brief, spatially uniform beta episodes
    "healthy": dict(beta_modulation=0.10, beta_floor=0.10, burst_rate_hz=0.5, burst_duration_s=0.15,
                    synapse_scale=1.0),
    # parkinsonian: frequent, long beta bursts, deepest at the upper-middle axis (~0.65); stronger GPe->STN synapses
    "parkinsonian": dict(beta_modulation=1.00, beta_floor=0.20, beta_center=0.65, beta_width=0.15,
                         burst_rate_hz=2.0, burst_duration_s=0.40, synapse_scale=1.5),
}


class GPeInputGenerator:
    """Rate-based GPe groups with topographic connections to and from the STN."""

    def __init__(
        self,
        n_stn: int,
        dt_s: float,
        stn_axis: np.ndarray | list[float],
        n_gpe_groups: int = 12,
        target_probability: float = 1.0,
        topo_sigma: float = 0.10,
        gpe_baseline_hz: float = 25.0,
        stn_reference_hz: float = 10.0,
        stn_to_gpe_gain: float = 0.5,
        gpe_rate_tau_s: float = 0.020,
        feedback_delay_s: float = 0.005,
        gpe_rate_max_hz: float = 100.0,
        beta_hz: float = 20.0,
        beta_modulation: float = 0.0,
        beta_floor: float = 0.0,
        beta_center: float = 0.5,
        beta_width: float = 1e9,
        burst_rate_hz: float = 0.0,
        burst_duration_s: float = 0.0,
        beta_phase_jitter: float = 0.3,
        synapse_scale: float = 1.0,  # used by stn_inputs.GPeFeedback; accepted here so STATES can be passed whole
        seed: int = 1,
    ) -> None:
        stn_axis = np.asarray(stn_axis, dtype=float)
        if n_stn <= 0 or dt_s <= 0 or n_gpe_groups <= 0:
            raise ValueError("n_stn, dt_s, and n_gpe_groups must be positive")
        if len(stn_axis) != n_stn:
            raise ValueError("stn_axis must give one position per STN cell")
        if not 0.0 <= target_probability <= 1.0 or topo_sigma <= 0:
            raise ValueError("target_probability must be in [0, 1] and topo_sigma positive")
        if not 0.0 <= beta_floor <= 1.0 or not 0.0 <= beta_modulation <= 1.0:
            raise ValueError("beta_modulation and beta_floor must be in [0, 1]")
        if gpe_rate_tau_s <= 0 or feedback_delay_s < 0:
            raise ValueError("rate time constant must be positive; delay nonnegative")

        self.n_stn = n_stn
        self.dt_s = dt_s
        self.n_groups = n_gpe_groups
        self.rng = np.random.default_rng(seed)
        self.gpe_baseline_hz = gpe_baseline_hz
        self.stn_reference_hz = stn_reference_hz
        self.stn_to_gpe_gain = stn_to_gpe_gain
        self.gpe_rate_tau_s = gpe_rate_tau_s
        self.delay_steps = int(round(feedback_delay_s / dt_s))
        self.gpe_rate_max_hz = gpe_rate_max_hz
        self.step_index = 0

        # Group home positions along the STN axis and Gaussian affinity to each STN cell.
        self.group_center = (np.arange(self.n_groups) + 0.5) / self.n_groups
        affinity = np.exp(-(stn_axis[None, :] - self.group_center[:, None]) ** 2 / (2 * topo_sigma ** 2))  # groups x stn

        # GPe -> STN: sparse, local connections; every STN cell gets at least its nearest group.
        connected = self.rng.random(affinity.shape) < target_probability * affinity
        nearest = np.argmax(affinity, axis=0)
        connected[nearest, np.arange(n_stn)] |= ~connected.any(axis=0)
        self.targets_by_gpe_group = [np.flatnonzero(row).astype(np.int32) for row in connected]
        self.inputs_per_stn = connected.sum(axis=0)

        # STN -> GPe feedback: each group averages the STN rates with the same Gaussian weights.
        self.feedback_weights = affinity / affinity.sum(axis=1, keepdims=True)

        # Pallidal beta drive, depth graded along the axis by group position.
        self.beta_hz = beta_hz
        self.beta_depth = beta_floor + (beta_modulation - beta_floor) * np.exp(
            -(self.group_center - beta_center) ** 2 / (2 * beta_width ** 2))
        self.burst_start_prob = burst_rate_hz * dt_s
        self.burst_steps = int(round(burst_duration_s / dt_s))
        self.burst_steps_left = 0
        self.beta_phase = self.rng.uniform(0.0, 2.0 * np.pi) + self.rng.normal(0.0, beta_phase_jitter, self.n_groups)
        self.in_burst = False

        self.rate_hz = np.full(self.n_groups, gpe_baseline_hz, dtype=np.float64)
        self.emitted_rate_hz = self.rate_hz.copy()
        self.history = np.full((self.delay_steps + 1, self.n_groups), stn_reference_hz, dtype=np.float64)

    def step(self, stn_spiked_ids: np.ndarray | list[int]) -> np.ndarray:
        """Advance one time step and return STN IDs receiving GPe arrivals (see gpe_feedback.py)."""
        ids = np.asarray(stn_spiked_ids, dtype=np.int64)
        if ids.size and (ids.min() < 0 or ids.max() >= self.n_stn):
            raise ValueError("STN spike ID outside [0, n_stn)")

        if self.step_index >= self.delay_steps:
            delayed_stn_rate = self.history[(self.step_index - self.delay_steps) % len(self.history)]
        else:
            delayed_stn_rate = np.full(self.n_groups, self.stn_reference_hz)
        target_rate = self.gpe_baseline_hz + self.stn_to_gpe_gain * (delayed_stn_rate - self.stn_reference_hz)
        target_rate = np.clip(target_rate, 0.0, self.gpe_rate_max_hz)
        self.rate_hz += (1.0 - np.exp(-self.dt_s / self.gpe_rate_tau_s)) * (target_rate - self.rate_hz)

        # Beta burst envelope (shared by all groups) and group-specific modulation depth.
        if self.burst_steps_left > 0:
            self.burst_steps_left -= 1
        elif self.burst_steps > 0 and self.rng.random() < self.burst_start_prob:
            self.burst_steps_left = self.burst_steps
        self.in_burst = self.burst_steps_left > 0
        modulation = 1.0
        if self.in_burst:
            t = self.step_index * self.dt_s
            modulation = np.maximum(0.0, 1.0 + self.beta_depth * np.cos(2.0 * np.pi * self.beta_hz * t + self.beta_phase))
        self.emitted_rate_hz = self.rate_hz * modulation

        arrival_chunks = []
        for group_id, targets in enumerate(self.targets_by_gpe_group):
            counts = self.rng.poisson(self.emitted_rate_hz[group_id] * self.dt_s, size=len(targets))
            if np.any(counts):
                arrival_chunks.append(np.repeat(targets, counts))
        arrivals = np.concatenate(arrival_chunks).astype(np.int32, copy=False) if arrival_chunks else np.empty(0, dtype=np.int32)

        # Topographic feedback: Gaussian-weighted mean STN rate per group in this bin.
        stn_rates = np.bincount(ids, minlength=self.n_stn) / self.dt_s if ids.size else np.zeros(self.n_stn)
        self.history[self.step_index % len(self.history)] = self.feedback_weights @ stn_rates
        self.step_index += 1
        return arrivals
