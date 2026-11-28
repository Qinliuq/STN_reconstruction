"""Cortical->STN topography variants (results/ctx_topography), selected with cfg.ctx_topo_variant (None = off).

Each variant gives a per-cell arrival-rate scale and an optional per-cell conductance scale, both normalised to mean 1
over the model cells, so the population-mean cortical input is unchanged. The same variant applies in healthy and PD.
When cfg.ctx_topo_variant is set it replaces cfg.ctx_topography; it cannot be combined with cfg.ctx_gscale.

  T0  rate (0.65 - 0.5 x)                          the current default rule (bit-identical to cfg.ctx_topography)
  T1  uniform rate
  T2  uniform rate; conductance (0.568 - 0.273 x)  re-digitised Jeon 2022 Fig. 7H.2 Gria1 linear fit
                                                   (params/jeon2022_profiles.json; results/emergent_peak/step01.md)
  T3  rate max(T_motor(x), 0.3)                    T_motor: Plantinga 2018 motor territory (stn_drive_topography,
                                                   human intervals, smoothed, territories summing to 1); 0.3 = floor
                                                   for non-motor cortex
"""
import numpy as np

T3_FLOOR = 0.3
VARIANTS = {
    'T0': dict(rate=('linear', 0.65, -0.5), g=None),
    'T1': dict(rate=None, g=None),
    'T2': dict(rate=None, g=('linear', 0.568, -0.273)),
    'T3': dict(rate=('motor', T3_FLOOR), g=None),
}


def _shape(axis, spec):
    axis = np.asarray(axis, float)
    if spec is None:
        return np.ones(len(axis))
    if spec[0] == 'linear':
        v = spec[1] + spec[2] * axis
    elif spec[0] == 'motor':
        import stn_drive_topography as D
        v = np.maximum(D.territories(axis)['motor'], spec[1])
    else:
        raise ValueError(f'unknown topography spec {spec}')
    if np.any(v <= 0):
        raise ValueError(f'topography {spec} gives a non-positive scale')
    return v / v.mean()


def rate_scale(axis, variant):
    return _shape(axis, VARIANTS[variant]['rate'])


def g_scale(axis, variant):
    spec = VARIANTS[variant]['g']
    return None if spec is None else _shape(axis, spec)
