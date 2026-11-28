# Supplementary result: K_ATP is not sufficient for PV+ burst/pause firing (2026-11-18)

**PD-upregulated K_ATP is not sufficient for burst/pause firing in the data-constrained PV+ cell.**

In configuration D (PV+ Kv3.1 1.2, Kir 0.3), K_ATP at 0.5-5e-5 S/cm2 (k_Ca 20 uM) only slows
pacemaking, from 22 to 5-6 Hz, and it stays tonic (CV <= 0.18). Across the Kv3.1 x K_ATP map
(Kv3.1 0.9-1.2 at Kir 0.3; `results/katp_kv31_map/`), no cell meets the pre-registered burst/pause
rule. At most, K_ATP produces regular doublets (CV 0.7-0.9) at Kv3.1 <= 1.15 and gbar >= 3e-5.

K_ATP bursting needs **both** PV+ Kv3.1 <= 1.0 and Kir >= 0.5, the pre-configuration-D cell. There it
gives 3-spike bursts at ~5 Hz with CV 1.1, though even these fall short of the Poisson-surprise
criterion. That cell conflicts with depolarization-block protection, which needs Kv3.1 >= 1.10-1.15.
At Kv3.1 1.0 its spikes fail to repolarize at +0.8 nA (trough -36 mV). Lowering Kir alone to 0.3 at
Kv3.1 1.0 also abolishes the bursting, so the conflict is not with Kv3.1 alone. The configuration-D
Kir equalization removes it too.

K_ATP also abolishes the -0.16 nA rhythmic bursting from 3e-5 at every Kv3.1, and leaves the
rebound burst intact.

Consequences:
- No network screen was run. Configuration D is unchanged, and K_ATP stays an optional,
  default-off, PD-only mechanism (`cfg.katp_gbar`, `cfg.katp_state_scale`; bit-identical when off).
- The held-out H1 phenotype mismatch stands: the model PD state has no burst/pause phenotype, and
  this intrinsic mechanism does not supply one within the constrained cell.

Details: `README.md` (configuration-D sweep), `../katp_kv31_map/README.md` (diagnostic and map).
