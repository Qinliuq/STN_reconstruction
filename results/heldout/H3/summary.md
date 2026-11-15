# H3

**H3: MATCH — slow PPC match (model median 0.0862 vs rat 0.0806); beta PPC to cortical input not testable (model 0.00337 vs rat activated-epoch nan, 0 rat units qualify); beta to GPe drive not testable (model 0.00367). Secondary, slow-wave epochs: beta mismatch (model PD+SO 0.0043 vs rat all-epoch -4.62e-05).**

Prediction written 2026-11-15T13:20:22 (git 63ec255d); compared 2026-11-15T13:20:31.

## Summary table

| source | group | n | median_ppc | q25 | q75 | frac_rayleigh_p01 | circ_mean_pref_rad |
|---|---|---|---|---|---|---|---|
| rat slow | all | 37 | 0.0806 | 0.0234 | 0.13 | 0.973 | 2.21 |
| rat beta (activated) | all | 0 | nan | nan | nan | nan | nan |
| rat beta (all epochs) | all | 37 | -4.62e-05 | -0.000249 | 0.000201 | 0.108 | -0.787 |
| rat slow | bursty | 13 | 0.194 | 0.105 | 0.299 | 1 | -0.386 |
| rat beta (activated) | bursty | 0 | nan | nan | nan | nan | nan |
| rat beta (all epochs) | bursty | 13 | 0.000179 | -0.000196 | 0.000814 | 0.308 | -2.45 |
| rat slow | regular | 24 | 0.0402 | 0.0189 | 0.0947 | 0.958 | 2.53 |
| rat beta (activated) | regular | 0 | nan | nan | nan | nan | nan |
| rat beta (all epochs) | regular | 24 | -0.000109 | -0.000266 | 0.000125 | 0 | 0.323 |
| model slow (PD+SO, cortical input) | all | 290 | 0.0862 | 0.0497 | 0.117 | 0.993 | -0.16 |
| model slow (PD+SO, cortical input) | PV+ | 60 | 0.0367 | 0.032 | 0.0461 | 1 | -0.106 |
| model slow (PD+SO, cortical input) | PV- | 230 | 0.0981 | 0.0738 | 0.129 | 0.991 | -0.174 |
| model slow (PD+SO, cortical input) | bursty | 213 | 0.0941 | 0.0709 | 0.12 | 0.991 | -0.177 |
| model slow (PD+SO, cortical input) | regular | 77 | 0.0413 | 0.0335 | 0.0591 | 1 | -0.113 |
| model beta (cortical input) | all | 290 | 0.00337 | 0.000466 | 0.00813 | 0.169 | 1.36 |
| model beta (cortical input) | PV+ | 60 | 0.0011 | -0.00034 | 0.00205 | 0 | 1.41 |
| model beta (cortical input) | PV- | 230 | 0.0052 | 0.0011 | 0.00943 | 0.213 | 1.35 |
| model beta (cortical input) | bursty | 213 | 0.00527 | 0.000916 | 0.00957 | 0.202 | 1.36 |
| model beta (cortical input) | regular | 77 | 0.00139 | -7.98e-05 | 0.00295 | 0.0779 | 1.34 |
| model beta (GPe drive) | all | 290 | 0.00367 | 8.03e-05 | 0.00839 | 0.224 | -1.32 |
| model beta (GPe drive) | PV+ | 60 | 0.00436 | 0.00217 | 0.00682 | 0.383 | -1.18 |
| model beta (GPe drive) | PV- | 230 | 0.00338 | -0.000803 | 0.00907 | 0.183 | -1.37 |
| model beta (GPe drive) | bursty | 213 | 0.00362 | -0.000389 | 0.00943 | 0.192 | -1.4 |
| model beta (GPe drive) | regular | 77 | 0.0037 | 0.00167 | 0.00618 | 0.312 | -1.11 |
| model beta (GPe drive, in bursts) | all | 290 | 0.0148 | 0.00326 | 0.0331 | 0.428 | -1.35 |
| model beta (GPe drive, in bursts) | PV+ | 60 | 0.0204 | 0.0142 | 0.0307 | 0.817 | -1.24 |
| model beta (GPe drive, in bursts) | PV- | 230 | 0.0108 | 0.00223 | 0.0345 | 0.326 | -1.38 |
| model beta (GPe drive, in bursts) | bursty | 213 | 0.0135 | 0.00246 | 0.0412 | 0.347 | -1.41 |
| model beta (GPe drive, in bursts) | regular | 77 | 0.0159 | 0.00608 | 0.0266 | 0.649 | -1.19 |
| model slow, no SO (control) | all | 290 | 0.000661 | -0.00091 | 0.00305 | 0.0241 | -0.872 |
| model slow, no SO (control) | PV+ | 60 | 0.00115 | -6.06e-05 | 0.00248 | 0.0167 | -1.01 |
| model slow, no SO (control) | PV- | 230 | 0.000536 | -0.00137 | 0.00315 | 0.0261 | -0.835 |
| model slow, no SO (control) | bursty | 213 | 0.000231 | -0.00149 | 0.0029 | 0.0235 | -0.761 |
| model slow, no SO (control) | regular | 77 | 0.00146 | -1.24e-05 | 0.00312 | 0.026 | -1.15 |
| model beta in PD+SO runs (cortical input) | all | 290 | 0.0043 | 0.00089 | 0.0082 | 0.193 | 1.17 |
| model beta in PD+SO runs (cortical input) | PV+ | 60 | 0.0012 | -0.000565 | 0.0029 | 0.0167 | 1.26 |
| model beta in PD+SO runs (cortical input) | PV- | 230 | 0.00525 | 0.00171 | 0.00927 | 0.239 | 1.15 |
| model beta in PD+SO runs (cortical input) | bursty | 213 | 0.00509 | 0.00152 | 0.00934 | 0.211 | 1.17 |
| model beta in PD+SO runs (cortical input) | regular | 77 | 0.00218 | -0.000325 | 0.00412 | 0.143 | 1.17 |
| model beta in PD+SO runs (GPe drive) | all | 290 | 0.00274 | -0.00022 | 0.00729 | 0.19 | -1.5 |
| model beta in PD+SO runs (GPe drive) | PV+ | 60 | 0.00341 | 0.000949 | 0.00614 | 0.283 | -1.28 |
| model beta in PD+SO runs (GPe drive) | PV- | 230 | 0.00242 | -0.000874 | 0.00869 | 0.165 | -1.58 |
| model beta in PD+SO runs (GPe drive) | bursty | 213 | 0.00262 | -0.000903 | 0.00936 | 0.178 | -1.59 |
| model beta in PD+SO runs (GPe drive) | regular | 77 | 0.00283 | 0.000347 | 0.00594 | 0.221 | -1.27 |

## Activated-cortex epochs (rat)

Rule: 4 s windows with slow fraction < 0.4. Per session (windows within the unit spike span):

| file | units | windows | activated |
|---|---|---|---|
| NPR040d.mat | 6 | 23 | 0 |
| NPR046.mat | 4 | 23 | 0 |
| NPR052e.mat | 4 | 24 | 0 |
| NPR057.mat | 9 | 24 | 0 |
| NPR059.mat | 5 | 24 | 0 |
| pjx310c.mat | 9 | 24 | 0 |

0 of 37 units with EEG have >= 50 spikes in activated windows. 8 units (pjx286c, pjx289c) have no EEG. The unit trains are 100 s epochs chosen during slow-wave activity; the activated windows of NPR040d, NPR057 and pjx310c fall outside them, so the activated-epoch beta test cannot be run on these data.

## Tests

- slow: match; MW p 0.396, KS p 0.0198, medians model 0.0862 / rat 0.0806
- beta (activated epochs), cortical input: not testable
- beta (activated epochs), GPe drive: not testable
- secondary: rat beta all epochs vs model beta in PD+SO runs, cortical input: mismatch; MW p 8.49e-09, KS p 7.47e-15, medians model 0.0043 / rat -4.62e-05
- secondary: same, GPe drive: mismatch; MW p 9.23e-05, KS p 8.87e-10, medians model 0.00274 / rat -4.62e-05

## Caveats

The model has no slow-wave cortex: the slow oscillation is an imposed ~1 Hz up/down envelope on the cortical input (priors.json: depth 0.9, duty 0.5, mean preserved), so slow locking is largely set by that choice. The model reference is the cortical input rate (phase 0 = peak input) or the GPe beta drive (phase 0 = peak GPe rate); the rat reference is surface EEG, whose polarity relative to cortical firing is not established here, so preferred phases are compared for concentration, not absolute angle. Rat EEG is mostly slow-wave dominated; activated epochs exist in 3 of 6 sessions. Anesthesia: urethane (rat).

