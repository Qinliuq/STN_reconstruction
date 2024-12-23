TITLE Non-inactivating inwardly rectifying potassium current (Kir2.3)


NEURON {
    SUFFIX Kir
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik, shift
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gk = 0.0001 (mho/cm2)
    shift = 0
    q = 1 	: body temperature 35 C
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    winf
    wtau (ms)
}

STATE { w }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gk*w*(v-ek)
}

DERIVATIVE states {
    rates()
    w' = (winf-w)/wtau*q
}

INITIAL {
    rates()
    w = winf
}

PROCEDURE rates() {
    UNITSOFF
    winf = 1/(1+exp((v-(-82)-shift)/13))
    wtau = 1/(exp((v-(-103)-shift)/(-14.5))+0.125/(1+exp((v-(-35)-shift)/(-19))))
    UNITSON
}



COMMENT

Original model by Wolf et al (2005) [1] for the rat MSN cells from the
nucleus accumbens.  The activation curve was fitted to a mouse Kir2.1
channel expressed in HEK cells [2] and shifted to match extracellular
concentration of K in rat.  Time constants were derived from Aplysia data
[3] and adjusted to match the rat experiments [1]. Time constant was
further tuned [4] to fit the rat data below -80 mV [5].  Kinetics is
corrected to the body temperature 35 C.

Non-inactivating Kir current was observed in cells expressing Kir2.2
and/or Kir2.3 [5]. Activation variable with m^1 kinetics is used [1,4].
Smooth fit of the time constants by Alexander Kozlov <akozlov@kth.se>.

[1] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

[2] Kubo Y, Murata Y (2001) Control of rectification and permeation by
two distinct sites after the second transmembrane region in Kir2.1 K+
channel. J Physiol 531, 645-660.

[3] Hayashi H, Fishman HM (1988) Inward rectifier K+ channel kinetics
from analysis of the complex conductance of aplysia neuronal membrane.
Biophys J 53, 747-757.

[4] Steephen JE, Manchanda R (2009) Differences in biophysical properties
of nucleus accumbens medium spiny neurons emerging from inactivation of
inward rectifying potassium currents. J Comput Neurosci 27(3):453-70

[5] Mermelstein PG, Song WJ, Tkatch T, Yan Z, Surmeier DJ (1998)
Inwardly rectifying potassium (IRK) currents are correlated with IRK
subunit expression in rat nucleus accumbens medium spiny neurons. J
Neurosci 18(17):6650-61.

ENDCOMMENT
