TITLE katp_stn.mod : Ca2+/NO-gated ATP-sensitive K+ current (Kir6.2/SUR1) for STN neurons

COMMENT
-----------------------------------------------------------------------------
ATP-sensitive K+ current (K_ATP) for subthalamic nucleus projection neurons.

RATIONALE FOR THE GATING CHOICE
    In STN neurons, K_ATP is not usefully modelled as a pure function of bulk
    [ATP]. The direct experimental evidence is that it is recruited by a
    Ca2+ -> nNOS -> NO -> cGMP cascade downstream of NMDA receptor activation:

      Shen K-Z, Johnson SW. Ca2+ influx through NMDA-gated channels activates
      ATP-sensitive K+ currents through a nitric oxide-cGMP pathway in
      subthalamic neurons. J Neurosci 30:1882-1893 (2010).

      Shen K-Z, Johnson SW. Chronic dopamine depletion augments the functional
      expression of K-ATP channels in the rat subthalamic nucleus.
      Neurosci Lett 531:104-108 (2012). doi:10.1016/j.neulet.2012.10.030

      Shen K-Z, Johnson SW. AMP kinase regulates K-ATP currents evoked by
      NMDA receptor stimulation in rat subthalamic nucleus neurons.
      Neuroscience 274:138-152 (2014).

    This mechanism therefore follows the approach of Knowlton et al. (2018,
    J Neurophysiol, doi:10.1152/jn.00351.2017; ModelDB 227678), who drove
    K_ATP open probability from intracellular Ca2+ rather than from bulk ATP
    in midbrain dopamine neurons.

    Subunit note: the pore-forming subunit modelled here is Kir6.2 (KCNJ11),
    which with SUR1 forms the neuronal K_ATP channel. Kir6.1 (KCNJ8) is the
    principal pore-forming subunit of ASTROCYTIC, not neuronal, plasma-membrane
    K_ATP (Thomzig et al. 2001) and should not be used as a somatic neuronal
    conductance.

STRUCTURE
    o            slow activation gate, driven by cai through the NO/cGMP lag
    fatp         metabolic availability factor (ATP block relieved by ADP)
    pdscale      multiplier for the parkinsonian condition. Shen & Johnson
                 (2012) report larger NMDA-evoked outward current and greater
                 conductance on the 6-OHDA-lesioned side, so pdscale > 1
                 represents chronic dopamine depletion. 1.0 = control.

    g  = gbar * pdscale * fatp * o
    ik = g * (v - ek)

    Note that fatp multiplies gbar, so the effective conductance is
    gbar * fatp * o, not gbar. With the defaults below fatp ~ 0.53, i.e. about
    half the channel population is available for Ca2+-dependent recruitment.

    The current is ohmic. Kir6.2 is only weakly inward-rectifying, so this is
    the usual simplification; if rectification matters for your question, add a
    Boltzmann or Mg2+-block factor to g.

CAVEAT ON THE SIGN OF THE EFFECT
    In STN, opening K_ATP SUPPRESSES firing, and Shen & Johnson (2012) propose
    K_ATP activation as potentially beneficial in parkinsonism. This is the
    opposite direction from a generic inward-rectifier (Kir2.x-like) block.
    Do not treat this mechanism as a drop-in reinterpretation of a "Kir block"
    result without re-running the simulations.

DEFAULT PARAMETERS ARE PLACEHOLDERS. kca, tau and gbar in particular are not
tightly constrained by STN data and must be fitted or swept.

Units: gbar S/cm2, cai and nucleotides mM, tau ms.
-----------------------------------------------------------------------------
ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (mM) = (milli/liter)
}

NEURON {
    SUFFIX katp_stn
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gbar, pdscale, g
    RANGE omin, kca, nca, tau
    RANGE atp, adp, kati, kadp, hatp
    RANGE oinf, fatp
}

PARAMETER {
    gbar    = 1e-4  (S/cm2)  : peak conductance density; fit or sweep
    pdscale = 1.0            : >1 for chronic dopamine depletion

    : --- Ca2+/NO-driven gate ---
    omin    = 0.005          : basal open probability at resting Ca2+
    kca     = 5e-4  (mM)     : half-activating [Ca]i (0.5 uM); placeholder
    nca     = 2              : cooperativity of Ca2+ drive
    tau     = 300   (ms)     : lumped NO/cGMP cascade lag; placeholder

    : --- metabolic availability (ATP block, relieved by ADP) ---
    atp     = 2.5   (mM)     : intracellular ATP
    adp     = 0.05  (mM)     : intracellular ADP
    kati    = 1.0   (mM)     : apparent intact-cell K1/2 for ATP block. The
                             : reported intact-cell range is ~0.1-1 mM.
                             : NB excised-patch values (10-30 uM) shut the
                             : channel completely at physiological ATP: at
                             : ATP = 2.5 mM, kati = 0.02 gives fatp ~ 5e-4.
    kadp    = 0.03  (mM)     : [ADP] that doubles the apparent ATP K1/2
    hatp    = 2              : Hill coefficient for ATP block
}

ASSIGNED {
    v       (mV)
    ek      (mV)
    cai     (mM)
    ik      (mA/cm2)
    g       (S/cm2)
    oinf
    fatp
}

STATE { o }

INITIAL {
    rates()
    o = oinf
    fatp = avail()
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    fatp = avail()
    g  = gbar * pdscale * fatp * o
    ik = g * (v - ek)
}

DERIVATIVE states {
    rates()
    o' = (oinf - o) / tau
}

PROCEDURE rates() {
    LOCAL drive
    : Hill activation by intracellular Ca2+ (safe at cai = 0)
    drive = (cai^nca) / (cai^nca + kca^nca)
    oinf  = omin + (1 - omin) * drive
}

FUNCTION avail() {
    : Fraction of channels not blocked by ATP. ADP competitively relieves
    : block by raising the apparent K1/2. For an explicit MgADP-activated
    : formulation instead, substitute the Magnus & Keizer (1998) expression.
    avail = 1 / (1 + (atp / (kati * (1 + adp/kadp)))^hatp)
}
