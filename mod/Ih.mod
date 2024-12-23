TITLE Potassium Ih channel for STh (HCN)

COMMENT
 Ih from Huguenard & McCormick 92 LGN relay
 Magee (1998) shows Ih to have a reversal potential that is ionic non
 specific with a reversal potential of around -30mV.  
 Huguenard recordings at 35.5degC. 

 How the q10 works: There is a q10 for the rates (alpha and beta's)
 called Q10 and a Q10 for the maximum conductance called gmaxQ10.  The
 q10s should have been measured at specific temperatures temp1 and
 temp2 (that are 10degC apart). Ideally, as Q10 is temperature
 dependant, we should know these two temperatures.  We used to
 follow the more formal Arrhenius derived Q10 approach.  The
 temperature at which this channel's kinetics were recorded is tempb
 (base temperature).  What we then need to calculate is the desired
 rate scale for now working at temperature celsius (rate_k).  This was
 given by the empirical Arrhenius equation, using the Q10, but now is 
 using the quick Q10 approximation. 
ENDCOMMENT


UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Ih
	NONSPECIFIC_CURRENT ih
	RANGE gk
	GLOBAL eih,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gk    = 0.001 (mho/cm2)
	eih   = -5.611047394e+01 (mV)
	celsius

	activate_Q10 = 1
	Q10 = 2.0
	gmaxQ10 = 2.0
	temp1 = 25.0 (degC)
	temp2 = 35.0 (degC)
	tempb = 35.5 (degC)
}

STATE {
        f
}

ASSIGNED { 
        ih (mA/cm2)
	finf  
	ftau (ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE integrate METHOD cnexp
	ih = (gk*gmax_k)*f*(v-eih)
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  rate_k = Q10^((celsius-tempb)/10)
          gmax_k = gmaxQ10^((celsius-tempb)/10)
	}else{
	  rate_k = 1.0
	  gmax_k = 1.0
	}
        setinf(v)
	f = finf
}

DERIVATIVE integrate {
        setinf(v)
	f' = (finf - f)/ftau
}

PROCEDURE setinf(v) {
                    :Voltage shift (for temp effects) of 5.0.
	TABLE finf, ftau DEPEND celsius FROM -100 TO 100 WITH 400

	finf = 1.0/(1+exp((v + 80 )/ 5.5))
        ftau = (1.0/(exp(-15.02 - 0.086*v)+exp(-1.5195 + 0.0701*v))) /rate_k
}

UNITSON
