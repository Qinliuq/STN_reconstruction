"""
Territory map x denervation weights for the STN PD-change gradient.
Produces: T_c(x), D(x), drive scaling s(x), and the conjunction with PV+ density.
x: principal STN axis (axis1), 0 = ventromedial/anteromedial, 1 = dorsolateral/posterolateral.
"""
import numpy as np

x = np.arange(0.05, 1.0, 0.10)
gauss = lambda x, mu, sd: np.exp(-0.5*((x-mu)/sd)**2)
hill  = lambda D, n, K: D**n/(D**n + K**n)

# --- territory layout (section 3 of the doc) -------------------------------
VARIANTS = {
    "primate_human": dict(limbic=(0.08, 0.14), assoc=(0.38, 0.20), motor=(0.78, 0.22)),
    "mouse_jeon":    dict(limbic=(0.10, 0.16), assoc=(0.35, 0.22), motor=(0.62, 0.26)),
}
# --- channel denervation weights (section 2) -------------------------------
W = {"advanced": dict(limbic=0.45, assoc=0.80, motor=0.965),
     "early":    dict(limbic=0.20, assoc=0.35, motor=0.60)}
# --- PV+ density along axis1 (Jeon: dorsolateral AND middle) ---------------
PV = {"jeon_dl_mid": gauss(x, 0.68, 0.22),
      "monotonic":   1/(1+np.exp(-(x-0.5)/0.12))}       # control

def territories(variant):
    T = {c: gauss(x, *VARIANTS[variant][c]) for c in VARIANTS[variant]}
    S = sum(T.values())
    return {c: T[c]/S for c in T}                        # sum_c T_c(x) = 1

def drive(variant="primate_human", stage="advanced", n=2, K=0.70):
    T = territories(variant)
    D = sum(W[stage][c]*T[c] for c in T)
    s = hill(D, n, K)
    return D, s/s.mean()                                 # mean 1 == same total as uniform

def ratio(p, lo=0.25, hi=0.65):
    f = lambda t: np.interp(t, x, p)
    return f(hi)/f(lo)

if __name__ == "__main__":
    for variant in VARIANTS:
        for stage in W:
            for n, K in [(1, .70), (2, .70), (4, .85)]:
                D, s = drive(variant, stage, n, K)
                for pvname, pv in PV.items():
                    pvn = pv/pv.mean()
                    prod = s*pvn; prod /= prod.mean()
                    print(f"{variant:<14}{stage:<9}n={n} K={K:.2f} PV={pvname:<12}"
                          f"drive peak x={x[s.argmax()]:.2f} r={ratio(s):.2f} | "
                          f"product peak x={x[prod.argmax()]:.2f} r={ratio(prod):.2f}")
