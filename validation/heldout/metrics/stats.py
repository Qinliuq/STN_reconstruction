"""Distribution comparisons and contingency statistics shared by the comparison scripts."""
import numpy as np
from scipy import stats as sstats


def compare_distributions(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    mw = sstats.mannwhitneyu(a, b, alternative='two-sided')
    ks = sstats.ks_2samp(a, b)
    return dict(n_a=len(a), n_b=len(b), median_a=float(np.median(a)), median_b=float(np.median(b)),
                mw_U=float(mw.statistic), mw_p=float(mw.pvalue), ks_D=float(ks.statistic), ks_p=float(ks.pvalue))


def ecdf(x):
    x = np.sort(np.asarray(x, float))
    x = x[np.isfinite(x)]
    return x, np.arange(1, len(x) + 1) / len(x)


def cramers_v(a, b):
    """Cramér's V and chi-square / Fisher p for two categorical vectors; returns dict incl. the table."""
    import pandas as pd
    tab = pd.crosstab(pd.Series(a, name='a'), pd.Series(b, name='b'))
    T = tab.values
    out = dict(table=tab.to_dict())
    if min(T.shape) < 2:
        return dict(out, V=np.nan, p=np.nan)
    chi2, p, _, _ = sstats.chi2_contingency(T, correction=False)
    n = T.sum()
    out.update(V=float(np.sqrt(chi2 / (n * (min(T.shape) - 1)))), chi2_p=float(p))
    if T.shape == (2, 2):
        out['fisher_p'] = float(sstats.fisher_exact(T).pvalue)
    return out


def circ_diff(a, b):
    return float(np.angle(np.exp(1j * (a - b))))
