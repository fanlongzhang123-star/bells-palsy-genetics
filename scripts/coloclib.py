"""coloclib.py — 共定位与精细定位的共享实现 (coloc.abf 的 Python 移植)"""
import os
import math

import numpy as np

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def logsumexp(a):
    a = np.asarray(a, dtype=float)
    m = np.max(a)
    if not np.isfinite(m):
        return -np.inf
    return float(m + np.log(np.sum(np.exp(a - m))))


def abf_log(beta, se, sd_prior):
    """近似贝叶斯因子 (log)。coloc: r = V/(V+se^2), lABF = 0.5*(log(1-r) + r*z^2)"""
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    V = sd_prior ** 2
    r = V / (V + se ** 2)
    z2 = (beta / se) ** 2
    return 0.5 * (np.log(1 - r) + r * z2)


def coloc_abf(d1, d2, p1=1e-4, p2=1e-4, p12=1e-5, sd1=0.2, sd2=0.15):
    """d1/d2: dict(beta=..., se=...) 同一位点集。返回五个后验概率 + 每位点 H4 后验

    正确形式 (Giambartolomei 2014):
      H1 = Sigma_i a_i                -> A
      H2 = Sigma_j b_j                -> B
      H4 = Sigma_i a_i b_i            -> C = logsumexp(l1+l2)
      H3 = Sigma_{i!=j} a_i b_j       -> A + B + log(1 - exp(C-A-B))
    """
    l1 = abf_log(d1["beta"], d1["se"], sd1)
    l2 = abf_log(d2["beta"], d2["se"], sd2)
    A = logsumexp(l1)
    B = logsumexp(l2)
    C = logsumexp(l1 + l2)
    d = C - A - B                      # <= 0
    if d > -1e-12:
        l3_core = -np.inf              # 完全同一个变异 -> H3 = 0
    else:
        l3_core = A + B + math.log1p(-math.exp(d))
    ll = np.array([
        0.0,                                                    # H0
        A + math.log(p1) + math.log(1 - p2),                    # H1
        B + math.log(p2) + math.log(1 - p1),                    # H2
        l3_core + math.log(p1) + math.log(p2) + math.log(1 - p12),   # H3
        C + math.log(p1) + math.log(p2) + math.log(p12),         # H4
    ])
    mx = ll.max()
    pp = np.exp(ll - mx)
    pp /= pp.sum()
    w = np.exp(l1 + l2 - C)
    tot = w.sum()
    w = w / tot if tot > 0 else w
    return dict(PP_H0=float(pp[0]), PP_H1=float(pp[1]), PP_H2=float(pp[2]),
                PP_H3=float(pp[3]), PP_H4=float(pp[4]),
                nsnps=len(d1["beta"]), snp_pp=w)
