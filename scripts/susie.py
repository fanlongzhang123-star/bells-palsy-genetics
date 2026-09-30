"""susie.py — SuSiE (Sum of Single Effects) 的 IBSS 实现, 用于汇总统计 (susie_rss 风格)

模型: z = sqrt(n) * R * b + e,  e ~ N(0, R)
      b = sum_l b_l, 每个 b_l 只有一个非零分量, 先验 b_lj ~ N(0, V_l)

SER 对数边际似然 (去掉与 j 无关的常数):
   lbf_j = 0.5 * n * z_j^2 / (n * R_jj + 1/V) - 0.5 * log(1 + n * V * R_jj)
后验: alpha_j ∝ exp(lbf_j);  mu_j = sqrt(n) * z_j * V / (1 + n * V * R_jj)
"""
import os
import math

import numpy as np

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def ser_update(z, R, n, V, Rjj=None):
    """单效应回归 (SER) 的变分更新"""
    p = len(z)
    if Rjj is None:
        Rjj = np.diag(R)
    denom = n * Rjj + 1.0 / V
    lbf = 0.5 * n * z ** 2 / denom - 0.5 * np.log1p(n * V * Rjj)
    lbf -= lbf.max()
    alpha = np.exp(lbf)
    s = alpha.sum()
    alpha = alpha / s if s > 0 else alpha
    mu = math.sqrt(n) * z * V / (1.0 + n * V * Rjj)
    return alpha, mu, float(np.log(s) + lbf.max())


def ser_loglik(z, R, n, V, Rjj=None):
    p = len(z)
    if Rjj is None:
        Rjj = np.diag(R)
    denom = n * Rjj + 1.0 / V
    lbf = 0.5 * n * z ** 2 / denom - 0.5 * np.log1p(n * V * Rjj)
    m = lbf.max()
    return float(m + np.log(np.sum(np.exp(lbf - m))))


def susie_rss(z, R, n, L=10, V_grid=None, max_iter=100, tol=1e-4, verbose=False):
    """IBSS 主循环。V 在每个成分上按网格最大化 SER 边际似然 (逐成分经验贝叶斯)。"""
    z = np.asarray(z, float)
    p = len(z)
    Rjj = np.diag(R).copy()
    if V_grid is None:
        V_grid = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2,
                           0.1, 0.2, 0.5])
    alpha = np.zeros((L, p))
    mu = np.zeros((L, p))
    V = np.zeros(L)
    for l in range(L):
        V[l] = V_grid[len(V_grid) // 2]

    # 初始化: 第一个成分用原始 z
    for l in range(L):
        best = None
        for Vv in V_grid:
            ll = ser_loglik(z, R, n, Vv, Rjj)
            if best is None or ll > best[0]:
                best = (ll, Vv)
        V[l] = best[1]
        a, m, _ = ser_update(z, R, n, V[l], Rjj)
        alpha[l], mu[l] = a, m
        break

    for it in range(max_iter):
        alpha_old = alpha.copy()
        for l in range(L):
            others = np.sum(alpha * mu, axis=0) - alpha[l] * mu[l]
            z_res = z - math.sqrt(n) * (R @ others)
            best = None
            for Vv in V_grid:
                ll = ser_loglik(z_res, R, n, Vv, Rjj)
                if best is None or ll > best[0]:
                    best = (ll, Vv)
            V[l] = best[1]
            a, m, _ = ser_update(z_res, R, n, V[l], Rjj)
            alpha[l], mu[l] = a, m
        delta = np.max(np.abs(alpha - alpha_old))
        if verbose:
            print(f"    iter {it+1}: max|dalpha| = {delta:.2e}, V = {np.round(V,5)}")
        if delta < tol:
            break
    return dict(alpha=alpha, mu=mu, V=V, Rjj=Rjj, n=n,
                pip=1 - np.prod(1 - alpha, axis=0))


def credible_sets(alpha, R, coverage=0.95, purity=0.5, min_purity=True):
    """从 alpha 提取可信集 (标准 SuSiE 规则: 累积到 coverage, 且纯度 >= purity)"""
    L, p = alpha.shape
    sets = []
    for l in range(L):
        a = alpha[l]
        if a.max() < 1e-6:
            continue
        order = np.argsort(-a)
        cum = np.cumsum(a[order])
        k = int(np.searchsorted(cum, coverage) + 1)
        cs = order[:k]
        # 纯度: 该集合内最小 r^2 (用 |r| 表示)
        sub = R[np.ix_(cs, cs)]
        rho = np.abs(sub)
        pur = rho[np.triu_indices(len(cs), 1)]
        minp = pur.min() if len(pur) else 1.0
        ok = (minp >= purity) if min_purity else True
        sets.append(dict(index=l, snps=cs, size=len(cs), purity=float(minp),
                         cumsum=float(cum[k - 1]), alpha_sum=float(a[cs].sum()),
                         ok=bool(ok)))
    return sets
