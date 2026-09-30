"""
mrlib.py — 纯 numpy/标准库实现的孟德尔随机化工具集
（本环境无法安装 pandas/scipy，故自行实现统计分布与 MR 估计量）

包含:
  - 特殊函数: 正则化不完全 gamma / beta (用于 chi2 与 t 检验 p 值)
  - qnorm (Acklam 近似)
  - 等位基因协调 (harmonise)
  - MR 估计量: Wald ratio, IVW(固定/随机), MR-Egger, 加权中位数, 加权众数
  - 敏感性: Cochran Q 异质性, Egger 截距(多效性), 留一法, 离群点检测
  - Steiger 方向性
"""
import os
import math
import numpy as np

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---------------------------------------------------------------- 特殊函数
def _gser(a, x, itmax=1000, eps=3e-16):
    ap, s, d = a, 1.0 / a, 1.0 / a
    for _ in range(itmax):
        ap += 1.0
        d *= x / ap
        s += d
        if abs(d) < abs(s) * eps:
            break
    return s * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a, x, itmax=1000, eps=3e-16, fpmin=1e-300):
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def gammq(a, x):
    """正则化上不完全 gamma 函数 Q(a,x)"""
    if x < 0 or a <= 0:
        raise ValueError("gammq domain")
    if x == 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)


def chi2_sf(x, df):
    """chi2 生存函数 (=p 值, 单侧上尾)"""
    if df <= 0:
        return float("nan")
    return gammq(df / 2.0, x / 2.0)


def _betacf(a, b, x, itmax=500, eps=3e-16, fpmin=1e-300):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def betai(a, b, x):
    """正则化不完全 beta 函数 I_x(a,b)"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t, df):
    if df <= 0 or not np.isfinite(t):
        return float("nan")
    return betai(df / 2.0, 0.5, df / (df + t * t))


def norm_sf_two_sided(z):
    """标准正态双侧 p 值 (TwoSampleMR 对 IVW / Wald / 加权中位数 采用正态近似)"""
    if not np.isfinite(z):
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2.0))


def qnorm(p):
    """标准正态分位数 (Acklam 近似, 精度 ~1e-9)"""
    if p <= 0.0 or p >= 1.0:
        return float("nan")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


# ---------------------------------------------------------------- 数据协调
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonise(bx, by, ea_x, oa_x, ea_y, oa_y, sx=None, sy=None,
              rsid=None, action=2):
    """
    按效应等位基因统一方向。
    返回 dict of numpy arrays (仅保留可协调的 SNP), 并统计各 action 类别。
    action=2 (TwoSampleMR 默认): 丢弃回文且 MAF 不可判断的 SNP。
    """
    n = len(bx)
    keep, act = [], []
    for i in range(n):
        ex, ox = str(ea_x[i]).upper(), str(oa_x[i]).upper()
        ey, oy = str(ea_y[i]).upper(), str(oa_y[i]).upper()
        if ex == ey and ox == oy:
            a = 0
        elif ex == oy and ox == ey:
            a = 1          # 需要翻转 outcome
        elif ex in _COMP and _COMP[ex] == ey:
            a = 2          # 回文 (palindromic)
        else:
            a = 3          # 等位基因不兼容
        keep.append(a)
        act.append(a)
    keep = np.array(keep)
    act = np.array(act)
    ok = keep == 0
    flip = keep == 1
    bx, by = np.asarray(bx, float).copy(), np.asarray(by, float).copy()
    by[flip] = -by[flip]
    ea_y_arr, oa_y_arr = np.asarray(ea_y).copy(), np.asarray(oa_y).copy()
    ea_y_arr[flip], oa_y_arr[flip] = oa_y_arr[flip].copy(), ea_y_arr[flip].copy()
    use = (keep == 0) | (keep == 1)
    out = dict(bx=bx[use], by=by[use],
               ea=np.asarray(ea_x)[use], oa=np.asarray(oa_x)[use],
               ea_y=ea_y_arr[use], oa_y=oa_y_arr[use])
    if sx is not None:
        out["sx"] = np.asarray(sx, float)[use]
    if sy is not None:
        out["sy"] = np.asarray(sy, float)[use]
    if rsid is not None:
        out["rsid"] = np.asarray(rsid)[use]
    out["n_palindromic_dropped"] = int((act == 2).sum())
    out["n_incompatible_dropped"] = int((act == 3).sum())
    out["n_flipped"] = int(flip.sum())
    out["n_used"] = int(use.sum())
    return out


# ---------------------------------------------------------------- MR 估计量
def _ratio(by, bx, sy):
    b = by / bx
    se = np.abs(sy / bx)
    return b, se


def mr_wald(by, bx, sy):
    b, se = _ratio(np.array([by]), np.array([bx]), np.array([sy]))
    p = norm_sf_two_sided(b[0] / se[0])
    return dict(method="Wald ratio", nsnp=1, b=float(b[0]), se=float(se[0]),
                p=float(p), lo=float(b[0] - 1.96 * se[0]),
                hi=float(b[0] + 1.96 * se[0]))


def mr_ivw(by, bx, sy, random=True):
    b, se = _ratio(by, bx, sy)
    w = 1.0 / se ** 2
    b_ivw = float(np.sum(w * b) / np.sum(w))
    se_fe = float(math.sqrt(1.0 / np.sum(w)))
    q = float(np.sum(w * (b - b_ivw) ** 2))
    n = len(b)
    df = n - 1
    q_p = chi2_sf(q, df) if df > 0 else float("nan")
    phi = max(1.0, q / df) if df > 0 else 1.0
    se_use = se_fe * math.sqrt(phi) if random else se_fe
    p = norm_sf_two_sided(b_ivw / se_use)
    return dict(method="IVW" + (" (random)" if random else " (fixed)"),
                nsnp=n, b=b_ivw, se=se_use, se_fixed=se_fe, p=float(p),
                lo=b_ivw - 1.96 * se_use, hi=b_ivw + 1.96 * se_use,
                Q=q, Q_df=df, Q_p=q_p)


def mr_egger(by, bx, sy):
    """加权最小二乘 by ~ bx (权重 1/sy^2)，截距=多效性"""
    y = np.asarray(by, float)
    x = np.asarray(bx, float)
    w = 1.0 / np.asarray(sy, float) ** 2
    n = len(x)
    X = np.column_stack([np.ones(n), x])
    W = np.diag(w)
    XtW = X.T @ W
    cov = np.linalg.pinv(XtW @ X)
    beta = cov @ (XtW @ y)
    resid = y - X @ beta
    dof = n - 2
    sigma2 = float(resid @ W @ resid / dof) if dof > 0 else float("nan")
    se = np.sqrt(np.diag(cov) * sigma2)
    return dict(method="MR Egger", nsnp=n,
                b=float(beta[1]), se=float(se[1]),
                p=float(t_sf_two_sided(beta[1] / se[1], dof)),
                lo=float(beta[1] - 1.96 * se[1]), hi=float(beta[1] + 1.96 * se[1]),
                intercept=float(beta[0]), intercept_se=float(se[0]),
                intercept_p=float(t_sf_two_sided(beta[0] / se[0], dof)))


def _weighted_median(b, w):
    order = np.argsort(b)
    b, w = b[order], w[order]
    w = w / np.sum(w)
    cw = np.cumsum(w) - 0.5 * w
    below = np.where(cw < 0.5)[0]
    if len(below) == 0:
        return float(b[0])
    k = int(below[-1])
    if k >= len(b) - 1:
        return float(b[-1])
    d = cw[k + 1] - cw[k]
    if d <= 0:
        return float(b[k])
    return float(b[k] + (b[k + 1] - b[k]) * (0.5 - cw[k]) / d)


def mr_weighted_median(by, bx, sy, nboot=2000, seed=20260914):
    b, se = _ratio(by, bx, sy)
    w = 1.0 / se ** 2
    est = _weighted_median(b.copy(), w.copy())
    rng = np.random.default_rng(seed)
    n = len(b)
    boot = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        se_b = se[idx]
        if np.sum(se_b == 0) > 0 or np.any(~np.isfinite(se_b)):
            boot[i] = np.nan
            continue
        wb = 1.0 / se_b ** 2
        boot[i] = _weighted_median(b[idx], wb)
    se_est = float(np.nanstd(boot, ddof=1))
    return dict(method="Weighted median", nsnp=n, b=est, se=se_est,
                p=float(norm_sf_two_sided(est / se_est)),
                lo=est - 1.96 * se_est, hi=est + 1.96 * se_est)


def _weighted_mode(b, w, phi=1.0):
    n = len(b)
    sd = np.std(b, ddof=1)
    iqr = np.subtract(*np.percentile(b, [75, 25]))
    denom = min(sd, iqr / 1.349) if iqr > 0 else sd
    bw = phi * 0.9 * denom / (n ** 0.2)
    if bw <= 0 or not np.isfinite(bw):
        bw = max(sd, 1e-8)
    dens = np.array([np.sum(w * np.exp(-0.5 * ((bj - b) / bw) ** 2)) for bj in b])
    return float(b[int(np.argmax(dens))])


def mr_weighted_mode(by, bx, sy, nboot=2000, seed=20260914):
    b, se = _ratio(by, bx, sy)
    w = 1.0 / se ** 2
    est = _weighted_mode(b.copy(), w.copy())
    rng = np.random.default_rng(seed)
    n = len(b)
    boot = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        boot[i] = _weighted_mode(b[idx], 1.0 / se[idx] ** 2)
    se_est = float(np.nanstd(boot, ddof=1))
    return dict(method="Weighted mode", nsnp=n, b=est, se=se_est,
                p=float(norm_sf_two_sided(est / se_est)),
                lo=est - 1.96 * se_est, hi=est + 1.96 * se_est)


def mr_rucker(by, bx, sy, sx=None):
    """Rücker 框架: 不含截距的加权回归 (by - b*bx 的残差方差), 用于 Q'/I²"""
    out = mr_ivw(by, bx, sy, random=False)
    q = out["Q"]
    df = out["Q_df"]
    ip = max(0.0, (q - df) / q) if q > 0 else 0.0
    return dict(Q=q, Q_df=df, Q_p=out["Q_p"], I2=ip)


def mr_presso_like(by, bx, sy, nsim=2000, seed=20260914):
    """MR-PRESSO 风格的全局检验 + 离群点检测(简化: 全局用 Cochran Q; 离群用学生化残差)"""
    b, se = _ratio(by, bx, sy)
    w = 1.0 / se ** 2
    b_ivw = float(np.sum(w * b) / np.sum(w))
    res = b - b_ivw
    stud = res / se
    n = len(b)
    # 全局检验: Cochran Q
    q = float(np.sum(w * res ** 2))
    p_global = chi2_sf(q, n - 1)
    # 离群点: |学生化残差| 超过 Bonferroni 校正的 t 临界值
    alpha = 0.05 / n
    # 双侧临界 t
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if t_sf_two_sided(mid, n - 2) > alpha:
            lo = mid
        else:
            hi = mid
    crit = (lo + hi) / 2
    outlier_idx = np.where(np.abs(stud) > crit)[0]
    return dict(global_Q=q, global_p=float(p_global),
                outlier_idx=outlier_idx.tolist(),
                outlier_crit=float(crit),
                studentised=stud.tolist())


def leave_one_out(by, bx, sy):
    """留一法 IVW"""
    n = len(by)
    rows = []
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        r = mr_ivw(np.asarray(by)[m], np.asarray(bx)[m], np.asarray(sy)[m],
                   random=True)
        rows.append((i, r["b"], r["se"], r["p"]))
    return rows


def steiger_r(b, se, p, n):
    """TwoSampleMR 的 get_r_from_bsen / get_r_from_pn: 返回 r 及 r2"""
    out = []
    for bi, sei, pi, ni in zip(b, se, p, n):
        if pi is not None and np.isfinite(pi) and 0 < pi < 1 and ni:
            t = abs(qnorm(pi / 2.0))
            r2 = t * t / (t * t + ni - 2)
            out.append(math.sqrt(max(r2, 0.0)))
        elif sei and ni:
            r = bi / math.sqrt(bi * bi + sei * sei * ni)
            out.append(abs(r))
        else:
            out.append(float("nan"))
    return np.array(out)


def f_statistic(bx, sx, n, eaf=None):
    """工具变量强度 F; 若给 eaf 用 R2 公式, 否则用 (b/se)^2"""
    bx = np.asarray(bx, float)
    sx = np.asarray(sx, float)
    t2 = (bx / sx) ** 2
    if eaf is not None:
        eaf = np.asarray(eaf, float)
        r2 = 2 * eaf * (1 - eaf) * bx ** 2
        k = len(bx)
        f = (r2 / k) / ((1 - r2) / (n - k - 1))
        return f, r2
    return t2, None


def format_result(r):
    return (f"{r['method']:<22s} n={r['nsnp']:>4d}  "
            f"b={r['b']:+.4f} (se={r['se']:.4f})  "
            f"OR={math.exp(r['b']):.3f} [{math.exp(r['lo']):.3f}, {math.exp(r['hi']):.3f}]  "
            f"p={r['p']:.4g}")
