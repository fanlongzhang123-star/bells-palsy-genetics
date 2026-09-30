"""coloc_mr3.py — 优化版: 每个表型文件只扫一次, 同时收集
  (a) 贝尔面瘫位点区域 chr6:44.0-46.0 Mb
  (b) MHC 区 chr6:30.5-33.5 Mb
  (c) SuSiE 可信集位点
  (d) 全基因组显著 lead SNP (用于跨性状 MR)
然后做: 共定位 (BP vs 各表型) + 可信集查询 + 跨性状 MR
"""
import csv
import gzip
import math
import os
import pickle
import sys
import time

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT + r"\data"
W = _ROOT + r"\work"
OUT = _ROOT + r"\out"
CACHE = os.path.join(W, "trait_cache.pkl")

TRAITS = {
    "G6_BELLPA": ("贝尔面瘫", 5062),
    "G6_NERPLEX": ("神经/神经根疾患", 68559),
    "G6_CARPTU": ("腕管综合征", 32951),
    "G6_EPLEPSY": ("癫痫", 16295),
    "AUTOIMMUNE_NONTHYROID": ("自身免疫复合", 72412),
    "E4_HYTHY_AI_STRICT": ("自身免疫性甲减", 56098),
    "T1D_WIDE": ("1型糖尿病", 11197),
    "M13_RHEUMA": ("类风湿关节炎", 16775),
    "L12_PSORIASIS": ("银屑病", 13309),
    "G6_MS": ("多发性硬化", 3002),
    "T2D": ("2型糖尿病", 89727),
}
BP_REG = ("6", 44_000_000, 46_000_000)
MHC_REG = ("6", 30_500_000, 33_500_000)
CS_POS = {"44478351", "44478937", "44479182", "44479267", "45600864", "45602126"}
P1 = P2 = 1e-4
P12 = 1e-5
WVAR = 0.15 ** 2


def scan_trait(path, do_leads=True):
    bp, mhc, cs = {}, {}, {}
    leads, used = [], {}
    with gzip.open(path, "rt", errors="replace") as f:
        next(f)
        for line in f:
            c = line.split("\t", 13)
            if len(c) < 13:
                continue
            try:
                pv = float(c[6]); b = float(c[8]); se = float(c[9])
            except ValueError:
                continue
            if se <= 0:
                continue
            if c[0] == "6":
                try:
                    pos = int(c[1])
                except ValueError:
                    pos = -1
                if BP_REG[1] <= pos <= BP_REG[2] or MHC_REG[1] <= pos <= MHC_REG[2]:
                    rec = dict(pos=pos, ref=c[2].upper(), alt=c[3].upper(),
                               beta=b, se=se, af=float(c[10]), p=pv)
                    if c[4]:
                        rs = c[4].split(",")[0].strip()
                        (bp if pos <= BP_REG[2] else mhc)[rs] = rec
                    if c[1] in CS_POS:
                        cs[c[1]] = (b, pv)
                if do_leads and pv < 5e-8 and c[4]:
                    rs = c[4].split(",")[0].strip()
                    lk = used.get(c[0])
                    if lk is None or abs(pos - lk) >= 1_000_000:
                        used[c[0]] = pos
                        leads.append((c[0], pos, rs, c[2].upper(), c[3].upper(), b, se))
    return dict(bp=bp, mhc=mhc, cs=cs, leads=leads)


def logsumexp(v):
    if not v:
        return -math.inf
    m = max(v)
    if m == -math.inf:
        return -math.inf
    return m + math.log(sum(math.exp(x - m) for x in v))


def logdiffexp(a, b):
    return -math.inf if b >= a else a + math.log1p(-math.exp(b - a))


def log_abf(beta, se):
    v = se * se
    r = WVAR / (WVAR + v)
    return 0.5 * math.log(1 - r) + 0.5 * r * (beta / se) ** 2


def align(rec, a_alt, a_ref):
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
    e, o = rec["alt"], rec["ref"]
    if e == a_alt and o == a_ref:
        return rec["beta"]
    if e == a_ref and o == a_alt:
        return -rec["beta"]
    if COMP.get(e) == a_alt and COMP.get(o) == a_ref:
        return rec["beta"]
    if COMP.get(e) == a_ref and COMP.get(o) == a_alt:
        return -rec["beta"]
    return None


def coloc(r1, r2, label, quiet=False):
    la1, la2 = [], []
    for rs, a in r1.items():
        b = r2.get(rs)
        if not b or a["af"] < 0.01 or a["af"] > 0.99:
            continue
        b2 = align(b, a["alt"], a["ref"])
        if b2 is None:
            continue
        la1.append(log_abf(a["beta"], a["se"]))
        la2.append(log_abf(b2, b["se"]))
    n = len(la1)
    if n < 20:
        if not quiet:
            print(f"  {label:<30} 可比 SNP {n}, 跳过", flush=True)
        return None
    L1 = logsumexp(la1) - math.log(n)
    L2 = logsumexp(la2) - math.log(n)
    L4 = logsumexp([x + y for x, y in zip(la1, la2)]) - math.log(n)
    L3 = logdiffexp(L1 + L2, L4)
    prior = [1 - P1 - P2 - P12, P1, P2, P12, P12]
    terms = [math.log(p) + l for p, l in zip(prior, [0.0, L1, L2, L3, L4])]
    z = logsumexp(terms)
    pp = [math.exp(t - z) for t in terms]
    if not quiet:
        print(f"  {label:<30} n={n:<6} H0={pp[0]:.3f} H1={pp[1]:.3f} "
              f"H2={pp[2]:.3f} H3={pp[3]:.3f} H4={pp[4]:.3f}", flush=True)
    return dict(label=label, n=n, H0=pp[0], H1=pp[1], H2=pp[2], H3=pp[3], H4=pp[4])


if __name__ == "__main__":
    if os.path.exists(CACHE):
        cache = pickle.load(open(CACHE, "rb"))
        print(f"从缓存载入 {len(cache)} 个表型")
    else:
        cache = {}
        for code in TRAITS:
            p = os.path.join(D, f"finngen_R13_{code}.gz")
            if not os.path.exists(p):
                continue
            t0 = time.time()
            cache[code] = scan_trait(p)
            d = cache[code]
            print(f"  {code}: BP区 {len(d['bp'])}, MHC区 {len(d['mhc'])}, "
                  f"lead {len(d['leads'])}  ({time.time()-t0:.0f}s)", flush=True)
        pickle.dump(cache, open(CACHE, "wb"))
        print("缓存已保存\n")

    # ---- 验证: 自己 vs 自己 ----
    print("=" * 96)
    print("校验 (自己 vs 自己, 应得 H4≈1)")
    print("=" * 96)
    coloc(cache["G6_BELLPA"]["bp"], cache["G6_BELLPA"]["bp"], "BP vs BP [应当 H4]")

    rows = []
    for rname, key in [("BP_locus(chr6:44-46Mb)", "bp"), ("MHC(chr6:30.5-33.5Mb)", "mhc")]:
        print("\n" + "=" * 96)
        print(f"共定位: 贝尔面瘫 vs 各表型  [{rname}]")
        print("=" * 96)
        for code, (name, nc) in TRAITS.items():
            if code == "G6_BELLPA":
                continue
            r = coloc(cache["G6_BELLPA"][key], cache[code][key], f"BP vs {name}")
            if r:
                r["region"] = rname
                rows.append(r)
    with open(os.path.join(OUT, "coloc_bp_vs_traits.csv"), "w", newline="",
              encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["region", "pair", "n_snp", "PPH0", "PPH1", "PPH2", "PPH3", "PPH4"])
        for r in rows:
            wr.writerow([r["region"], r["label"], r["n"], r["H0"], r["H1"],
                         r["H2"], r["H3"], r["H4"]])
    print(f"\n-> {os.path.join(OUT, 'coloc_bp_vs_traits.csv')}")

    # ---- 可信集查询 ----
    print("\n" + "=" * 96)
    print("FinnGen SuSiE 可信集 SNP 在各表型中的效应")
    print("  CS1 = 4 个基因间 SNP (chr6:44,478,351-44,479,267); CS2 = RUNX2 内含子")
    print("=" * 96)
    codes = list(TRAITS)
    print(f"{'位点':<12}{'CS':>4}" + "".join(f"{TRAITS[c][0][:7]:>10}" for c in codes))
    for pos in ["44478351", "44478937", "44479182", "44479267", "45600864", "45602126"]:
        csname = "CS1" if int(pos) < 45_000_000 else "CS2"
        cells = []
        for c in codes:
            v = cache.get(c, {}).get("cs", {}).get(pos)
            cells.append(f"{v[0]:>10.4f}" if v else f"{'—':>10}")
        print(f"{pos:<12}{csname:>4}" + "".join(cells), flush=True)
    print("\n(p 值)")
    print(f"{'位点':<12}{'CS':>4}" + "".join(f"{TRAITS[c][0][:7]:>10}" for c in codes))
    for pos in ["44478351", "44478937", "44479182", "44479267", "45600864", "45602126"]:
        csname = "CS1" if int(pos) < 45_000_000 else "CS2"
        cells = []
        for c in codes:
            v = cache.get(c, {}).get("cs", {}).get(pos)
            cells.append(f"{v[1]:>10.3g}" if v else f"{'—':>10}")
        print(f"{pos:<12}{csname:>4}" + "".join(cells), flush=True)

    # ---- 跨性状 MR ----
    print("\n" + "=" * 96)
    print("跨性状孟德尔随机化: 各表型独立 lead SNP -> 贝尔面瘫")
    print("=" * 96)
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
    bp_idx = {}
    for rs, rec in cache["G6_BELLPA"]["bp"].items():
        bp_idx[rs] = rec
    print(f"贝尔面瘫区域索引 {len(bp_idx)} 个 SNP (仅 chr6:44-46Mb; MR 需全基因组索引)")
    # 用缓存里的 lead SNP 信息不完整(需要 BP 全基因组), 改为单独构建
    bp_full = {}
    with gzip.open(os.path.join(D, "finngen_R13_G6_BELLPA.gz"), "rt",
                   errors="replace") as f:
        next(f)
        for line in f:
            c = line.split("\t", 13)
            if len(c) < 13 or not c[4]:
                continue
            for rs in c[4].split(","):
                bp_full[rs.strip()] = (float(c[8]), float(c[9]), c[2].upper(), c[3].upper())
    print(f"贝尔面瘫全基因组索引: {len(bp_full):,}")
    res = []
    print(f"\n{'暴露(表型)':<26}{'nSNP':>6}{'IVW_b':>10}{'se':>9}{'OR':>8}{'p':>11}")
    for code, (name, nc) in TRAITS.items():
        if code == "G6_BELLPA":
            continue
        tools = cache[code]["leads"]
        bx, by, sx, sy = [], [], [], []
        for ch, pos, rs, ref, alt, b, se in tools:
            f2 = bp_full.get(rs)
            if not f2:
                continue
            bb, sb, refB, altB = f2
            if alt == altB and ref == refB:
                pass
            elif alt == refB and ref == altB:
                bb = -bb
            elif COMP.get(alt) == altB and COMP.get(ref) == refB:
                pass
            elif COMP.get(alt) == refB and COMP.get(ref) == altB:
                bb = -bb
            else:
                continue
            bx.append(b); sx.append(se); by.append(bb); sy.append(sb)
        if len(bx) < 2:
            print(f"{name:<26}{len(bx):>6}  工具不足")
            continue
        w = [1.0 / (s / abs(b)) ** 2 for b, s in zip(bx, sy)]
        ratios = [yy / xx for xx, yy in zip(bx, by)]
        bm = sum(wi * r for wi, r in zip(w, ratios)) / sum(w)
        sem = math.sqrt(1 / sum(w))
        z = bm / sem
        pm = math.erfc(abs(z) / math.sqrt(2))
        print(f"{name:<26}{len(bx):>6}{bm:>10.4f}{sem:>9.4f}{math.exp(bm):>8.3f}{pm:>11.3g}",
              flush=True)
        res.append([name, code, len(bx), bm, sem, math.exp(bm), pm])
    with open(os.path.join(OUT, "cross_trait_mr_to_BP.csv"), "w", newline="",
              encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["exposure", "phenocode", "nsnp", "ivw_b", "se", "OR", "p"])
        wr.writerows(res)
    print(f"\n-> {os.path.join(OUT, 'cross_trait_mr_to_BP.csv')}")
