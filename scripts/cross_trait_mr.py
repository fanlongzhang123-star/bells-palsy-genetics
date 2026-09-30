"""mr_sensitivity2.py — 修正版跨性状 MR 敏感性分析
剪枝: 按基因组顺序, 每个 win 窗口只保留最显著的 SNP (正确做法)
"""
import csv
import gzip
import math
import os

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT + r"\data"
W = _ROOT + r"\work"
OUT = _ROOT + r"\out"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
TRAITS = {
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
bp = {}
with gzip.open(os.path.join(D, "finngen_R13_G6_BELLPA.gz"), "rt", errors="replace") as f:
    next(f)
    for line in f:
        c = line.split("\t", 13)
        if len(c) < 13 or not c[4]:
            continue
        for rs in c[4].split(","):
            bp[rs.strip()] = (float(c[8]), float(c[9]), c[2].upper(), c[3].upper())
print(f"贝尔面瘫索引: {len(bp):,}")


def leads_win(code, win):
    """按基因组顺序, 每个 win 窗口保留最显著 SNP"""
    p = os.path.join(D, f"finngen_R13_{code}.gz")
    cand = []
    with gzip.open(p, "rt", errors="replace") as f:
        next(f)
        for line in f:
            c = line.split("\t", 13)
            if len(c) < 13 or not c[4]:
                continue
            try:
                pv = float(c[6])
            except ValueError:
                continue
            if pv < 5e-8:
                cand.append((c[0], int(c[1]), pv, c[4].split(",")[0].strip(),
                             c[2].upper(), c[3].upper(), float(c[8]), float(c[9])))
    cand.sort(key=lambda x: (x[0], x[1]))
    # 先按窗口分组, 每组取最小 p
    best = {}
    for ch, pos, pv, rs, ref, alt, b, se in cand:
        key = (ch, pos // win)
        if key not in best or pv < best[key][2]:
            best[key] = (ch, pos, pv, rs, ref, alt, b, se)
    out = sorted(best.values(), key=lambda x: (x[0], x[1]))
    return out


def mr(tools, exclude_mhc, n_win):
    bx, by, sx, sy = [], [], [], []
    n_ex = 0
    for ch, pos, pv, rs, ref, alt, b, se in tools:
        if exclude_mhc and ch == "6" and 25_000_000 <= pos <= 35_000_000:
            n_ex += 1
            continue
        f2 = bp.get(rs)
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
        return None
    w = [1.0 / (s / abs(b)) ** 2 for b, s in zip(bx, sy)]
    ratios = [yy / xx for xx, yy in zip(bx, by)]
    bm = sum(wi * r for wi, r in zip(w, ratios)) / sum(w)
    sem = math.sqrt(1 / sum(w))
    zz = bm / sem
    pm = math.erfc(abs(zz) / math.sqrt(2))
    return (n_win, rs and len(bx), n_ex, bm, sem, math.exp(bm), pm)


res = []
for win_mb in (1, 5):
    win = win_mb * 1_000_000
    print(f"\n{'='*104}")
    print(f"跨性状 MR —— 距离剪枝窗口 {win_mb} Mb (每窗口取最显著 SNP)")
    print("=" * 104)
    print(f"{'暴露(表型)':<24}{'全部SNP':>8}{'OR':>8}{'p':>11}   |"
          f"{'剔MHC后SNP':>11}{'OR':>8}{'p':>11}")
    for code, (name, nc) in TRAITS.items():
        tools = leads_win(code, win)
        r1 = mr(tools, False, win_mb)
        r2 = mr(tools, True, win_mb)
        if not r1:
            print(f"{name:<24} 工具不足")
            continue
        s2 = (f"{r2[1]:>11}{r2[5]:>8.3f}{r2[6]:>11.3g}" if r2
              else f"{'—':>11}{'—':>8}{'—':>11}")
        print(f"{name:<24}{r1[1]:>8}{r1[5]:>8.3f}{r1[6]:>11.3g}   |{s2}")
        if r1:
            res.append([win_mb, name, len(tools), r1[1], r1[5], r1[6], ""])
        if r2:
            res.append([win_mb, name + " [剔MHC]", len(tools), r2[1], r2[5], r2[6], r2[2]])

with open(os.path.join(OUT, "cross_trait_mr_sensitivity.csv"), "w", newline="",
          encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["prune_Mb", "exposure", "n_loci", "nsnp_used", "OR", "p", "n_excluded_MHC"])
    wr.writerows(res)
print(f"\n-> {os.path.join(OUT, 'cross_trait_mr_sensitivity.csv')}")
