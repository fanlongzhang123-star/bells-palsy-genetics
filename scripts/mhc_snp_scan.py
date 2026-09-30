"""mhc_scan_all.py — 统计各表型 MHC 区(chr6:25-35Mb) 的全基因组显著 SNP 数与占比"""
import csv
import gzip
import os
import time

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT + r"\data"
OUT = _ROOT + r"\out"
MHC = ("6", 25_000_000, 35_000_000)

TRAITS = [
    ("G6_BELLPA", "贝尔面瘫", 5062, "目标"),
    ("G6_NERPLEX", "神经/神经根/神经丛疾患", 68559, "神经-近"),
    ("G6_CARPTU", "腕管综合征(卡压性)", 32951, "神经-卡压"),
    ("G6_EPLEPSY", "癫痫", 16295, "神经-对照"),
    ("AUTOIMMUNE_NONTHYROID", "自身免疫病(非甲状腺)复合", 72412, "自身免疫"),
    ("E4_HYTHY_AI_STRICT", "自身免疫性甲减", 56098, "自身免疫"),
    ("T1D_WIDE", "1型糖尿病", 11197, "自身免疫"),
    ("M13_RHEUMA", "类风湿关节炎", 16775, "自身免疫"),
    ("L12_PSORIASIS", "银屑病", 13309, "自身免疫"),
    ("G6_MS", "多发性硬化", 3002, "自身免疫"),
    ("T2D", "2型糖尿病", 89727, "代谢"),
]

rows = []
print(f"{'表型':<30}{'病例':>8}{'显著SNP':>10}{'MHC内':>9}{'占比':>9}{'MHC最强p':>13}{'MHC外最强p':>14}")
for code, name, nc, cat in TRAITS:
    p = os.path.join(D, f"finngen_R13_{code}.gz")
    if not os.path.exists(p):
        print(f"{name:<30} 文件缺失")
        continue
    t0 = time.time()
    n_gws = n_mhc = 0
    mhc_best = None
    non_best = None
    with gzip.open(p, "rt", errors="replace") as f:
        next(f)
        for line in f:
            c = line.split("\t", 13)
            if len(c) < 13:
                continue
            try:
                pv = float(c[6])
            except ValueError:
                continue
            in_mhc = (c[0] == "6" and MHC[1] <= int(c[1]) <= MHC[2])
            if pv < 5e-8:
                n_gws += 1
                if in_mhc:
                    n_mhc += 1
            if in_mhc:
                if mhc_best is None or pv < mhc_best:
                    mhc_best = pv
            else:
                if non_best is None or pv < non_best:
                    non_best = pv
    frac = n_mhc / n_gws if n_gws else 0
    print(f"{name:<30}{nc:>8}{n_gws:>10}{n_mhc:>9}{frac*100:>8.1f}%{mhc_best:>13.3g}{non_best:>14.3g}")
    rows.append([name, code, cat, nc, n_gws, n_mhc, round(frac, 4), mhc_best, non_best])

with open(os.path.join(OUT, "mhc_signal_all_traits.csv"), "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["trait", "phenocode", "category", "n_cases", "n_gws_snp",
                 "n_gws_in_MHC", "frac_in_MHC", "mhc_best_p", "nonmhc_best_p"])
    wr.writerows(rows)
print(f"\n-> {os.path.join(OUT, 'mhc_signal_all_traits.csv')}")
