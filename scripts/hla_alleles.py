"""hla_analysis.py — HLA 等位基因层面关联对比 (自身免疫最敏感的检测)
FinnGen 提供各表型的 HLA 插补关联结果; 若贝尔面瘫有自身免疫基础,
其经典 HLA 等位基因 (DRB1*15:01, DQB1*02:01, A*02:01 等) 应显示关联。
"""
import csv
import glob
import gzip
import os

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT + r"\data\finngen_aux"
OUT = _ROOT + r"\out"

LABEL = {
    "G6_BELLPA": ("贝尔面瘫", 5062, "目标"),
    "G6_NERPLEX": ("神经/神经根/神经丛疾患", 68559, "神经-近"),
    "G6_CARPTU": ("腕管综合征", 32951, "神经-卡压"),
    "G6_EPLEPSY": ("癫痫", 16295, "神经-对照"),
    "AUTOIMMUNE_NONTHYROID": ("自身免疫病(非甲状腺)复合", 72412, "自身免疫"),
    "E4_HYTHY_AI_STRICT": ("自身免疫性甲减", 56098, "自身免疫"),
    "T1D_WIDE": ("1型糖尿病", 11197, "自身免疫"),
    "M13_RHEUMA": ("类风湿关节炎", 16775, "自身免疫"),
    "L12_PSORIASIS": ("银屑病", 13309, "自身免疫"),
    "G6_MS": ("多发性硬化", 3002, "自身免疫"),
    "T2D": ("2型糖尿病", 89727, "代谢"),
}

results = []
print(f"{'表型':<28}{'类别':<10}{'检验数':>7}{'最强等位':>14}{'最强p':>11}"
      f"{'p<5e-8':>8}{'p<Bonf':>8}")
for code, (name, nc, cat) in LABEL.items():
    p = os.path.join(D, f"HLA_{code}.gz")
    if not os.path.exists(p):
        print(f"{name:<28} 文件缺失 ({code})")
        continue
    rows = []
    with gzip.open(p, "rt", errors="replace") as f:
        hdr = next(f).rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(hdr):
                continue
            try:
                pv = float(c[idx["pval"]])
            except ValueError:
                continue
            rows.append((pv, c[idx["alt"]], float(c[idx["beta"]])))
    n = len(rows)
    if n == 0:
        continue
    rows.sort()
    bonf = 0.05 / n
    n5 = sum(1 for r in rows if r[0] < 5e-8)
    nb = sum(1 for r in rows if r[0] < bonf)
    best = rows[0]
    print(f"{name:<28}{cat:<10}{n:>7}{best[1]:>14}{best[0]:>11.3g}{n5:>8}{nb:>8}")
    results.append([name, code, cat, nc, n, best[1], best[0], best[2], n5, nb])

with open(os.path.join(OUT, "hla_allele_association.csv"), "w", newline="",
          encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["trait", "phenocode", "category", "n_cases", "n_alleles",
                 "top_allele", "top_p", "top_beta", "n_p_lt_5e8", "n_p_lt_bonferroni"])
    wr.writerows(results)
print(f"\n-> {os.path.join(OUT, 'hla_allele_association.csv')}")

# 经典自身免疫 HLA 等位基因在各表型中的表现
KEY = ["DRB1*15:01", "DRB1*03:01", "DRB1*04:01", "DQB1*02:01", "DQB1*03:02",
        "DQA1*05:01", "A*02:01", "B*57:01", "C*06:02", "DPB1*04:01", "B*27:05",
        "DRB1*07:01", "DRB1*09:01", "DQB1*06:02"]
print("\n" + "=" * 100)
print("经典自身免疫 HLA 等位基因在各表型中的 p 值")
print("=" * 100)
hdr = f"{'等位基因':<14}" + "".join(f"{LABEL[c][0][:8]:>11}" for c in LABEL
                                     if os.path.exists(os.path.join(D, f"HLA_{c}.gz")))
print(hdr)
table = {}
codes = [c for c in LABEL if os.path.exists(os.path.join(D, f"HLA_{c}.gz"))]
for c in codes:
    table[c] = {}
    with gzip.open(os.path.join(D, f"HLA_{c}.gz"), "rt", errors="replace") as f:
        h = next(f).rstrip("\n").split("\t")
        ia = h.index("alt"); ip = h.index("pval")
        for line in f:
            x = line.rstrip("\n").split("\t")
            if len(x) > max(ia, ip):
                table[c][x[ia]] = float(x[ip])
for a in KEY:
    cells = []
    for c in codes:
        pv = table[c].get(a)
        cells.append(f"{pv:>11.3g}" if pv is not None else f"{'—':>11}")
    print(f"{a:<14}" + "".join(cells))
