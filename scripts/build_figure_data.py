"""fig_data.py — 生成绘图用的 tidy 数据表 (英文标注, 供 R 读取)
输出到 figures/data/
"""
import csv, json, math, os

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT
OUT = D + r"\figures\data"
os.makedirs(OUT, exist_ok=True)

# phenocode -> 英文名 (与正文/表格一致)
NAME = {
    "G6_BELLPA": "Bell's palsy",
    "G6_NERPLEX": "Nerve, nerve root and plexus disorders",
    "G6_CARPTU": "Carpal tunnel syndrome",
    "G6_EPLEPSY": "Epilepsy",
    "AUTOIMMUNE_NONTHYROID": "Autoimmune diseases, non-thyroid",
    "E4_HYTHY_AI_STRICT": "Autoimmune hypothyroidism",
    "T1D_WIDE": "Type 1 diabetes",
    "M13_RHEUMA": "Rheumatoid arthritis",
    "L12_PSORIASIS": "Psoriasis",
    "G6_MS": "Multiple sclerosis",
    "T2D": "Type 2 diabetes",
}
CAT = {
    "G6_BELLPA": "target", "G6_NERPLEX": "entrapment", "G6_CARPTU": "entrapment",
    "G6_EPLEPSY": "neuro", "T2D": "metabolic",
    "AUTOIMMUNE_NONTHYROID": "autoimmune", "E4_HYTHY_AI_STRICT": "autoimmune",
    "T1D_WIDE": "autoimmune", "M13_RHEUMA": "autoimmune",
    "L12_PSORIASIS": "autoimmune", "G6_MS": "autoimmune",
}
ORDER = ["G6_BELLPA", "G6_NERPLEX", "G6_CARPTU", "G6_EPLEPSY", "AUTOIMMUNE_NONTHYROID",
         "E4_HYTHY_AI_STRICT", "T1D_WIDE", "M13_RHEUMA", "L12_PSORIASIS", "G6_MS", "T2D"]


def rd(p):
    with open(D + p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def wr(name, header, rows):
    with open(f"{OUT}\\{name}", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header); w.writerows(rows)
    print(f"  {name}: {len(rows)} 行")


# ================= 正态分位数 (Wichura AS241) =================
# p 可小到 1e-300, 故 z 不能由 1-p/2 反解 (会下溢为 1.0)
def _horner(x, c):
    v = c[0]
    for k in c[1:]:
        v = v * x + k
    return v


_A = [2.5090809287301226727e3, 3.3430575583588128105e4, 6.7265770927008700853e4,
      4.5921953931549871457e4, 1.3731693765509461125e4, 1.9715909503065514427e3,
      1.3314166789178437745e2, 3.3871328727963666080]
_B = [5.2264952788528545610e3, 2.8729085735721942674e4, 3.9307895800092710610e4,
      2.1213794301586595867e4, 5.3941960214247511077e3, 6.8718700749205790830e2,
      4.2313330701600911252e1, 1.0]
_C = [7.7454501427834140764e-4, 2.27238449892691845833e-2, 2.4178072517745061177e-1,
      1.27045825245236838258, 3.64784832476320460504, 5.7694972214606914055,
      4.6303378461565452959, 1.42343711074968357734]
_Dd = [1.05075007164441684324e-9, 5.475938084995344946e-4, 1.51986665636164571966e-2,
       1.4810397642748007459e-1, 6.8976733498510000455e-1, 1.6763848301838038494,
       2.05319162663775882187, 1.0]
_E = [2.01033439929228813265e-7, 2.71155556874348757815e-5, 1.2426609473880784386e-3,
      2.6532189526576123093e-2, 2.9656057182850489123e-1, 1.7848265399172913358,
      5.4637849111641143699, 6.6579046435011037772]
_F = [2.04426310338993978564e-15, 1.4215117583164458887e-7, 1.8463183175100546818e-5,
      7.868691311456132591e-4, 1.48753612908506148525e-2, 1.3692988092273580531e-1,
      5.9983220655588793769e-1, 1.0]


def _ppnd16(p):
    q = p - 0.5
    if abs(q) <= 0.425:
        r = 0.180625 - q * q
        return q * _horner(r, _A) / _horner(r, _B)
    r = p if q < 0 else 1.0 - p
    r = math.sqrt(-math.log(r))
    if r <= 5.0:
        val = _horner(r - 1.6, _C) / _horner(r - 1.6, _Dd)
    else:
        val = _horner(r - 5.0, _E) / _horner(r - 5.0, _F)
    return -val if q < 0 else val


def z_from_p(p):
    """双侧 p -> |z|"""
    return -_ppnd16(p / 2.0)


for _p, _z in [(0.05, 1.959963985), (0.5, 0.674489750), (0.3173105, 1.0),
               (0.001, 3.290526731), (1e-10, 6.466951087)]:
    assert abs(z_from_p(_p) - _z) / _z < 1e-5, (_p, z_from_p(_p), _z)


def _log_upper_tail(z):
    return -z * z / 2 - math.log(z) - 0.5 * math.log(2 * math.pi) \
        + math.log1p(-1 / z**2 + 3 / z**4 - 15 / z**6)


for _p in [1e-10, 4.2e-25, 1e-100, 1e-300]:
    _rel = abs(_log_upper_tail(z_from_p(_p)) - math.log(_p / 2)) / abs(math.log(_p / 2))
    assert _rel < 5e-4, (_p, _rel)
print(f"  [自校验] z(0.05)={z_from_p(0.05):.6f} z(0.5)={z_from_p(0.5):.6f} "
      f"z(1e-10)={z_from_p(1e-10):.4f} z(4.2e-25)={z_from_p(4.2e-25):.4f}")

# ================= Fig 3a: MHC 比例 =================
mhc = {r["phenocode"]: r for r in rd(r"\out\mhc_signal_all_traits.csv")}
wr("fig3a_mhc.tsv", ["trait", "phenocode", "category", "n_cases", "n_gws", "n_mhc", "pct_mhc"],
   [[NAME[pc], pc, CAT[pc], int(mhc[pc]["n_cases"]), int(mhc[pc]["n_gws_snp"]),
     int(mhc[pc]["n_gws_in_MHC"]), float(mhc[pc]["frac_in_MHC"]) * 100] for pc in ORDER])

# ================= Fig 3b: HLA 等位基因 =================
hla = {r["phenocode"]: r for r in rd(r"\out\hla_allele_association.csv")}
wr("fig3b_hla.tsv", ["trait", "phenocode", "category", "n_cases", "n_alleles",
                     "top_allele", "top_p", "n_lt_5e8"],
   [[NAME[pc], pc, CAT[pc], int(hla[pc]["n_cases"]), int(hla[pc]["n_alleles"]),
     hla[pc]["top_allele"], float(hla[pc]["top_p"]), int(hla[pc]["n_p_lt_5e8"])] for pc in ORDER])

# ================= Fig 4: 可信集位点效应 =================
# 来源 work/fig_cs_lookup.json (导出自 trait_cache.pkl; 与 LDSC与共定位_报告.md 第 4 节一致)
CS1_V = ["44478351", "44478937", "44479182", "44479267"]
CS2_V = "45602126"     # FinnGen 官方 SuSiE 的 CS2 lead (6:45,602,126:A:G)
P2N = {"44478351": "6:44478351", "44478937": "6:44478937",
       "44479182": "6:44479182", "44479267": "6:44479267", "45602126": "6:45602126"}
_look = json.load(open(D + r"\work\fig_cs_lookup.json", encoding="utf-8"))

# 交叉核对正文引用的关键数值
assert abs(_look["G6_BELLPA"]["44479182"][0] - 0.2108) < 1e-3
assert abs(_look["G6_BELLPA"]["44478351"][1] / 4.2e-25 - 1) < 0.05
assert abs(_look["G6_BELLPA"]["45602126"][0] - 0.1920) < 1e-3
assert abs(_look["G6_BELLPA"]["45602126"][1] / 4.13e-11 - 1) < 0.05
_cs1_other = [abs(_look[pc][v][0]) for pc in ORDER[1:] for v in CS1_V]
assert max(_cs1_other) <= 0.021, max(_cs1_other)          # 正文: 其余 |beta| <= 0.021
# 没有任何对照表型在该位点达到全基因组显著
assert min(_look[pc][v][1] for pc in ORDER[1:] for v in CS1_V) > 5e-8
# 逐一核对正文第 4a/4b 段引用的 beta 与 p 区间 (四舍五入到 4 位/2 位)
RANGE = {
    "G6_NERPLEX": ((0.0038, 0.0048), (0.41, 0.52)),
    "G6_CARPTU": ((-0.0191, -0.0185), (0.020, 0.025)),
    "G6_EPLEPSY": ((-0.0156, -0.0139), (0.17, 0.22)),
    "AUTOIMMUNE_NONTHYROID": ((-0.0038, -0.0032), (0.52, 0.58)),
    "E4_HYTHY_AI_STRICT": ((-0.0086, -0.0062), (0.21, 0.36)),
    "T1D_WIDE": ((-0.0197, -0.0186), (0.16, 0.19)),
    "M13_RHEUMA": ((-0.0034, -0.0028), (0.77, 0.81)),
    "L12_PSORIASIS": ((-0.0111, -0.0093), (0.37, 0.46)),
    "G6_MS": ((0.0200, 0.0207), (0.43, 0.44)),
    "T2D": ((-0.0032, -0.0027), (0.58, 0.64)),
}
for pc, ((bl, bh), (pl, ph)) in RANGE.items():
    for v in CS1_V:
        b, p = _look[pc][v]
        assert bl - 5e-5 <= b <= bh + 5e-5, (pc, v, b, bl, bh)
        assert pl - 5e-3 <= p <= ph + 5e-3, (pc, v, p, pl, ph)
_bp2 = _look["G6_BELLPA"]
assert all(0.208 <= _bp2[v][0] <= 0.211 for v in CS1_V)                  # 0.2081-0.2108
assert all(9.6e-26 <= _bp2[v][1] <= 4.2e-25 for v in CS1_V)              # 9.67e-26 - 4.16e-25
_oth2 = [abs(_look[pc][CS2_V][0]) for pc in ORDER[1:]]
assert max(_oth2) <= 0.0252, max(_oth2)          # 实测最大 0.0251 (多发性硬化)
assert all(0.129 <= _look[pc][CS2_V][1] <= 0.9525 for pc in ORDER[1:])   # 正文: 0.13 到 0.95
print("  [核对] BP CS1 beta 0.208-0.211 (p 9.7e-26 - 4.2e-25); "
      f"CS2 {P2N[CS2_V]} beta {_bp2[CS2_V][0]:.4f} p={_bp2[CS2_V][1]:.2e}; "
      f"其余 |beta|<={max(_cs1_other):.4f}; 10 个对照的 beta/p 区间全部与正文一致")

rows = []
for pc in ORDER:
    for v in CS1_V:
        b, p = _look[pc][v]
        rows.append([P2N[v], "CS1", NAME[pc], pc, CAT[pc], b,
                     abs(b) / z_from_p(p), p, -math.log10(p)])
for pc in ORDER:
    b, p = _look[pc][CS2_V]
    rows.append([P2N[CS2_V], "CS2", NAME[pc], pc, CAT[pc], b,
                 abs(b) / z_from_p(p), p, -math.log10(p)])
wr("fig4_credibleset.tsv", ["variant", "cs", "trait", "phenocode", "category",
                            "beta", "se", "p", "mlogp"], rows)

# ================= Fig 5: 表型组扫描 (来源 Table 4) =================
PH = [
    (1, "Bell's palsy", 6.10e-24, 0.2016, "target"),
    (2, "Spondylopathies", 5.30e-11, 0.0442, "skeletal"),
    (3, "Spinal stenosis", 1.23e-10, 0.0596, "skeletal"),
    (4, "Dorsopathies", 1.66e-09, 0.0268, "skeletal"),
    (5, "Other intervertebral disc disorders", 6.51e-06, 0.0303, "skeletal"),
    (6, "Coxarthrosis", 9.40e-06, 0.0387, "skeletal"),
    (7, "Early lumbar disc prolapse, operated", 1.70e-05, 0.0897, "skeletal"),
    (8, "Lumbar disc prolapse", 1.93e-05, 0.0337, "skeletal"),
    (9, "Dorsalgia", 2.27e-05, 0.0230, "skeletal"),
    (10, "Cervical disc disorders", 2.30e-05, 0.0458, "skeletal"),
    (11, "Coxarthrosis, primary, with hip surgery", 2.79e-05, 0.0419, "skeletal"),
    (12, "Autoimmune diseases, non-thyroid", 0.893, -0.0008, "autoimmune"),
    (13, "Nerve, nerve root and plexus disorders", 0.338, 0.0057, "entrapment"),
    (14, "Type 2 diabetes", 0.471, -0.0041, "metabolic"),
]
wr("fig5_phenome.tsv", ["rank", "trait", "p", "beta", "mlogp", "kind"],
   [[r, n, p, b, -math.log10(p), k] for r, n, p, b, k in PH])

# ================= Fig 6a: z 相关 =================
cor = {r["phenocode"]: r for r in rd(r"\out\autoimmune_boundary_summary.csv")}
wr("fig6a_zcorr.tsv", ["trait", "phenocode", "category", "r"],
   [[NAME[pc], pc, CAT[pc], float(cor[pc]["z_corr_with_BP"])] for pc in ORDER])

# ================= Fig 6b: 跨性状 MR (1 Mb 剪枝) =================
# 该文件只含 10 个对照暴露 (不含贝尔面瘫自身), 顺序稳定, 对应 ORDER[1:]
with open(D + r"\out\cross_trait_mr_sensitivity.csv", encoding="utf-8") as f:
    sens = list(csv.DictReader(f))
block = [r for r in sens if r["prune_Mb"] == "1" and "[剔MHC]" not in r["exposure"]]
assert len(block) == 10, len(block)
rows = []
for i, r in enumerate(block):
    pc = ORDER[i + 1]
    orv, p = float(r["OR"]), float(r["p"])
    se = abs(math.log(orv)) / z_from_p(p)
    rows.append([NAME[pc], pc, CAT[pc], int(r["nsnp_used"]), orv,
                 math.exp(math.log(orv) - 1.96 * se), math.exp(math.log(orv) + 1.96 * se),
                 se, p, -math.log10(p)])
wr("fig6b_crosstrait_mr.tsv", ["trait", "phenocode", "category", "nsnp", "OR",
                               "lo", "hi", "se_logOR", "p", "mlogp"], rows)

# 交叉核对正文数值
EXP = {"G6_NERPLEX": (35, 1.61, 1.8e-07), "G6_CARPTU": (85, 1.21, 2.9e-06),
       "T2D": (372, 1.09, 1.0e-04), "AUTOIMMUNE_NONTHYROID": (90, 1.12, 0.0016),
       "M13_RHEUMA": (53, 1.06, 0.0055), "G6_MS": (21, 0.95, 0.013)}
for row in rows:
    pc = row[1]
    if pc in EXP:
        n, o, p = EXP[pc]
        assert row[3] == n, (pc, row[3], n)
        assert abs(row[4] - o) < 0.006, (pc, row[4], o)
        assert abs(row[8] / p - 1) < 0.06, (pc, row[8], p)
print("  [核对] Fig6b 6 个正文引用的 nSNP/OR/p 全部一致")
# ================= 基因 / eQTL 轨道 =================
genes = json.load(open(D + r"\work\fig_genes.json", encoding="utf-8"))
eqtl = json.load(open(D + r"\work\fig_eqtl.json", encoding="utf-8"))
grows = []
for sym, g in sorted(genes.items(), key=lambda kv: kv[1]["start"]):
    e = eqtl.get(sym, {})
    grows.append([sym, g["start"], g["end"], g["strand"],
                  e.get("snp", "NA"), e.get("z", "NA"), e.get("pos", "NA")])
wr("fig_genes.tsv", ["symbol", "start", "end", "strand",
                     "eqtl_snp", "eqtl_z", "eqtl_pos"], grows)
print(f"  [核对] 蛋白编码基因 {len(grows)} 个, 其中有 eQTL 的 "
      f"{sum(1 for r in grows if r[5] != 'NA')} 个 (正文: 16 个基因)")

# ================= 可信集成员位置 =================
CS_POS = {"CS1": [44478351, 44478937, 44479182, 44479267], "CS2": [45602126]}
wr("fig_crediblesets.tsv", ["cs", "pos", "lead"],
   [["CS1", 44478351, 0], ["CS1", 44478937, 0], ["CS1", 44479182, 1],
    ["CS1", 44479267, 0], ["CS2", 45602126, 1], ["CS3", 45960517, 1]])

print("\n=== Fig 6b ===")
for row in rows:
    print(f"  {row[0]:42s} nSNP={row[3]:4d} OR={row[4]:.3f} "
          f"({row[5]:.3f}-{row[6]:.3f}) p={row[8]:.3g}")
