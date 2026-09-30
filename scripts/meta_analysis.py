"""run_r13_and_meta.py — 三步走
  (A) FinnGen R13 G6_BELLPA (5,062 例)  — 最大单队列窄表型
  (B) FinnGen R13 G6_NERPLEX (68,559 例) — 宽表型 (神经/神经根/神经丛疾患)
  (C) 跨队列 meta: R13 BELLPA + GCST90018843 (1,894 例欧洲)
每个都输出: 全基因组显著位点 -> 基因集 -> CAM 通路富集 -> PTPRC/ICAM1 检查
"""
import csv
import gzip
import json
import os
import time
import urllib.request

import build_bp_geneset as B

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT + r"\data"
W = _ROOT + r"\work"
OUT = _ROOT + r"\out"
OLD_VCF = os.environ.get(
    "GCST90018843_VCF",
    os.path.join(_ROOT, "data", "GCST90018843.vcf.gz"))

CAM14 = ["CADM3", "CDH4", "CLDN16", "CNTN2", "CNTNAP2", "CTLA4", "GLG1", "ITGA9",
         "NRXN1", "NRXN3", "PTPRC", "PTPRD", "PTPRM", "PVRL3"]
WATCH = CAM14 + ["ICAM1", "ITGB2", "ITGAL", "SELE", "VCAM1", "PECAM1", "NECTIN3",
                 "CDC5L", "CLIC5", "TCF12"]

_g38 = None


def genes38():
    global _g38
    if _g38 is None:
        _g38 = B.load_genes(os.path.join(D, "Homo_sapiens.GRCh38.112.gtf.gz"))
    return _g38


def scan_finngen(path, pcut, want_dict=None, label="", keep_p=None):
    """扫描 FinnGen 文件。返回 (低于 pcut 的行, {rsid:(beta,se,p)} if keep_p)"""
    rows = []
    d = {} if keep_p is not None else None
    n = 0
    t0 = time.time()
    with gzip.open(path, "rt", errors="replace") as f:
        hdr = None
        for line in f:
            if hdr is None:
                hdr = line.rstrip("\n").split("\t")
                continue
            n += 1
            c = line.split("\t", 13)
            try:
                pv = float(c[6])
            except (ValueError, IndexError):
                continue
            if pv < pcut:
                rows.append(c[:13])
            if d is not None and pv < keep_p:
                try:
                    rs = c[4].split(",")[0].strip()
                    d[rs] = (float(c[8]), float(c[9]), pv, c[2], c[3])
                except Exception:
                    pass
    print(f"  {label}: 扫描 {n:,} 行 ({time.time()-t0:.0f}s); "
          f"p<{pcut:g} 有 {len(rows)} 位点; 缓存 {len(d) if d is not None else 0}")
    return rows, d


def loci(rows, win=500_000):
    ps = sorted(int(r[1]) for r in rows)
    out = []
    for x in ps:
        if not out or x - out[-1] > win:
            out.append(x)
    return out


def gprofiler_cam(genes, label):
    body = json.dumps({"query": genes, "organism": "hsapiens", "sources": ["KEGG", "GO:BP"],
                       "user_threshold": 1.0,
                       "significance_threshold_method": "gSCS"}).encode()
    req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
                                 data=body,
                                 headers={"Content-Type": "application/json",
                                          "Accept": "application/json",
                                          "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        res = json.load(r).get("result", [])
    sig = [x for x in res if x["p_value"] < 0.05]
    cam = [x for x in res if str(x.get("native", "")).endswith("04514")]
    cam_p = cam[0]["p_value"] if cam else None
    cam_n = cam[0]["intersection_size"] if cam else 0
    print(f"   富集: 显著条目 {len(sig)} (KEGG {len([x for x in sig if x['source']=='KEGG'])})")
    print(f"   >>> CAM: p={cam_p if cam_p is not None else '未返回'}"
          f"{'  (显著)' if cam_p and cam_p < 0.05 else '  (不显著)'}  交集={cam_n}")
    return cam_p, sig


def analyse(name, rows, pcut_label, genes_cache=True):
    print(f"\n{'='*80}\n### {name}: p<{pcut_label}\n{'='*80}")
    print(f"  位点 {len(rows)}; 按 500 kb 分窗约 {len(loci(rows))} 个区域")
    snps = [(r[0], int(r[1]), float(r[6])) for r in rows]
    res, un = B.map_snps(genes38(), snps)
    genes = sorted(res)
    print(f"  基因 {len(genes)} 个 (未映射 {un})")
    top = sorted(res.items(), key=lambda x: x[1])[:10]
    print("  最强基因:", ", ".join(f"{g}({v:.1e})" for g, v in top))
    tag = name.replace(" ", "_")
    with open(os.path.join(OUT, f"geneset_{tag}_{pcut_label}.txt"), "w", encoding="utf-8") as f:
        for g in genes:
            f.write(g + "\n")
    print(f"  关键基因: ", {g: ("在" if g in res else "—") for g in WATCH if g in res} or "无")
    ncam = len([g for g in CAM14 if g in res])
    print(f"  CAM 14 基因命中 {ncam}/14")
    if len(genes) >= 5:
        gprofiler_cam(genes, tag)
    return genes


if __name__ == "__main__":
    # ---- (A) R13 BELLPA ----
    print("=" * 80, "\n(A) FinnGen R13 G6_BELLPA (5,062 例)\n", "=" * 80)
    rows13, keep13 = scan_finngen(os.path.join(D, "finngen_R13_G6_BELLPA.gz"),
                                  1e-4, label="R13_BELLPA", keep_p=1e-3)
    for pcut, lab in [(5e-8, "p5e8"), (1e-5, "p1e5"), (1e-4, "p1e4")]:
        sub = [r for r in rows13 if float(r[6]) < pcut]
        if pcut > 5e-8:
            continue
        analyse("R13BELLPA", sub, lab)

    # ---- (B) R13 NERPLEX ----
    print("\n" + "=" * 80, "\n(B) FinnGen R13 G6_NERPLEX (68,559 例, 宽表型)\n", "=" * 80)
    rowsN, _ = scan_finngen(os.path.join(D, "finngen_R13_G6_NERPLEX.gz"),
                            1e-5, label="R13_NERPLEX")
    for pcut, lab in [(5e-8, "p5e8")]:
        sub = [r for r in rowsN if float(r[6]) < pcut]
        analyse("R13NERPLEX", sub, lab)

    # ---- (C) meta ----
    print("\n" + "=" * 80, "\n(C) 跨队列 meta: R13 BELLPA + GCST90018843\n", "=" * 80)
    print("  扫描旧 GWAS 缓存 p<1e-3 ...", flush=True)
    old = {}
    t0 = time.time()
    with gzip.open(OLD_VCF, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            keys = c[8].split(":"); vals = c[9].split(":")
            try:
                lp = float(vals[keys.index("LP")])
            except Exception:
                continue
            if lp > 3.0:
                old[c[2]] = (float(vals[keys.index("ES")]), float(vals[keys.index("SE")]),
                             c[3], c[4])
    print(f"  旧 GWAS 缓存 {len(old)} 个 (p<1e-3), 用时 {time.time()-t0:.0f}s")
    COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
    meta = []
    for rs, (b1, s1, p1, ref1, alt1) in keep13.items():
        o = old.get(rs)
        if not o:
            continue
        b2, s2, ref2, alt2 = o
        a1, o1 = alt1.upper(), ref1.upper()
        a2, o2 = alt2.upper(), ref2.upper()
        if a1 == a2 and o1 == o2:
            pass
        elif a1 == o2 and o1 == a2:
            b2 = -b2
        elif COMP.get(a1) == a2 and COMP.get(o1) == o2:
            pass
        elif COMP.get(a1) == o2 and COMP.get(o1) == a2:
            b2 = -b2
        else:
            continue
        w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
        bm = (b1 * w1 + b2 * w2) / (w1 + w2)
        sm = (1 / (w1 + w2)) ** 0.5
        z = abs(bm / sm)
        import math
        pm = math.erfc(z / math.sqrt(2))
        meta.append((rs, bm, sm, pm))
    print(f"  可比对并元分析: {len(meta)} 个位点")
    sig = [m for m in meta if m[3] < 5e-8]
    print(f"  meta p<5e-8: {len(sig)} 个位点")
    for m in sorted(sig, key=lambda x: x[3])[:20]:
        print(f"     {m[0]:<14} b={m[1]:+.4f} p={m[3]:.3g}")
    with open(os.path.join(OUT, "meta_bp_r13_plus_gcst90018843.csv"), "w", newline="",
              encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["rsid", "meta_beta", "meta_se", "meta_p"])
        for m in sorted(meta, key=lambda x: x[3]):
            wr.writerow(m)
    print(f"  -> {os.path.join(OUT, 'meta_bp_r13_plus_gcst90018843.csv')}")
