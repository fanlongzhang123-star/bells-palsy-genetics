"""build_bp_geneset.py — 构建贝尔麻痹 GWAS 的"真正"基因集

流程:
  1. 解析 Ensembl GTF -> 蛋白编码基因区间
  2. 把 SNP (FinnGen GRCh38 / 旧 GWAS GRCh37) 映射到基因
     (SNP 落在基因内, 或距最近基因 <= 10 kb —— 与上一轮 MPGK 复现相同协议)
  3. 输出各阈值下的基因集
"""
import csv
import gzip
import os
import sys

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

W = _ROOT + r"\work"
D = _ROOT + r"\data"
OUT = _ROOT + r"\out"


def load_genes(gtf_path, chrom_style="ensembl"):
    """返回 {chrom: [(start, end, gene_name, gene_id)]}, 仅蛋白编码基因"""
    genes = {}
    n = 0
    with gzip.open(gtf_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            if 'gene_biotype "protein_coding"' not in c[8]:
                continue
            chrom = c[0]
            try:
                s, e = int(c[3]), int(c[4])
            except ValueError:
                continue
            name = None
            gid = None
            for kv in c[8].split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    name = kv.split(" ", 1)[1].strip('"')
                elif kv.startswith("gene_id "):
                    gid = kv.split(" ", 1)[1].strip('"')
            if not name:
                continue
            genes.setdefault(chrom, []).append((s, e, name, gid))
            n += 1
    for ch in genes:
        genes[ch].sort()
    print(f"  蛋白编码基因 {n} 个 ({len(genes)} 条染色体)", flush=True)
    return genes


def map_snps(genes, snps, window=10_000):
    """snps: list of (chrom, pos, pval)。返回 {gene: min_p}, 以及未映射数"""
    import bisect
    res = {}
    unmapped = 0
    for chrom, pos, pv in snps:
        ch = str(chrom)
        g = genes.get(ch) or genes.get("chr" + ch)
        if not g:
            unmapped += 1
            continue
        starts = [x[0] for x in g]
        i = bisect.bisect_right(starts, pos)
        cands = []
        for j in (i - 1, i):
            if 0 <= j < len(g):
                cands.append(g[j])
        best = None
        for s, e, name, gid in cands:
            if s <= pos <= e:
                best = (0, name)
                break
            d = min(abs(pos - s), abs(pos - e))
            if d <= window and (best is None or d < best[0]):
                best = (d, name)
        if best is None:
            unmapped += 1
        else:
            res[best[1]] = min(res.get(best[1], 1.0), pv)
    return res, unmapped


def read_snps(path, pcol, chrom_col="chrom", pos_col="pos"):
    out = []
    with open(path, encoding="utf-8") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            try:
                out.append((d[chrom_col], int(float(d[pos_col])), float(d[pcol])))
            except Exception:
                continue
    return out


if __name__ == "__main__":
    print("=== 1) 载入 GRCh38 注释 (FinnGen 用) ===")
    g38 = load_genes(os.path.join(D, "Homo_sapiens.GRCh38.112.gtf.gz"))
    print("=== 2) 载入 GRCh37 注释 (旧 GWAS 用) ===")
    g37 = load_genes(os.path.join(D, "Homo_sapiens.GRCh37.87.gtf.gz"))

    sets = {}
    for th, key in [("5e-8", "p5e8"), ("1e-5", "p1e5"), ("1e-4", "p1e4")]:
        p = os.path.join(W, f"finngen_bp_{key}.tsv")
        snps = read_snps(p, "pval", "chrom", "pos")
        res, un = map_snps(g38, snps)
        sets[f"FinnGen_{th}"] = res
        print(f"\nFinnGen R11 p<{th}: {len(snps)} SNP -> {len(res)} 基因 (未映射 {un})")
        top = sorted(res.items(), key=lambda x: x[1])[:15]
        print("   最强基因:", ", ".join(f"{g}(p={v:.1e})" for g, v in top))

    # 旧 GWAS: 用本地提取的全区域? 没有全量, 用已有 657 基因集作对照
    old = set(l.strip() for l in open(
        r"C:\Users\12137\Desktop\os数据\GCST90018843\bells_palsy_MPGK_analysis"
        r"\GO_KEGG\基因映射\genes_p1e-4.txt", encoding="utf-8") if l.strip())
    old5 = set(l.strip() for l in open(
        r"C:\Users\12137\Desktop\os数据\GCST90018843\bells_palsy_MPGK_analysis"
        r"\GO_KEGG\基因映射\genes_p1e-5.txt", encoding="utf-8") if l.strip())
    print(f"\n上一轮 (旧 GWAS, GRCh37, Ensembl 87): p<1e-4 -> {len(old)} 基因; "
          f"p<1e-5 -> {len(old5)} 基因")

    print("\n=== 3) 两套基因集的重叠 ===")
    for key in ["FinnGen_1e-4", "FinnGen_1e-5", "FinnGen_5e-8"]:
        g = set(sets[key])
        print(f"  {key}: {len(g)} 基因; 与上一轮 657 集重叠 {len(g & old)}; "
              f"与 98 集重叠 {len(g & old5)}")

    print("\n=== 4) 关键基因在两套基因集中的地位 ===")
    key_g = ["PTPRC", "ICAM1", "ITGB2", "ITGAL", "SELE", "VCAM1", "PECAM1",
             "CTLA4", "PTPRM", "CADM3", "CDH4", "CLDN16", "CNTN2", "CNTNAP2",
             "GLG1", "ITGA9", "NRXN1", "NRXN3", "PTPRD", "NECTIN3", "PVRL3"]
    print(f"  {'基因':<10}{'旧p1e-4':>10}{'旧p1e-5':>10}{'FG p1e-4':>10}{'FG p1e-5':>10}{'FG p5e-8':>10}")
    for g in key_g:
        print(f"  {g:<10}{('是' if g in old else '—'):>10}{('是' if g in old5 else '—'):>10}"
              f"{('是' if g in sets['FinnGen_1e-4'] else '—'):>10}"
              f"{('是' if g in sets['FinnGen_1e-5'] else '—'):>10}"
              f"{('是' if g in sets['FinnGen_5e-8'] else '—'):>10}")

    os.makedirs(OUT, exist_ok=True)
    for k, v in sets.items():
        with open(os.path.join(OUT, f"geneset_{k}.txt"), "w", encoding="utf-8") as f:
            for g in sorted(v):
                f.write(g + "\n")
    print(f"\n基因集已写出到 {OUT}/geneset_*.txt")
