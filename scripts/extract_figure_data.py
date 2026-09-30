"""extract_fig_data.py — 为 6 张图抽取真实数据
输出:
  work/fig_manhattan.tsv     chrom, pos, mlogp   (全部 p<1e-3 + 其余 1/20 抽稀)
  work/fig_region.tsv        chr6:42.9-46.0Mb 全部变异 (pos, ref, alt, rsid, pval, beta, se, af)
  work/fig_gws.tsv           全基因组显著变异明细
  work/fig_gws_summary.txt   汇总
"""
import os
import gzip, io, sys, time

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT
SRC = D + r"\data\finngen_R13_G6_BELLPA.gz"

t0 = time.time()
n_tot = n_kept = n_gws = 0
with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f, \
     open(D + r"\work\fig_manhattan.tsv", "w", encoding="utf-8") as fm, \
     open(D + r"\work\fig_region.tsv", "w", encoding="utf-8") as fr, \
     open(D + r"\work\fig_gws.tsv", "w", encoding="utf-8") as fg:
    hdr = next(f)
    fm.write("chrom\tpos\tmlogp\n")
    fr.write("chrom\tpos\tref\talt\trsids\tpval\tbeta\tsebeta\taf_alt\n")
    fg.write("chrom\tpos\trsid\tpval\tbeta\tsebeta\taf_alt\n")
    for line in f:
        c = line.split("\t")
        n_tot += 1
        chrom, pos, ref, alt, rsid, gene, pval, mlogp = c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]
        p = float(pval)
        if p < 1e-3:
            fm.write(f"{chrom}\t{pos}\t{mlogp}\n"); n_kept += 1
        elif n_tot % 20 == 0:
            fm.write(f"{chrom}\t{pos}\t{mlogp}\n"); n_kept += 1
        if chrom == "6":
            ip = int(pos)
            if 42_900_000 <= ip <= 46_000_000:
                fr.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{rsid}\t{pval}\t{c[8]}\t{c[9]}\t{c[10]}\n")
        if p < 5e-8:
            n_gws += 1
            fg.write(f"{chrom}\t{pos}\t{rsid}\t{pval}\t{c[8]}\t{c[9]}\t{c[10]}\n")
        if n_tot % 2_000_000 == 0:
            print(f"  {n_tot:,} rows  {time.time()-t0:.0f}s", flush=True)

print(f"总变异 {n_tot:,}  曼哈顿保留 {n_kept:,}  全基因组显著 {n_gws}")
print(f"耗时 {time.time()-t0:.0f}s")
