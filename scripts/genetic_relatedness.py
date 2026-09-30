"""run_autoimmune.py — 贝尔面瘫 vs 其他神经疾病 / 自身免疫病 的遗传边界

分析内容:
  1) 每个表型: 全基因组显著位点数、MHC 区(chr6:25-35Mb)最强信号、MHC 外最强信号
  2) 用统一的距离剪枝 SNP 面板算 贝尔面瘫 与各表型的 z 分数相关 (近似遗传相关)
  3) 各表型的独立 lead SNP 在贝尔面瘫中的效应方向一致性
  4) 贝尔面瘫 lead (rs9357446) 在各表型中的效应
  5) 一组公认自身免疫风险位点在贝尔面瘫中的效应
"""
import csv
import gzip
import math
import os
import statistics
import sys
import time

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D = _ROOT + r"\data"
W = _ROOT + r"\work"
OUT = _ROOT + r"\out"
MHC = ("6", 25_000_000, 35_000_000)

TRAITS = [
    ("G6_BELLPA", "贝尔面瘫", 5062),
    ("G6_NERPLEX", "神经/神经根/神经丛疾患", 68559),
    ("G6_CARPTU", "腕管综合征(卡压性)", 32951),
    ("G6_EPLEPSY", "癫痫(非自身免疫神经对照)", 16295),
    ("AUTOIMMUNE_NONTHYROID", "自身免疫病(非甲状腺)复合", 72412),
    ("E4_HYTHY_AI_STRICT", "自身免疫性甲减", 56098),
    ("T1D_WIDE", "1型糖尿病", 11197),
    ("M13_RHEUMA", "类风湿关节炎", 16775),
    ("L12_PSORIASIS", "银屑病", 13309),
    ("G6_MS", "多发性硬化", 3002),
    ("T2D", "2型糖尿病(代谢对照)", 89727),
]

# 公认自身免疫风险位点 (rsid, 风险等位基因, 疾病)
AI_VARIANTS = [
    ("rs2476601", "A", "PTPN22 R620W (RA/T1D/SLE/…多病共同)"),
    ("rs3087243", "A", "CTLA4 (T1D/RA/甲状腺)"),
    ("rs2104286", "A", "IL2RA (MS/T1D)"),
    ("rs6897932", "C", "IL7R (MS)"),
    ("rs763361", "T", "CD226 (MS/T1D)"),
    ("rs34536443", "G", "TYK2 P1104A (保护性, MS/RA)"),
    ("rs7574865", "T", "STAT4 (RA/SLE/Sjögren)"),
    ("rs10499194", "T", "TNFAIP3 (RA)"),
    ("rs11209026", "A", "IL23R R381Q (保护性, IBD/AS/银屑病)"),
    ("rs2066844", "T", "NOD2 (克罗恩)"),
    ("rs2241880", "G", "ATG16L1 T300A (克罗恩)"),
    ("rs2004640", "T", "IRF5 (SLE/Sjögren)"),
    ("rs2542151", "G", "PTPN2 (T1D/克罗恩)"),
    ("rs12708716", "A", "CLEC16A (T1D/MS)"),
    ("rs6822844", "G", "IL2/IL21 (T1D/乳糜泻/RA)"),
    ("rs917997", "A", "IL18RAP (乳糜泻)"),
    ("rs2187668", "A", "HLA-DQ2.5 tag (乳糜泻)"),
    ("rs9277534", "G", "HLA-DPB1 (MS/乳糜泻)"),
    ("rs3135388", "A", "HLA-DRB1*15:01 tag (MS)"),
    ("rs6457617", "T", "HLA-DQB1 (RA)"),
    ("rs10774671", "A", "OAS1 (T1D)"),
    ("rs17445836", "T", "IRF8 (MS)"),
    ("rs17810546", "A", "IL12A (乳糜泻/MS)"),
    ("rs10903129", "A", "TMEM57/… (T1D)"),
    ("rs2281808", "C", "SIRPG (T1D)"),
    ("rs6679677", "A", "PHTF1-PTPN22 (T1D/RA)"),
    ("rs4788084", "A", "IL27 (T1D/IBD)"),
    ("rs917997b", "A", "占位"),
]


def scan_trait(path, panel_rs=None):
    """返回 dict: gws, mhc_best, nonmhc_best, panel_z"""
    gws = []
    mhc_best = None
    non_best = None
    panel_z = {}
    n = 0
    with gzip.open(path, "rt", errors="replace") as f:
        hdr = next(f).rstrip("\n").split("\t")
        for line in f:
            n += 1
            c = line.split("\t", 13)
            if len(c) < 13:
                continue
            try:
                pv = float(c[6])
            except ValueError:
                continue
            chrom, pos = c[0], int(c[1])
            if pv < 5e-8:
                gws.append((chrom, pos, c[2], c[3], float(c[8]), float(c[9]), pv, c[5]))
            if chrom == "6" and MHC[1] <= pos <= MHC[2]:
                if mhc_best is None or pv < mhc_best[0]:
                    mhc_best = (pv, pos, c[2], c[3])
            else:
                if non_best is None or pv < non_best[0]:
                    non_best = (pv, pos, c[2], c[3])
            if panel_rs is not None and c[4]:
                rs = c[4].split(",")[0].strip()
                if rs in panel_rs:
                    try:
                        b, se = float(c[8]), float(c[9])
                        panel_z[rs] = b / se
                    except Exception:
                        pass
    return dict(n=n, gws=gws, mhc=mhc_best, nonmhc=non_best, panel=panel_z)


def loci(gws, win=1_000_000):
    ps = sorted((p for _, p, *_ in gws))
    out = []
    for x in ps:
        if not out or x - out[-1] > win:
            out.append(x)
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    # ---------- 先建统一面板 (来自贝尔面瘫) ----------
    print("建面板: 从贝尔面瘫 R13 取 MAF>0.05、间距>100kb 的独立 SNP ...", flush=True)
    panel = {}
    last_pos = {}
    t0 = time.time()
    with gzip.open(os.path.join(D, "finngen_R13_G6_BELLPA.gz"), "rt", errors="replace") as f:
        next(f)
        for line in f:
            c = line.split("\t", 13)
            if len(c) < 13:
                continue
            try:
                af = float(c[10]); b = float(c[8]); se = float(c[9])
            except ValueError:
                continue
            if not (0.05 <= af <= 0.95) or se <= 0:
                continue
            if not c[4]:
                continue
            rs = c[4].split(",")[0].strip()
            ch, pos = c[0], int(c[1])
            if pos - last_pos.get(ch, -10**9) < 100_000:
                continue
            last_pos[ch] = pos
            panel[rs] = b / se
    print(f"  面板 SNP 数: {len(panel)} ({time.time()-t0:.0f}s)", flush=True)

    avail = []
    for code, name, nc in TRAITS:
        p = os.path.join(D, f"finngen_R13_{code}.gz")
        if os.path.exists(p):
            avail.append((code, name, nc, p))
        else:
            print(f"  [缺] {code}")
    print(f"可分析表型: {len(avail)}")

    rows = []
    for code, name, nc, p in avail:
        t0 = time.time()
        r = scan_trait(p, set(panel.keys()))
        nl = len(loci(r["gws"]))
        mhc_sig = "是" if (r["mhc"] and r["mhc"][0] < 5e-8) else "否"
        common = [k for k in r["panel"] if k in panel]
        rho = float("nan")
        if len(common) > 100:
            xs = [panel[k] for k in common]
            ys = [r["panel"][k] for k in common]
            mx, my = statistics.mean(xs), statistics.mean(ys)
            num = sum((a-mx)*(b-my) for a, b in zip(xs, ys))
            den = (sum((a-mx)**2 for a in xs) * sum((b-my)**2 for b in ys)) ** 0.5
            rho = num/den if den else float("nan")
        print(f"\n{name} ({code}, {nc} 例)  用时 {time.time()-t0:.0f}s")
        print(f"   全基因组显著 SNP {len(r['gws'])}, 独立区域 ~{nl}")
        print(f"   MHC 区最强 p = {r['mhc'][0]:.3g} (pos {r['mhc'][1]})  显著={mhc_sig}"
              if r["mhc"] else "   MHC 区无数据")
        print(f"   MHC 外最强 p = {r['nonmhc'][0]:.3g} (pos {r['nonmhc'][1]})"
              if r["nonmhc"] else "")
        print(f"   与贝尔面瘫 z 相关 (n={len(common)}) = {rho:+.4f}")
        rows.append([name, code, nc, len(r["gws"]), nl,
                     r["mhc"][0] if r["mhc"] else "", mhc_sig,
                     r["nonmhc"][0] if r["nonmhc"] else "", len(common), rho])

    with open(os.path.join(OUT, "autoimmune_boundary_summary.csv"), "w", newline="",
              encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["trait", "phenocode", "n_cases", "n_gws_snp", "n_loci",
                     "mhc_best_p", "mhc_significant", "nonmhc_best_p", "n_panel", "z_corr_with_BP"])
        wr.writerows(rows)
    print(f"\n-> {os.path.join(OUT, 'autoimmune_boundary_summary.csv')}")

    # ---------- rs9357446 在各表型中的效应 ----------
    print("\n" + "=" * 90)
    print("贝尔面瘫 lead 位点 rs9357446 (chr6:44,479,861 G>A) 在各表型中的效应")
    print("=" * 90)
    print(f"{'表型':<28}{'beta':>10}{'se':>9}{'p':>12}{'af':>8}")
    for code, name, nc, p in avail:
        found = None
        with gzip.open(p, "rt", errors="replace") as f:
            next(f)
            for line in f:
                c = line.split("\t", 13)
                if len(c) > 4 and c[4] and "rs9357446" in c[4]:
                    found = c
                    break
        if found:
            print(f"{name[:26]:<28}{float(found[8]):>10.4f}{float(found[9]):>9.4f}"
                  f"{float(found[6]):>12.3g}{float(found[10]):>8.3f}")
        else:
            print(f"{name[:26]:<28}{'未找到':>10}")

    # ---------- 公认自身免疫位点在贝尔面瘫中的效应 ----------
    print("\n" + "=" * 90)
    print("公认自身免疫风险位点在贝尔面瘫(R13, 5,062 例)中的效应")
    print("=" * 90)
    p = os.path.join(D, "finngen_R13_G6_BELLPA.gz")
    want = {v[0]: v for v in AI_VARIANTS}
    hits = {}
    with gzip.open(p, "rt", errors="replace") as f:
        next(f)
        for line in f:
            c = line.split("\t", 13)
            if len(c) > 4 and c[4]:
                for rs in c[4].split(","):
                    rs = rs.strip()
                    if rs in want:
                        hits[rs] = c
    print(f"{'rsid':<14}{'疾病/基因':<40}{'风险等位':>9}{'BP_beta':>10}{'BP_p':>11}{'一致':>7}")
    n_conc = n_tot = 0
    for rs, ra, desc in AI_VARIANTS:
        c = hits.get(rs)
        if not c:
            print(f"{rs:<14}{desc[:38]:<40}{ra:>9}{'—':>10}{'—':>11}")
            continue
        try:
            b, se, pv = float(c[8]), float(c[9]), float(c[6])
        except Exception:
            continue
        alt = c[3].upper()
        # beta 相对 alt; 若风险等位=alt 则同向
        same = (alt == ra)
        conc = "同向" if (same and b > 0) or ((not same) and b < 0) else "反向"
        n_tot += 1
        if conc == "同向":
            n_conc += 1
        sig = " *" if pv < 0.05/len(AI_VARIANTS) else ""
        print(f"{rs:<14}{desc[:38]:<40}{ra:>9}{b:>10.4f}{pv:>11.3g}  {conc}{sig}")
    print(f"\n方向一致: {n_conc}/{n_tot}")
