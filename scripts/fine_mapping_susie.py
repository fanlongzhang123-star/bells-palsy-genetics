"""finngen_susie2.py — 用 FinnGen 官方 SuSiE (SISu 芬蘭 LD, n=3775) 做最終判斷
  欄位 (66): trait,region,v,rsid,chromosome,position,allele1,allele2,maf,beta,se,p,
             mean,sd,prob,cs,cs_specific_prob,low_purity,lead_r2, [99 版同構], alpha1..10,
             mean1..10, sd1..10, lbf_variable1..10
"""
import csv
import os

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

W = _ROOT + r"\work"
OUT = _ROOT + r"\out"
LO, HI = 44_000_000, 45_800_000
COL = ["trait", "region", "v", "rsid", "chromosome", "position", "allele1", "allele2",
       "maf", "beta", "se", "p", "mean", "sd", "prob", "cs", "cs_specific_prob",
       "low_purity", "lead_r2", "mean_99", "sd_99", "prob_99", "cs_99",
       "cs_specific_prob_99", "low_purity_99", "lead_r2_99",
       "alpha1", "alpha2", "alpha3", "alpha4", "alpha5", "alpha6", "alpha7", "alpha8",
       "alpha9", "alpha10",
       "mean1", "mean2", "mean3", "mean4", "mean5", "mean6", "mean7", "mean8", "mean9",
       "mean10",
       "sd1", "sd2", "sd3", "sd4", "sd5", "sd6", "sd7", "sd8", "sd9", "sd10",
       "lbf1", "lbf2", "lbf3", "lbf4", "lbf5", "lbf6", "lbf7", "lbf8", "lbf9", "lbf10"]

PH = [("G6_BELLPA", "贝尔面瘫"), ("H8_OTOSCLE", "耳硬化症")]
data = {}
for ph, name in PH:
    fp = os.path.join(W, f"finngen_R13_{ph}.SUSIE.snp.region.tsv")
    if not os.path.exists(fp):
        print(f"缺少 {fp}")
        continue
    d = {}
    for line in open(fp, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 20 or f[0] == "trait":
            continue
        r = dict(zip(COL, f))
        try:
            r["pos"] = int(float(r["position"]))
            r["prob"] = float(r["prob"])
            r["cs"] = int(float(r["cs"]))
            r["p"] = float(r["p"])
        except Exception:
            continue
        if not (LO <= r["pos"] <= HI):
            continue
        d[r["pos"]] = r
        d[r["rsid"]] = r
    data[ph] = d
    print(f"{name}: 區域位點 {sum(1 for k in d if isinstance(k, str) and k.startswith('chr6_'))}")

# ---------- 1) FinnGen 官方 credible set (chr6) ----------
print("\n" + "=" * 100)
print("1) FinnGen 官方 SuSiE credible set (芬蘭 SISu LD) —— chr6p21.1 區段")
print("=" * 100)
for ph, name in PH:
    fp = os.path.join(W, f"finngen_R13_{ph}.SUSIE.cred.summary.tsv")
    rows = list(csv.DictReader(open(fp, encoding="utf-8"), delimiter="\t"))
    chr6 = [r for r in rows if r["region"].split(":")[0] == "chr6"
            and 42_000_000 <= int(r["region"].split(":")[1].split("-")[0]) <= 46_000_000]
    print(f"\n--- {name} ({ph}) 在 chr6:42.98-45.98 Mb 的 credible sets ---")
    print(f"{'cs':>4}{'size':>6}{'log10BF':>10}{'avg_r2':>8}{'min_r2':>8}{'good':>7}"
          f"{'lead_v':>22}{'lead_p':>10}  lead_rsid")
    for r in sorted(chr6, key=lambda x: int(x["cs"])):
        print(f"{r['cs']:>4}{r['cs_size']:>6}{float(r['cs_log10bf']):>10.2f}"
              f"{float(r['cs_avg_r2']):>8.3f}{float(r['cs_min_r2']):>8.3f}{r['good_cs']:>7}"
              f"{r['v']:>22}{float(r['p']):>10.2e}  {r['rsid']}")
    allchr6 = [r for r in rows if r["region"].split(":")[0] == "chr6"]
    print(f"   (全 chr6 共 {len(allchr6)} 個 CS)")

# ---------- 2) SNP 級 PIP 比較 ----------
print("\n" + "=" * 100)
print("2) 貝爾面癱 vs 耳硬化症 的 PIP 與 CS 歸屬 (芬蘭 LD)")
print("=" * 100)
bp, ot = data["G6_BELLPA"], data["H8_OTOSCLE"]
pos_b = {k: v for k, v in bp.items() if isinstance(k, str)}
pos_o = {k: v for k, v in ot.items() if isinstance(k, str)}
common = [k for k in pos_b if k in pos_o]
print(f"共有位點 {len(common)}")

print("\n--- (a) 貝爾面癱 PIP 最高的 15 個位點, 及其在耳硬化症中的 PIP ---")
print(f"{'v (GRCh38)':<24}{'pos':>11}{'BP_p':>10}{'BP_PIP':>8}{'BP_CS':>7}"
      f"{'OTO_p':>10}{'OTO_PIP':>9}{'OTO_CS':>8}")
for k in sorted(common, key=lambda x: -pos_b[x]["prob"])[:15]:
    b, o = pos_b[k], pos_o[k]
    print(f"{b['v']:<24}{b['pos']:>11}{b['p']:>10.2e}{b['prob']:>8.3f}{b['cs']:>7}"
          f"{o['p']:>10.2e}{o['prob']:>9.3f}{o['cs']:>8}")

print("\n--- (b) 耳硬化症 PIP 最高的 15 個位點 (本區段), 及其在貝爾面癱中的 PIP ---")
print(f"{'v (GRCh38)':<24}{'pos':>11}{'OTO_p':>10}{'OTO_PIP':>9}{'OTO_CS':>8}"
      f"{'BP_p':>10}{'BP_PIP':>8}{'BP_CS':>7}")
for k in sorted(common, key=lambda x: -pos_o[x]["prob"])[:15]:
    b, o = pos_b[k], pos_o[k]
    print(f"{o['v']:<24}{o['pos']:>11}{o['p']:>10.2e}{o['prob']:>9.3f}{o['cs']:>8}"
          f"{b['p']:>10.2e}{b['prob']:>8.3f}{b['cs']:>7}")

print("\n--- (c) 貝爾面癱 credible set 成員 (PIP 為該 CS 的 alpha), 看耳硬化症 ---")
for k in common:
    b = pos_b[k]
    if b["cs"] in (1, 2) and b["prob"] > 0.005:
        o = pos_o[k]
        print(f"  BP_CS{b['cs']}  {b['v']:<24} pos={b['pos']} BP_PIP={b['prob']:.4f} "
              f"| OTO_p={o['p']:.2e} OTO_PIP={o['prob']:.4f} OTO_CS={o['cs']}")

print("\n--- (d) 耳硬化症在本區段是否有 credible set ---")
cs_o = sorted({pos_o[k]["cs"] for k in common if pos_o[k]["cs"] > 0})
print(f"  耳硬化症在 44.0-45.8 Mb 出現的 CS 編號: {cs_o if cs_o else '無'}")
for c in cs_o:
    mem = [pos_o[k] for k in common if pos_o[k]["cs"] == c]
    mem.sort(key=lambda r: r["pos"])
    print(f"  OTO CS{c}: {len(mem)} SNP, GRCh38 {mem[0]['pos']}-{mem[-1]['pos']}")
    for r in sorted(mem, key=lambda x: -x["prob"])[:5]:
        bb = pos_b.get(r["rsid"])
        print(f"      {r['v']:<24} OTO_PIP={r['prob']:.4f} p={r['p']:.2e} "
              f"| BP_PIP={(bb['prob'] if bb else float('nan')):.4f} BP_p={(bb['p'] if bb else float('nan')):.2e}")

print("\n--- (e) 關鍵位點 ---")
for rs in ["rs7770012", "rs2281690", "rs7770034", "rs9357446", "rs9472521",
           "rs13192457", "rs9463063"]:
    hit = None
    for k in common:
        if rs in k:
            hit = k
            break
    if hit:
        b, o = pos_b[hit], pos_o[hit]
        print(f"  {rs:<14} {b['v']:<24} BP_PIP={b['prob']:.4f} CS={b['cs']} "
              f"| OTO_PIP={o['prob']:.4f} CS={o['cs']} OTO_p={o['p']:.2e}")
    else:
        print(f"  {rs:<14} (區域檔案中無此 rsID; FinnGen 用 chr6_pos_ref_alt 命名)")

# ---------- 3) 手動定位 rs7770012 / rs9472521 的位置 ----------
print("\n--- (f) 直接按位置查找 ---")
for lab, p38 in [("貝爾面癱主 locus lead (rs7770012)", 44_479_182),
                 ("貝爾面癱第二 locus lead (rs9472521/RUNX2)", 45_602_126),
                 ("已發表耳硬化症 chr6 lead (rs13192457)", 44_887_654),
                 ("耳硬化症大塊核心 (rs9463063)", 45_063_426)]:
    bk = [k for k in common if pos_b[k]["pos"] == p38]
    if bk:
        b, o = pos_b[bk[0]], pos_o[bk[0]]
        print(f"  {lab}\n     {b['v']}  BP: p={b['p']:.2e} PIP={b['prob']:.4f} CS={b['cs']}"
              f"   ||  OTO: p={o['p']:.2e} PIP={o['prob']:.4f} CS={o['cs']}")
    else:
        print(f"  {lab}: 位置 {p38} 不在檔案中")
