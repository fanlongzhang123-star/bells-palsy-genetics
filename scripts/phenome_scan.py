"""phenome_scan.py -- 全表型组扫描 (对应正文 Table 4 / Fig. 4c)

对 Bell's palsy 的 lead 变异 rs9357446 (GRCh38 chr6:44,479,861 G>A) 查询
FinnGen release 13 pheweb 变异接口, 取该变异在全部 2,754 个表型上的关联,
按 P 值排序输出。

用法:  python phenome_scan.py [--variant 6-44479861-G-A] [--top 30]
输出:  out/phenome_scan_rs9357446.tsv  (phenocode, phenotype, pval, beta, sebeta, n_cases)
"""
import argparse
import csv
import json
import os
import ssl
import sys
import urllib.parse
import urllib.request

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

API = "https://r13.finngen.fi/api/variant/{}"
UA = {"User-Agent": "phenome-scan/1.0 (reproduction of the Bell's palsy analysis)"}


def fetch(variant):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(API.format(urllib.parse.quote(variant)), headers=UA)
    with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
        return json.load(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="6-44479861-G-A",
                    help="FinnGen 变异标识 chrom-pos-ref-alt (GRCh38)")
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "out"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    data = fetch(a.variant)
    allres = data.get("results", [])      # 接口返回数据库内全部表型
    rows = []
    for r in allres:
        pv = r.get("pval")
        if pv is None:
            continue
        rows.append({
            "phenocode": r.get("phenocode", ""),
            "phenotype": (r.get("phenostring") or "").replace("\t", " "),
            "pval": pv,
            "beta": r.get("beta"),
            "sebeta": r.get("sebeta"),
            "n_cases": r.get("num_cases") or r.get("n_case"),
            "mlogp": -__import__("math").log10(max(pv, 1e-320)),
        })
    rows.sort(key=lambda x: x["pval"])

    out = os.path.join(a.outdir, "phenome_scan_rs9357446.tsv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"变异 {a.variant}: 接口返回 {len(allres)} 个表型条目"
          f" (FinnGen R13 数据库表型总数), 其中 {len(rows)} 个对该变异有可用关联统计量")
    print(f"已写 {out}\n")
    print(f"{'排名':>4s}  {'表型':52s} {'P':>12s} {'beta':>10s}")
    for i, r in enumerate(rows[:a.top], 1):
        print(f"{i:>4d}  {r['phenotype'][:52]:52s} {r['pval']:12.2e} "
              f"{r['beta']:+10.4f}")


if __name__ == "__main__":
    main()
