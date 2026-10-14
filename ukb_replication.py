#!/usr/bin/env python3
"""UK Biobank (Pan-UKB ICD-10 G51) look-up of the two FinnGen Bell's palsy credible sets.

Reproduces Table 2 and Supplementary Table S2 of the manuscript:

    Genetic architecture of Bell's palsy: a 6p21.1 locus replicated in UK Biobank and
    no classical HLA association in 5,062 Finnish cases

Data
----
Pan-UK Biobank, phenotype ``icd10-G51-both_sexes`` (ICD-10 G51, facial nerve disorders;
European ancestry; 1,046 cases and 419,485 controls):

    https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/icd10-G51-both_sexes.tsv.bgz

The tabix index sits alongside the file as ``....tsv.bgz.tbi``. The file is distributed in
GRCh37 with the header

    chr  pos  ref  alt  af_cases_EUR  af_controls_EUR  beta_EUR  se_EUR  neglog10_pval_EUR  low_confidence_EUR

Reference builds
----------------
Pan-UK Biobank is GRCh37 and FinnGen release 13 is GRCh38. Within the two windows analysed
the two builds differ by a constant 32,263 bp; we established this empirically by indexing
both datasets on allele pairs and tabulating the position difference between matched pairs
(1,176 and 1,025 exact matches at 32,263 bp, about 23 at any other offset). Change ``OFFSET``
to re-derive it for a different region.

Usage
-----
    python scripts/ukb_replication.py          # writes out/ukb_G51_credibleset_lookup.csv

Environment
-----------
    BELLSPALSY_ROOT   workspace root; defaults to the parent of this directory.
    Only the standard library plus the bundled ``tabix.py`` are required.
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tabix import TabixHTTP  # noqa: E402

ROOT = os.environ.get("BELLSPALSY_ROOT",
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
VARIANTS = os.path.join(ROOT, "results", "credibleset_variants_finngen_r13.tsv")
OUTDIR = os.path.join(ROOT, "out")
OUT = os.path.join(OUTDIR, "ukb_G51_credibleset_lookup.csv")

URL = ("https://pan-ukb-us-east-1.s3.amazonaws.com/"
       "sumstats_flat_files/icd10-G51-both_sexes.tsv.bgz")
OFFSET = 32_263                     # GRCh37 (Pan-UKB) -> GRCh38 (FinnGen) in bp
WINDOWS = ((44_350_000, 44_600_000),   # credible set 1
           (45_480_000, 45_720_000),   # credible set 2 (RUNX2)
           (45_890_000, 45_990_000))   # credible set 3


def fetch_ukb(tb):
    """Return {GRCh38 position: row} for the analysed windows."""
    ukb = {}
    for start, end in WINDOWS:
        rows = tb.fetch("6", start - OFFSET, end - OFFSET)
        for line in rows:
            f = line.rstrip("\n").split("\t")
            if not f or f[0].startswith("#"):
                continue
            ukb.setdefault(int(f[1]) + OFFSET, []).append({
                "ref": f[2].upper(), "alt": f[3].upper(),
                "af_cases": f[4], "af_controls": f[5],
                "beta": f[6], "se": f[7], "low_conf": f[9].strip().lower()})
    return ukb


def main():
    if not os.path.exists(VARIANTS):
        raise SystemExit(f"missing variant list: {VARIANTS}")
    os.makedirs(OUTDIR, exist_ok=True)

    ukb = fetch_ukb(TabixHTTP(URL))
    print(f"Pan-UK Biobank rows retrieved: {len(ukb)}")

    out_rows = []
    with open(VARIANTS, encoding="utf-8") as fh:
        for v in csv.DictReader(fh, delimiter="\t"):
            pos = int(v["variant_GRCh38"].split(":")[1])
            a1, a2 = v["variant_GRCh38"].split(":")[2].upper(), v["variant_GRCh38"].split(":")[3].upper()
            row = {"credible_set": v["credible_set"], "variant_GRCh38": v["variant_GRCh38"],
                   "rsID": v["rsID"], "FinnGen_beta": v["FinnGen_beta"],
                   "FinnGen_SE": v["FinnGen_SE"], "FinnGen_P": v["FinnGen_P"],
                   "UKB_af_cases": "", "UKB_af_controls": "", "UKB_beta": "", "UKB_SE": "",
                   "UKB_z": "", "UKB_P": "", "same_direction": "not in UKB"}
            for cand in ukb.get(pos, []):
                if (cand["ref"], cand["alt"]) == (a1, a2):
                    sign = 1.0
                elif (cand["ref"], cand["alt"]) == (a2, a1):
                    sign = -1.0
                else:
                    continue
                try:
                    beta, se = float(cand["beta"]), float(cand["se"])
                except ValueError:
                    row["same_direction"] = "not usable"
                    break
                if se <= 0 or cand["low_conf"] == "true":
                    row["UKB_af_cases"], row["UKB_af_controls"] = cand["af_cases"], cand["af_controls"]
                    row["same_direction"] = "not usable"
                    break
                z = sign * beta / se
                p = math.erfc(abs(z) / math.sqrt(2))
                row.update({"UKB_af_cases": cand["af_cases"], "UKB_af_controls": cand["af_controls"],
                            "UKB_beta": f"{sign * beta:.4f}", "UKB_SE": cand["se"],
                            "UKB_z": f"{z:+.2f}", "UKB_P": f"{p:.3g}",
                            "same_direction": "yes" if (float(v["FinnGen_beta"]) > 0) == (z > 0) else "no"})
                break
            out_rows.append(row)

    cols = list(out_rows[0].keys())
    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out_rows)
    print("wrote", OUT, f"({len(out_rows)} variants)")

    for cs in ("1", "2", "3"):
        sub = [r for r in out_rows if r["credible_set"] == cs]
        same = [r for r in sub if r["same_direction"] == "yes"]
        print(f"  credible set {cs}: {len(sub)} variants, {len(same)} with concordant direction")


if __name__ == "__main__":
    main()
