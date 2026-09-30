"""gwasio.py — 通用 GWAS 汇总统计读取/提取工具

支持:
  - GWAS Catalog / FinnGen 风格 TSV(可 gzip)
  - OpenGWAS 统一 VCF (ES/SE/LP/AF/SS 在 FORMAT 字段)
功能:
  - 自动识别列名
  - 按 rsID 集合 和/或 染色体区域 提取
"""
import gzip
import os

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_CANON = {
    "chrom": ["#chrom", "chromosome", "chrom", "chr", "hg19chrc", "chr.x"],
    "pos": ["position", "pos", "base_pair_location", "bp", "hg19chrpos"],
    "ea": ["effect_allele", "alt", "allele1", "a1", "ea", "tested_allele"],
    "oa": ["other_allele", "ref", "allele0", "a2", "nea", "other_allele2"],
    "beta": ["beta", "effect", "b", "es", "beta1", "effect_size"],
    "se": ["standard_error", "se", "sebeta", "sebeta1", "stderr"],
    "p": ["p_value", "pval", "p-value", "p", "pvalue", "p_bolt_lmm_inf"],
    "neglogp": ["neg_log_10_p_value", "neg_log10_p", "neglog10p", "log_p", "logp"],
    "rsid": ["rs_id", "rsid", "rsids", "snp", "marker", "rs_number"],
    "vid": ["variant_id", "markername", "id"],
    "eaf": ["effect_allele_frequency", "eaf", "af_alt", "af", "freq", "a1freq",
            "effect_allele_freq", "alt_af", "frq_a1", "frq_a2", "freq_a1"],
    "n": ["n", "samplesize", "ss", "n_total", "totalsamplesize", "n_sample"],
    "info": ["info", "imputation_info", "info_score", "r2"],
    "nearest": ["nearest_genes", "gene", "nearest_gene"],
}


def _norm(s):
    return s.strip().lstrip("\ufeff").lower().replace(" ", "_")


def detect_columns(header):
    idx = {}
    hn = [_norm(h) for h in header]
    for key, cands in _CANON.items():
        for c in cands:
            if c in hn:
                idx[key] = hn.index(c)
                break
        if key not in idx:
            for i, h in enumerate(hn):
                if any(len(c) >= 3 and h.startswith(c) for c in cands):
                    idx[key] = i
                    break
    return idx


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "r", errors="replace")


def _split(line):
    """自动识别制表符 / 空白分隔"""
    if "\t" in line:
        return line.split("\t")
    return line.split()


def peek_header(path, nlines=3):
    """返回 (header_list, sample_lines)"""
    with _open(path) as f:
        lines = []
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            if header is None:
                header = _split(line.rstrip("\n"))
                continue
            lines.append(line.rstrip("\n"))
            if len(lines) >= nlines:
                break
    return header, lines


def scan(path, targets=None, region=None, keep_all_cols=False,
         max_rows=None, verbose=True, pos_targets=None):
    """
    targets: set of rsID (大小写敏感, 形如 rs123)
    region:  (chrom_str, start, end) 1-based 闭区间
    pos_targets: set of int, 按位置精确匹配(用于无 rsID 的数据集)
    返回 list[dict], 每个 dict 含 chrom,pos,ea,oa,beta,se,p,rsid,eaf,n,info
    """
    targets = targets or set()
    pos_targets = pos_targets or set()
    rows = []
    seen = 0
    header = None
    idx = None
    sep = "\t"
    with _open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if header is None:
                sep = "\t" if "\t" in line else None
                header = _split(line.rstrip("\n"))
                idx = detect_columns(header)
                if verbose:
                    print(f"    header[{len(header)}]: {header[:14]}  sep={'TAB' if sep else 'WHITESPACE'}")
                    print(f"    detected: { {k: header[v] for k, v in idx.items()} }")
                continue
            seen += 1
            if max_rows and seen > max_rows:
                break
            if sep is None:
                line = line.rstrip("\n")
            hit = False
            # 廉价预过滤: 只切到需要的列
            ri = idx.get("rsid")
            if region is not None:
                p = line.split(sep, 3) if sep else line.split(None, 3)
                try:
                    c = p[idx["chrom"]].strip()
                    pos = int(p[idx["pos"]])
                except Exception:
                    continue
                if c in (region[0], "chr" + region[0]) and region[1] <= pos <= region[2]:
                    hit = True
            if not hit and pos_targets:
                lim = idx["pos"] + 1
                p0 = line.split(sep, lim) if sep else line.split(None, lim)
                try:
                    if int(p0[idx["pos"]]) in pos_targets:
                        hit = True
                except Exception:
                    pass
            if not hit and targets and ri is not None:
                p = line.split(sep, ri + 1) if sep else line.split(None, ri + 1)
                if ri < len(p):
                    for r in p[ri].split(","):
                        if r.strip() in targets:
                            hit = True
                            break
            if not hit:
                continue
            p = _split(line)
            def g(key, cast=str):
                i = idx.get(key)
                if i is None or i >= len(p):
                    return None
                v = p[i].strip()
                if v in ("", ".", "NA", "nan"):
                    return None
                try:
                    return cast(v)
                except Exception:
                    return None
            pv = g("p", float)
            if pv is None:
                nl = g("neglogp", float)
                if nl is not None:
                    pv = 10.0 ** nl if nl < 0 else 10.0 ** (-nl)
            ea = g("ea")
            oa = g("oa")
            rows.append(dict(
                chrom=g("chrom"), pos=g("pos", int), rsid=g("rsid"), vid=g("vid"),
                ea=ea.upper() if ea else ea, oa=oa.upper() if oa else oa,
                beta=g("beta", float),
                se=g("se", float), p=pv, eaf=g("eaf", float),
                n=g("n", float), info=g("info", float), nearest=g("nearest"),
            ))
    if verbose:
        print(f"    scanned {seen:,} data rows -> {len(rows)} kept")
    return rows


def scan_vcf(path, targets):
    """解析 OpenGWAS 统一 VCF, 返回 ES(alt 等位基因效应)/SE/LP/AF/SS"""
    rows = []
    with _open(path) as f:
        fmt_keys = None
        for line in f:
            if line.startswith("##"):
                continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                continue
            if len(p) < 10:
                continue
            vid = p[2]
            if vid not in targets:
                continue
            keys = p[8].split(":")
            vals = p[9].split(":")
            d = dict(zip(keys, vals))
            es = d.get("ES")
            if es in (None, ".", ""):
                continue
            rows.append(dict(
                chrom=p[0], pos=int(p[1]), rsid=vid, ref=p[3], alt=p[4],
                beta=float(es),
                se=float(d["SE"]) if d.get("SE") not in (None, ".", "") else None,
                lp=float(d["LP"]) if d.get("LP") not in (None, ".", "") else None,
                eaf=float(d["AF"]) if d.get("AF") not in (None, ".", "") else None,
                n=float(d["SS"]) if d.get("SS") not in (None, ".", "") else None,
                nc=float(d["NC"]) if d.get("NC") not in (None, ".", "") else None,
            ))
    return rows


def write_tsv(rows, path, cols=None):
    cols = cols or ["chrom", "pos", "rsid", "ea", "oa", "beta", "se", "p", "eaf", "n", "info"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    print(f"    wrote {len(rows)} rows -> {path}")
