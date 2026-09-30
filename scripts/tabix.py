"""tabix.py — 极简 tabix over HTTP 读取器（BGZF 随机访问）

用法:
    from tabix import TabixHTTP
    tb = TabixHTTP(vcf_url)
    lines = tb.fetch("1", 198_200_000, 199_400_000)
"""
import os
import gzip
import io
import struct
import time
import urllib.request
import zlib

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

UA = {"User-Agent": "Mozilla/5.0"}


def bgzf_decompress(data: bytes) -> bytes:
    """解压由多个 BGZF(gzip) 成员拼接的数据 (高效实现)"""
    try:
        return gzip.decompress(data)
    except Exception:
        pass
    # 逐成员解压, 用 GzipFile + BytesIO (避免 O(n^2) 切片)
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as f:
            return f.read()
    except Exception:
        pass
    # 最后兜底: 增量喂给 decompressobj
    out = []
    pos = 0
    n = len(data)
    CH = 1 << 20
    while pos < n:
        d = zlib.decompressobj(31)
        res = []
        while pos < n:
            buf = data[pos:pos + CH]
            pos += len(buf)
            try:
                res.append(d.decompress(buf))
            except zlib.error:
                break
            if d.eof:
                extra = d.unused_data
                pos -= len(extra)
                break
        out.append(b"".join(res))
        if not d.eof:
            break
    return b"".join(out)


def reg2bins(beg, end):
    """tabix 标准 binning: [beg, end) 0-based"""
    if end <= beg:
        end = beg + 1
    end -= 1
    bins = [0]
    bins += list(range(1 + (beg >> 26), 1 + (end >> 26) + 1))
    bins += list(range(9 + (beg >> 23), 9 + (end >> 23) + 1))
    bins += list(range(73 + (beg >> 20), 73 + (end >> 20) + 1))
    bins += list(range(585 + (beg >> 17), 585 + (end >> 17) + 1))
    bins += list(range(4681 + (beg >> 14), 4681 + (end >> 14) + 1))
    return bins


class TabixHTTP:
    def __init__(self, url, tbi_url=None, timeout=90):
        self.url = url
        self.tbi_url = tbi_url or (url + ".tbi")
        self.timeout = timeout
        self.refs = []
        self.linear = []
        self.names = []
        self.header = {}
        self._load_index()

    def _http(self, url, start=None, end=None, retries=4):
        h = dict(UA)
        want = None
        if start is not None:
            h["Range"] = f"bytes={start}-{end}"
            want = end - start + 1
        last = None
        for a in range(retries):
            try:
                req = urllib.request.Request(url, headers=h)
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    data = r.read()
                if want is not None and len(data) < want:
                    # 服务端给了更少字节(可能是文件尾), 只要不是 0 就接受
                    if len(data) == 0:
                        raise IOError("空响应")
                return data
            except Exception as e:
                last = e
                time.sleep(2 + 3 * a)
        raise last

    def _load_index(self):
        raw = bgzf_decompress(self._http(self.tbi_url))
        if raw[:4] != b"TBI\x01":
            raise ValueError("不是有效的 tabix 索引")
        n_ref, = struct.unpack("<i", raw[4:8])
        # TBI 头部: magic, n_ref, format, col_seq, col_beg, col_end, meta, skip, l_nm, names
        fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<7i", raw[8:36])
        self.header = dict(format=fmt, col_seq=col_seq, col_beg=col_beg,
                           col_end=col_end, meta=meta, skip=skip)
        off = 36 + l_nm
        self.names = raw[36:36 + l_nm].split(b"\x00")
        for _ in range(n_ref):
            n_bin, = struct.unpack("<i", raw[off:off + 4]); off += 4
            bins = {}
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack("<Ii", raw[off:off + 8]); off += 8
                chunks = []
                for _ in range(n_chunk):
                    beg, end = struct.unpack("<QQ", raw[off:off + 16]); off += 16
                    chunks.append((beg, end))
                bins[bin_id] = chunks
            n_intv, = struct.unpack("<i", raw[off:off + 4]); off += 4
            self.linear.append(list(struct.unpack(f"<{n_intv}Q", raw[off:off + 8 * n_intv])))
            off += 8 * n_intv
            self.refs.append(bins)

    def _ref_index(self, chrom):
        c = str(chrom)
        for cand in (c, c.replace("chr", ""), "chr" + c):
            if cand.isdigit():
                i = int(cand)
                if 1 <= i <= len(self.refs):
                    return i - 1
        return 0

    def fetch(self, chrom, start, end, max_bytes=80 << 20):
        """start/end 为 1-based 闭区间; 返回该区域的原始文本行

        修正版: 逐 chunk 合并解压 (旧版把 bin 0 的整文件 chunk 并进来, 会读错数据)
        """
        ridx = self._ref_index(chrom)
        bins = self.refs[ridx]
        chunks = []
        for b in reg2bins(start - 1, end):
            if b == 0:
                continue                      # bin 0 覆盖整文件, 单独处理
            chunks.extend(bins.get(b, []))
        if not chunks:
            chunks = list(bins.get(0, []))
        if not chunks:
            return []
        # 合并重叠/相邻 chunk
        chunks.sort()
        merged = []
        for c in chunks:
            if merged and c[0] <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], c[1]))
            else:
                merged.append((c[0], c[1]))
        seq_c = max(0, self.header.get("col_seq", 1) - 1)
        beg_c = max(0, self.header.get("col_beg", 2) - 1)
        lim = max(seq_c, beg_c) + 1
        out = []
        total = 0
        for (lo, hi) in merged:
            b_start = lo >> 16
            b_end = (hi >> 16) + 65536
            if b_end - b_start > max_bytes:
                b_end = b_start + max_bytes
            data = self._http(self.url, b_start, b_end)
            total += len(data)
            txt = bgzf_decompress(data).decode("utf-8", "replace")
            for line in txt.split("\n"):
                if not line or line[0] == "#":
                    continue
                f = line.split("\t", lim)
                if len(f) <= beg_c:
                    continue
                try:
                    pos = int(f[beg_c])
                except ValueError:
                    continue
                if start <= pos <= end:
                    out.append(line)
            if total > max_bytes * 4:
                break
        return out
