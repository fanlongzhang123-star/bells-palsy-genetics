"""liftover.py — 极简 UCSC chain 文件 liftover (hg19 -> hg38)

hg19ToHg38.over.chain 中 tName = hg19 (target), qName = hg38 (query)。
因此要从 target(t) 映射到 query(q)。

用法:
    lf = LiftOver(chain_path)
    pos38 = lf.lift("6", 44446919)     # 1-based hg19 -> 1-based hg38, 失败返回 None
"""
import os
import gzip

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class LiftOver:
    def __init__(self, chain_path):
        # hg19 chrom -> list of (tStart, tEnd, qStart, qStrand, qSize)
        self.blocks = {}
        self._rev = {}
        self._load(chain_path)

    def _load(self, path):
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as f:
            cur = None
            for line in f:
                if line.startswith("chain"):
                    p = line.split()
                    tName, tStrand, tStart = p[2], p[4], int(p[5])
                    qSize, qStrand, qStart, qEnd = int(p[8]), p[9], int(p[10]), int(p[11])
                    assert tStrand == "+", "只支持 target 正链"
                    cur = dict(chrom=self._norm(tName), t=tStart, q=qStart,
                               qStrand=qStrand, qSize=qSize)
                    self.blocks.setdefault(cur["chrom"], [])
                elif line[0].isdigit():
                    parts = line.split()
                    size = int(parts[0])
                    dt = int(parts[1]) if len(parts) > 1 else 0
                    dq = int(parts[2]) if len(parts) > 2 else 0
                    if cur is None:
                        continue
                    self.blocks[cur["chrom"]].append(
                        (cur["t"], cur["t"] + size, cur["q"], cur["qStrand"], cur["qSize"]))
                    cur["t"] += size + dt
                    cur["q"] += size + dq
        for k in self.blocks:
            self.blocks[k].sort()

    @staticmethod
    def _norm(c):
        return c[3:] if c.startswith("chr") else c

    def lift(self, chrom, pos):
        """1-based hg19 -> 1-based hg38, 失败返回 None"""
        arr = self.blocks.get(self._norm(str(chrom)))
        if not arr:
            return None
        p0 = pos - 1                       # 0-based
        # 二分查找
        lo, hi, i = 0, len(arr) - 1, None
        while lo <= hi:
            mid = (lo + hi) // 2
            ta, tb = arr[mid][0], arr[mid][1]
            if p0 < ta:
                hi = mid - 1
            elif p0 >= tb:
                lo = mid + 1
            else:
                i = mid
                break
        if i is None:
            return None
        ta, tb, q, qs, qSize = arr[i]
        off = p0 - ta
        if qs == "+":
            return q + off + 1
        # query 负链: hg38 坐标 = qSize - (q + off) - 1  (0-based) -> +1
        return qSize - (q + off)

    def _rev_arr(self, chrom):
        key = ("rev", self._norm(str(chrom)))
        if key not in self._rev:
            arr = []
            for (ta, tb, q, qs, qSize) in self.blocks.get(self._norm(str(chrom)), []):
                size = tb - ta
                if qs == "+":
                    qa, qb = q, q + size          # 0-based hg38 区间
                    arr.append((qa, qb, ta, "+"))
                else:
                    # hg38 0-based 区间 = [qSize-(q+size), qSize-q)
                    qa, qb = qSize - (q + size), qSize - q
                    arr.append((qa, qb, ta, "-"))
            arr.sort()
            self._rev[key] = arr
        return self._rev[key]

    def lift_rev(self, chrom, pos38):
        """1-based hg38 -> 1-based hg19, 失败返回 None"""
        arr = self._rev_arr(chrom)
        if not arr:
            return None
        p0 = pos38 - 1
        lo, hi, i = 0, len(arr) - 1, None
        while lo <= hi:
            mid = (lo + hi) // 2
            qa, qb = arr[mid][0], arr[mid][1]
            if p0 < qa:
                hi = mid - 1
            elif p0 >= qb:
                lo = mid + 1
            else:
                i = mid
                break
        if i is None:
            return None
        qa, qb, base, qs = arr[i]
        off = p0 - qa
        size = qb - qa
        if qs == "+":
            return base + off + 1
        return base + (size - off - 1) + 1


if __name__ == "__main__":
    lf = LiftOver(_ROOT + r"\data\hg19ToHg38.over.chain.gz")
    print("chroms:", len(lf.blocks))
    checks = [("6", 44_446_919, 44_479_182), ("6", 44_446_674, 44_478_937),
              ("6", 44_447_004, 44_479_267), ("6", 44_447_598, 44_479_861),
              ("6", 44_877_054, 44_809_317), ("19", 10_300_000, 10_189_324)]
    for c, p19, want in checks:
        got = lf.lift(c, p19)
        print(f"  {c}:{p19} -> {got}  (期望 {want})  {'OK' if got == want else 'MISMATCH'}")
