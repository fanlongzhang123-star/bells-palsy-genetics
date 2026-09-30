"""tabix_local.py — 本地 BGZF/tabix 随机访问 (复用 TabixHTTP 的索引解析)
"""
import os

from tabix import TabixHTTP

_ROOT = os.environ.get("BELLSPALSY_ROOT",
                       os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TabixLocal(TabixHTTP):
    def __init__(self, path, tbi_path=None, timeout=90):
        self.path = path
        self.tbi_path = tbi_path or (path + ".tbi")
        super().__init__(path, tbi_url=self.tbi_path, timeout=timeout)

    def _http(self, url, start=None, end=None, retries=4):
        p = url
        with open(p, "rb") as f:
            if start is None:
                return f.read()
            f.seek(start)
            return f.read(end - start + 1)
