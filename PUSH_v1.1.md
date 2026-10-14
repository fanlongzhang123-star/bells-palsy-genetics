# 发布 v1.1：GitHub Release → Zenodo 自动归档（逐步操作）

> 目的：让 Zenodo 在**同一条记录**下生成 v1.1 新版本，从而拿到新版本 DOI（concept DOI 不变）。
> 本次 v1.1 的内容：新增 UK Biobank（Pan-UKB ICD-10 G51）复制查证脚本与随包位点清单，并用
> `.zenodo.json` 顺带修正 v1.0 的三处元数据（类型、许可、创建者）。

---

## 0. 前置检查（1 分钟）

- [ ] GitHub 上已有该仓库（v1.0 当时建的，例：`bells-palsy-genetics`）
- [ ] Zenodo 的 GitHub 归档开关仍为 **ON**：Zenodo → 头像 → Settings → GitHub →
      找到该仓库，开关在右边（若显示 `Sync now` 也可点）
- [ ] 仓库是 **Public**（私有仓库不会被自动归档）
- [ ] 仓库没有被删除或改名（改名会断掉 Zenodo 关联）

## 1. 本地准备（已为你做好）

`D:\dshwoek\code_deposit\` 内新增/改动：

| 文件 | 说明 |
|---|---|
| `scripts/ukb_replication.py` | **新增**：UKB G51 复制查证（tabix over HTTP 直读 Pan-UKB，GRCh37→GRCh38 偏移 32,263，按等位基因匹配，输出 71 行查证表） |
| `results/credibleset_variants_finngen_r13.tsv` | **新增**：随包的三可信集 71 个位点清单（来自 FinnGen R13 SuSiE 结果） |
| `.zenodo.json` | **新增**：Zenodo 元数据（software / MIT / Fanlong + ORCID / v1.1） |
| `README.md`、`MANIFEST.md`、`CITATION.cff` | 已更新：数据来源加 Pan-UKB、运行步骤加复制脚本、版本号 1.1 |

自检（已通过）：在该目录下运行 `python scripts/ukb_replication.py`，输出
`out/ukb_G51_credibleset_lookup.csv`，与稿件正文表 2 / 补充表 S2 **逐行一致（71 行 0 处差异）**。

## 2. 提交并推送

仓库里已有 v1.0 的历史，所以是**增量提交**：

```bash
cd D:\dshwoek\code_deposit
git add .
git commit -m "v1.1: add UK Biobank (ICD-10 G51) replication look-up; correct Zenodo metadata"
git push
```

**本机特别注意**：这台机器上 git 走 schannel 访问 GitHub 会报 `SEC_E_NO_CREDENTIALS`，
必须加 `-c http.sslBackend=openssl`，例如：

```bash
git -c http.sslBackend=openssl push
```

（也可以用 `push_v1.1.ps1`，脚本里已经带了这个开关。）

**不想用命令行**，两条替代路径：
1. **GitHub Desktop**：Add local repository → 选 `D:\dshwoek\code_deposit` → Commit → Push。
2. **网页上传**：仓库页 → `Add file` → `Upload files` → 拖入改动文件 → Commit（不需要 git）。

## 3. 发 Release（这一步才触发 Zenodo）

1. 打开仓库页 → 右侧 **Releases** → **Create a new release**
2. **Choose a tag**：填 `v1.1`（**必须是与 v1.0 不同的新 tag**；`git push` 之后这一步会自动创建 tag）
3. **Target**：`main`
4. **Release title**：`v1.1 - UK Biobank replication added`
5. 描述里写清差异（可粘贴）：
   > Adds scripts/ukb_replication.py and the credible-set variant list, reproducing the UK Biobank
   > (Pan-UKB ICD-10 G51) look-up reported in Table 2 and Supplementary Table S2. Zenodo metadata
   > (type, licence, creator ORCID) supplied via .zenodo.json.
6. 点 **Publish release**（⚠️ **不要存成 draft**：draft release 不会被归档）

## 4. 拿新 DOI

- 等 **10 秒到 2 分钟**（Zenodo 通过 webhook 抓取）
- Zenodo → 右上头像 → **Upload** → 找到原记录 → 页面上会出现 **Versions** 区块，里面是 `v1.1`
- 该页面会给出两个 DOI：
  - **Version DOI**：只指这一版（形如 `10.5281/zenodo.XXXXXXXX`）
  - **Concept DOI**：指"全部版本"（页面上的 `Cite all versions` 按钮给的就是它）——**论文请引这个**
- 顺便确认元数据已被 `.zenodo.json` 覆盖（Resource type = Software、License = MIT、Creator 带 ORCID）

## 5. 回来更新稿件

把新 DOI 发我，我一次性替换三处（英文主稿 Code availability、投稿信、中文镜像）并重建 docx。
若你自己改，三处原句见 `SUBMISSION_CHECKLIST_JMG.md`。

---

## 常见坑

| 现象 | 原因与处理 |
|---|---|
| Release 发了，Zenodo 没动 | 仓库的 Zenodo 开关没开，或发的是 draft；改正后**重新发新 tag**（如 v1.1.1） |
| 改了已发布的 release，Zenodo 没更新 | Zenodo 不跟踪 release 的后续修改，必须发新 release |
| 别名/改名/删库 | Zenodo 关联断裂，只能新建记录（DOI 会变） |
| `license` 字段被拒 | 网页上从下拉框选 MIT；`.zenodo.json` 里写 `"MIT"` 通常可用 |
| git push 报凭据错误 | 用 `git -c http.sslBackend=openssl push`；或改用 GitHub Desktop / 网页上传 |
