# push_v1.1.ps1 — 提交并推送 code_deposit，并打 v1.1 tag
#
# 用法：
#   powershell -ExecutionPolicy Bypass -File push_v1.1.ps1 -RepoUrl https://github.com/<用户名>/bells-palsy-genetics.git
#
# 本机注意：git 走 schannel 访问 GitHub 会失败，脚本内已统一加 -c http.sslBackend=openssl。
# 推完 tag 后，仍须到 GitHub 网页上 Create a new release（见 PUSH_v1.1.md 第 3 步）。

param(
    [Parameter(Mandatory = $true)][string]$RepoUrl,
    [string]$Root = "D:\dshwoek\code_deposit",
    [string]$Tag  = "v1.1",
    [string]$Message = "v1.1: add UK Biobank (ICD-10 G51) replication look-up; correct Zenodo metadata"
)

$ErrorActionPreference = "Stop"
Set-Location $Root
$git = { param($args_) & git -c http.sslBackend=openssl @args_ }

Write-Host "== repo: $Root"

if (-not (Test-Path (Join-Path $Root ".git"))) {
    Write-Host "== no git repo yet: init"
    & $git @("init")
    & $git @("branch", "-M", "main")
}

& $git @("add", ".")
$status = (& git -c http.sslBackend=openssl status --porcelain)
if ($status) {
    & $git @("commit", "-m", $Message)
} else {
    Write-Host "== nothing to commit"
}

$remotes = (& git -c http.sslBackend=openssl remote)
if ($remotes -contains "origin") {
    & $git @("remote", "set-url", "origin", $RepoUrl)
} else {
    & $git @("remote", "add", "origin", $RepoUrl)
}

& $git @("push", "-u", "origin", "main")

if ((& git -c http.sslBackend=openssl tag) -contains $Tag) {
    Write-Host "== tag $Tag already exists locally, skipping"
} else {
    & $git @("tag", "-a", $Tag, "-m", "v1.1: UK Biobank replication added")
}
& $git @("push", "origin", $Tag)

Write-Host ""
Write-Host "== done. Next (in the browser):"
Write-Host "   GitHub repo -> Releases -> Create a new release -> choose tag $Tag -> Publish release"
Write-Host "   (must be Published, not a draft, or Zenodo will not archive it)"
