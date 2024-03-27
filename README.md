# wow
何も考えずに作ったのでこれからなにか考えます。

S-box以外は全部根拠不明です。
MDS行列を使ってますがCauchy行列でもいいようです。

# 設計指針
S-box:非線形置換
$x^3+123 \pmod {257}$  

最大距離分離符号：Vandermonde Matrix over GF(256)

expand-key:2つの置換の共役を取り、鍵自身にGF(2)上の演算をして拡張鍵を生成する。
2つの置換が秘密鍵から生成される秘密置換である。
$\tau^{i+1}=\pi^{i-1}\tau^{i}\pi^{i-1}、key[i]\oplus=key[\tau[i]]$

addkey:加法
$m[i]+=key[i]$
