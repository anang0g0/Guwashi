# wow
何も考えずに作ったのでこれからなにか考えます。

S-box以外は全部根拠不明です。
MDS行列を使ってますがCauchy行列でもいいようです。

# 設計指針
S-box:非線形置換
$x^3+123 \pmod {257}$  

最大距離分離符号：Vandermonde Matrix over GF(256)

expand-key:合同
$\tau^{i+1}=\pi^{i-1}\tau^{i}\pi^{i-1}$

$k[i]\oplus=k[\tau[i]]$

addkey:加法
m[i]+=key[i]
