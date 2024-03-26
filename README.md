# wow
何も考えずに作ったのでこれからなにか考えます。

S-box以外は全部根拠不明です。
MDS行列を使ってますがCauchy行列でもいいようです。

# 設計指針
S-box:
$x^3+123 \pmod {257}$  

matrix:Vandermonde Matrix  

expand-key:
$\tau^{i+1}=\pi^{i-1}\tau^{i}\pi^{i-1}$

addkey:
$k[i]\oplus=k[\tau[i]]$
