import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
import math
from scipy.special import spherical_jn, spherical_yn
from scipy.special import jv, jvp, hankel1, h1vp

from RI import WLx, NumWLx
from RI import epAg, epAu, RIAu, RIAg

xmin = -100   # 計算範囲の設定
xmax = 100
ymin = -100
ymax = 100
zmin = -100
zmax = 100

numx = xmax - xmin + 1  # x 方向の計算する座標の数
numy = ymax - ymin + 1  # y 方向の計算する座標の数
numz = zmax - zmin + 1  # z 方向の計算する座標の数
num = numx * numy * numz  # すべての計算する座標の数

p = np.zeros([numx, numy, numz], dtype=int)  # フラグ p(x,y,z) の初期化

x = xmin   # x 方向 形状の重心 初期化
y = ymin   # y 方向 形状の重心 初期化
z = zmin   # z 方向 形状の重心 初期化

iii = 0
xorigin = 0
yorigin = 0
zorigin = 0

for z in range(zmin, zmax):
    for y in range(ymin, ymax):
        for x in range(xmin, xmax):
            if (x / 20)**2 + (y / 50)**2 + (z / 20)**2 <= 1:  # 形状を構成する座標か，の判断
                # 形状を構成する座標の場合には p=1
                p[x - xmin, y - ymin, z - zmin] = 1
                #  p の配列は 0 以上の整数なので，xmin だけずらしている
                xorigin = xorigin + x   # 形状の重心を求めるために x 座標の和をとる
                yorigin = yorigin + y   # 形状の重心を求めるために y 座標の和をとる
                zorigin = zorigin + z   # 形状の重心を求めるために z 座標の和をとる
                iii += 1
            else:
                p[x - xmin, y - ymin, z - zmin] = 0   # 形状を構成する座標ではない場合に は p=0

Xorigin = xorigin / iii  # 形状の重心の x 座標
Yorigin = yorigin / iii  # 形状の重心の y 座標
Zorigin = zorigin / iii  # 形状の重心の z 座標

l1 = "--- ddscat calc for FROM_FILE ---"

l3 = "1.000   0.000   0.000"   # a1 vecator
l4 = "1.000   1.000   0.000"   # a2 vecator
l5 = "1.      1.      1.   "   # d_x/d  d_y/d  d_z/d  (normally 1 1 1)
l7 = "J     JX      JY      JZ    ICOMPX    ICOMPY   ICOMPZ"

f = open("shape.dat", "w")   # "shape.dat"ファイルを書込みモードでオープン
f.write(l1 + "\n")  # 1 行目 書込み
f.write(str(iii) + "\n")  # 2 行目 書込み (双極子数)
f.write(l3 + "\n")  # 3 行目 書込み (双極子数)
f.write(l4 + "\n")  # 4 行目 書込み (双極子数)
f.write(l5 + "\n")  # 5 行目 書込み (双極子数)
f.write(str(Xorigin) + "   " + str(Yorigin) + "   " +
        str(Zorigin) + "\n")  # 6 行目 書込み (双極子数)
f.write(l7 + "\n")  # 7 行目 書込み (双極子数)

ii = 1                 # 座標データの書込み  x-xmin が配列中のアドレスで x が実際の座標に対応
for z in range(zmin, zmax):
    for y in range(ymin, ymax):
        for x in range(xmin, xmax):
            if p[x - xmin, y - ymin, z - zmin] == 1:
                f.write(str(ii) + "    " + str(x) + "    " + str(y) +
                        "    " + str(z) + "    1      1      1" + "\n")
                ii += 1
f.close()
