import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import zeros, array


num = 1

xmin = -100    # 計算範囲の設定
xmax = 100
ymin = -100
ymax = 100
zmin = -100
zmax = 100

numx = xmax - xmin + 1  # x 方向の計算する座標の数
numy = ymax - ymin + 1  # y 方向の計算する座標の数
numz = zmax - zmin + 1  # z 方向の計算する座標の数
num = numx * numy * numz      # すべての計算する座標の数

p = np.zeros([numx, numy, numz], dtype=int)

x = xmin
y = ymin
z = zmin

iii = 0
xorigin = 0       # x 方向 形状の重心 初期化
yorigin = 0       # y 方向 形状の重心 初期化
zorigin = 0       # z 方向 形状の重心 初期化

for z in range(zmin, zmax):
    for y in range(ymin, ymax):
        for x in range(xmin, xmax):
            if (x / 10)**2 + (y / 25)**2 + (z / 10)**2 <= 1:      # 形状を構成する座標か，の判断
                p[x - xmin, y - ymin, z - zmin] = 1      # 形状を構成する座標の場合には p=1
                xorigin = xorigin + x       # 形状の重心を求めるために x 座標の和をとる
                yorigin = yorigin + y       # 形状の重心を求めるために y 座標の和をとる
                zorigin = zorigin + z       # 形状の重心を求めるために z 座標の和をとる
                iii += 1
            else:
                p[x - xmin, y - ymin, z - zmin] = 0  # 形状を構成する座標ではない場合には p=0

Xorigin = xorigin / iii      # 形状の重心の x 座標
Yorigin = yorigin / iii      # 形状の重心の y 座標
Zorigin = zorigin / iii      # 形状の重心の z 座標

xx = zeros(iii, dtype=int)
yy = zeros(iii, dtype=int)
zz = zeros(iii, dtype=int)


i = 0
for z in range(zmin, zmax):
    for y in range(ymin, ymax):
        for x in range(xmin, xmax):
            if p[x - xmin, y - ymin, z - zmin] == 1:
                xx[i] = x
                yy[i] = y
                zz[i] = z
                i += 1

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter3D(np.ravel(xx), np.ravel(yy), np.ravel(zz))

plt.show()
