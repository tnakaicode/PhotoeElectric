
import scipy as sp
from scipy import pi, sin, cos, tan, meshgrid, arange
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

u = arange(0, 2 * pi, 0.1)  # メッシュ作成の変数方位角 φ 周り 0 から 2π，0.1 刻み
v = arange(0, 1 * pi, 0.1)  # メッシュ作成の変数方位角 θ 周り 0 から π，0.1 刻み

epz = 4              # z方向の誘電率
epx = 9              # x方向の誘電率

uu, vv = meshgrid(u, v)           # メッシュの作成

x = epz * cos(uu) * sin(vv)         # 誘電率楕円体の媒介変数表示  x 方向
y = epz * sin(uu) * sin(vv)         # 誘電率楕円体の媒介変数表示  y 方向
z = epx * cos(vv)                 # 誘電率楕円体の媒介変数表示  z 方向

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, z)           # ワイヤフレームのプロット
ax.set_xlabel('X')  # x 方向ラベル
ax.set_ylabel('Y')  # y 方向ラベル
ax.set_zlabel('Z')  # z 方向ラベル

ax.set_xlim3d(-10, 10)  # x 方向プロット範囲
ax.set_ylim3d(-10, 10)  # y 方向プロット範囲
ax.set_zlim3d(-10, 10)  # z 方向プロット範囲

plt.show()  # グラフの表示
