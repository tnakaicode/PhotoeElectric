import scipy as sp
from scipy import pi, sin, cos, tan, arcsin, meshgrid, linspace, sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def sec(x):
    return 1 / cos(x)    # 関数 sec を定義


u = linspace(0, 0.4 * pi, 20)    # メッシュ作成の変数方位角 θ 周り 0 から 0.4π，20 等分
v = linspace(0, 2 * pi, 20)      # メッシュ作成の変数方位角 φ 周り 0 から 2π，20 等分

epz = 3     # z 方向の誘電率(正の値)
epx = 5     # x 方向の誘電率の絶対値(負の値)

uu, vv = meshgrid(u, v)  # メッシュの作成

x1 = sqrt(epz) * tan(uu) * cos(vv)    # 屈折率楕円体 x 方向(z が正の範囲)
y1 = sqrt(epz) * tan(uu) * sin(vv)    # 屈折率楕円体 y 方向(z が正の範囲)
z1 = sqrt(epx) * sec(uu)            # 屈折率楕円体 z 方向(z が正の範囲)

x2 = sqrt(epz) * tan(uu) * cos(vv)    # 屈折率楕円体 x 方向(z が負の範囲)
y2 = sqrt(epz) * tan(uu) * sin(vv)    # 屈折率楕円体 y 方向(z が負の範囲)
z2 = -sqrt(epx) * sec(uu)           # 屈折率楕円体 z 方向(z が負の範囲)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x1, y1, z1)     # ワイヤフレームのプロット
ax.plot_wireframe(x2, y2, z2)     # ワイヤフレームのプロット

ax.set_xlabel('X')  # x 方向ラベル
ax.set_ylabel('Y')  # y 方向ラベル
ax.set_zlabel('Z')  # z 方向ラベル

ax.set_xlim3d(-6, 6)  # x 方向プロット範囲
ax.set_ylim3d(-6, 6)  # y 方向プロット範囲
ax.set_zlim3d(-7, 7)  # z 方向プロット範囲

plt.show()  # グラフの表示
