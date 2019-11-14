import scipy as sp
from scipy import pi, sin, cos, tan, arcsin, meshgrid, linspace, sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def sec(x):
    return 1 / cos(x)     # 関数 sec を定義


u = linspace(0, 2 * pi, 20)    # メッシュ作成の変数方位角 θ 周り 0 から 2π，20 等分
v = linspace(0, 2 * pi, 20)    # メッシュ作成の変数方位角 φ 周り 0 から 2π，20 等分

epz = 3  # z 方向の誘電率の絶対値 (負の値)
epx = 5  # x 方向の誘電率 (正の値)
uu, vv = meshgrid(u, v)     # メッシュの作成

x = sqrt(epz) * sec(uu) * cos(vv)     # 屈折率楕円体の媒介変数表示 x 方向
y = sqrt(epz) * sec(uu) * sin(vv)     # 屈折率楕円体の媒介変数表示 y 方向
z = sqrt(epx) * tan(uu)             # 屈折率楕円体の媒介変数表示 z 方向

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, z)       # ワイヤフレームのプロット
ax.set_xlabel('X')     # x 方向ラベル
ax.set_ylabel('Y')     # y 方向ラベル
ax.set_zlabel('Z')     # z 方向ラベル

ax.set_xlim3d(-20, 20)     # x 方向プロット範囲
ax.set_ylim3d(-20, 20)     # y 方向プロット範囲
ax.set_zlim3d(-30, 30)     # z 方向プロット範囲

plt.show()  # グラフの表示
