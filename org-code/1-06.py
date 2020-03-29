import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import pi, sin, cos, tan, arcsin, exp, linspace, arange, sqrt, zeros, array, matrix, asmatrix
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, legend, grid, axis


def mMATs(n1z, n2z):
    return (1 / (2 * n1z)) * matrix([[n1z + n2z, n1z - n2z], [n1z - n2z, n1z + n2z]])
    # s-偏光 Mij 行列の定義


def mMATp(n1z, n2z, n1, n2):
    return (1 / (2 * n1 * n2 * n1z)) *\
        matrix([[n1**2 * n2z + n2**2 * n1z, n1**2 * n2z - n2**2 * n1z],
                [n1**2 * n2z - n2**2 * n1z, n1**2 * n2z + n2**2 * n1z]])
    # p-偏光 Mij 行列の定義


def matFAI(n1z, d1, k0):
    return matrix([[exp(1j * n1z * k0 * d1), 0], [0, exp(-1j * n1z * k0 * d1)]])
    # Φ 行列の定義


n1 = 1.0            # 媒質 1 の屈折率
n2 = 1.5            # 媒質 2 の屈折率
n3 = 1.0            # 媒質 3 の屈折率
ep1 = n1**2         # 媒質 1 の誘電率
ep2 = n2**2         # 媒質 2 の誘電率
ep3 = n3**2         # 媒質 3 の誘電率
d2 = 100            # 媒質 2 の厚さ d2〔nm〕
WL = 500            # 真空中の波長 WL〔nm〕
k0 = 2 * pi / WL        # 真空中の波数

t1start = 0         # 始めの角度
t1end = 89           # 終わりの角度
t1points = 90        # プロット数

t1Deg = linspace(t1start, t1end, t1points)  # 入射角 t1 の配列の生成
t1 = t1Deg / 180 * pi   # 入射角をラジアンに直す
s1 = sin(t1)         # sin(t1)
c1 = cos(t1)         # cos(t1)
s2 = n1 / n2 * s1        # sin(t1)
c2 = sqrt(1 - s2**2)   # cos(t2)
s3 = n1 / n3 * s1        # sin(t3)
c3 = sqrt(1 - s3**2)   # cos(t3)

n1z = n1 * c1         # n1z=k1z/k0
n2z = n2 * c2         # n2z=k1z/k0
n3z = n3 * c3         # n2z=k1z/k0

mMats21 = zeros((t1points, 2, 2), dtype=complex)        # s-偏光 M21 行列初期化
mMats32 = zeros((t1points, 2, 2), dtype=complex)        # s-偏光 M32 行列初期化
mMatp21 = zeros((t1points, 2, 2), dtype=complex)        # p-偏光 M21 行列初期化
mMatp32 = zeros((t1points, 2, 2), dtype=complex)        # p-偏光 M32 行列初期化
matFAI2 = zeros((t1points, 2, 2), dtype=complex)        # Φ2 行列初期化
matTs = zeros((t1points, 2, 2), dtype=complex)          # s-偏光伝搬行列 Ts 初期化
matTp = zeros((t1points, 2, 2), dtype=complex)          # p-偏光伝搬行列 Tp 初期化
rs = zeros((t1points), dtype=complex)                 # rs 初期化
ts = zeros((t1points), dtype=complex)                 # ts 初期化
rp = zeros((t1points), dtype=complex)                 # rp 初期化
tp = zeros((t1points), dtype=complex)                 # tp 初期化

for i in range(t1points):
    # M21 作成
    # M32 作成
    # M21 作成
    # M32 作成
    mMats21[i] = mMATs(n2z[i], n1z[i])
    mMats32[i] = mMATs(n3z[i], n2z[i])
    mMatp21[i] = mMATp(n2z[i], n1z[i], n2, n1)
    mMatp32[i] = mMATp(n3z[i], n2z[i], n3, n2)

    # Φ2 行列
    matFAI2[i] = matFAI(n2z[i], d2, k0)

    # s-偏光伝搬行列 Ts 作成
    # p-偏光伝搬行列 Tp 作成
    matTs[i] = mMats32[i]@matFAI2[i]@mMats21[i]
    matTp[i] = mMatp32[i]@matFAI2[i]@mMatp21[i]

    # Rs: s-偏光反射係数
    # Ts: s-偏光透過係数
    # Rp: p-偏光反射係数
    # Tp: p-偏光透過係数
    rs[i] = -matTs[i, 1, 0] / matTs[i, 1, 1]
    ts[i] = matTs[i, 0, 0] - matTs[i, 0, 1] * matTs[i, 1, 0] / matTs[i, 1, 1]
    rp[i] = -matTp[i, 1, 0] / matTp[i, 1, 1]
    tp[i] = matTp[i, 0, 0] - matTp[i, 0, 1] * matTp[i, 1, 0] / matTp[i, 1, 1]

    RsAbs = abs(rs)**2        # s-偏光反射率
    RpAbs = abs(rp)**2        # p-偏光反射率

    plt.plot(t1Deg, RpAbs, label="Rp")        # p-偏光反射率のプロット
    plt.plot(t1Deg, RsAbs, label="Rs")        # s-偏光反射率のプロット
    plt.xlabel(r"Inject (deg.)", fontsize=20)   # x 軸のラベル
    plt.ylabel(r"Reflect Rate", fontsize=20)              # y 軸のラベル
    plt.title("Reflect Rate", fontsize=20)          # グラフタイトル
    plt.grid(True)                                 # グリッドを表示
    plt.legend(fontsize=16, loc='upper left')            # 凡例の表示とフォントサイズ
    plt.tick_params(labelsize=18)  # 軸の目盛表示とフォントサイズの指定
plt.show()                                     # グラフを表示
