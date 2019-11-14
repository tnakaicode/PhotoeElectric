import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy import pi, sin, cos, tan, arcsin, exp, linspace, arange, sqrt, zeros, array, matrix, asmatrix


def mMATs(n1z, n2z):
    # s-偏光 Mij 行列の定義
    return (1 / (2 * n1z)) * matrix([[n1z + n2z, n1z - n2z], [n1z - n2z, n1z + n2z]])


def mMATp(n1z, n2z, n1, n2):
    return (1 / (2 * n1 * n2 * n1z)) *\
        matrix([[n1**2 * n2z + n2**2 * n1z, n1**2 * n2z - n2**2 * n1z],
                [n1**2 * n2z - n2**2 * n1z, n1**2 * n2z + n2**2 * n1z]])
    # p-偏光 Mij 行列の定義


def matFAI(n1z, d1, k0):
    # Φ 行列の定義
    return matrix([[exp(1j * n1z * k0 * d1), 0], [0, exp(-1j * n1z * k0 * d1)]])


n1 = 1.86           # 媒質 1(プリズム)の屈折率
n2 = sqrt(-10.8 + 1j * 1.47)    # 媒質 2(金)の屈折率
n3 = 1.5            # 媒質 3(誘電体薄膜)の屈折率
n4 = 1.33           # 媒質 4(水)の屈折率
ep1 = n1**2         # 媒質 1 の誘電率
ep2 = n2**2         # 媒質 2 の誘電率
ep3 = n3**2         # 媒質 3 の誘電率
ep4 = n4**2         # 媒質 4 の誘電率
d2 = 47             # 媒質 2(金)の厚さ d2〔nm〕
d3 = 10             # 媒質 2(誘電体)の厚さ d2〔nm〕
WL = 633            # 真空中の波長 WL〔nm〕
k0 = 2 * pi / WL        # 真空中の波数

t1start = 40        # 始めの角度
t1end = 70          # 終わりの角度
t1points = 300      # プロット数

t1DegOut = linspace(t1start, t1end, t1points)    # 入射角 t1 の配列の生成
t1 = 0.25 * pi + (1 / n1) * arcsin((t1DegOut - 45) / 180 * pi)  # 入射角をラジアンに直す
s1 = sin(t1)         # sin(t1)
c1 = cos(t1)         # cos(t1)
s2 = n1 / n2 * s1        # sin(t2)
c2 = sqrt(1 - s2**2)   # cos(t2)
s3 = n1 / n3 * s1        # sin(t3)
c3 = sqrt(1 - s3**2)   # cos(t3)
s4 = n1 / n4 * s1        # sin(t4)
c4 = sqrt(1 - s4**2)   # cos(t4)

n1z = n1 * c1         # n1z=k1z/k0
n2z = n2 * c2         # n2z=k1z/k0
n3z = n3 * c3         # n2z=k1z/k0
n4z = n4 * c4         # n2z=k1z/k0

matT0 = zeros((t1points, 2, 2), dtype=complex)          # 伝搬行列 T0 の初期化
matT1 = zeros((t1points, 2, 2), dtype=complex)          # 伝搬行列 T1 の初期化
r0 = zeros((t1points), dtype=complex)                 # 反射係数(誘電体なし)の初期化
r1 = zeros((t1points), dtype=complex)                 # 反射係数(誘電体あり)の初期化

for i in range(t1points):

    matT0[i] = mMATp(n4z[i], n2z[i], n4, n2)@matFAI(
        n2z[i], d2, k0)@mMATp(n2z[i], n1z[i], n2, n1)
    # p-偏光伝搬行列 T0
    matT1[i] = mMATp(n4z[i], n3z[i], n4, n3)@matFAI(n3z[i], d3, k0)@mMATp(
        n3z[i], n2z[i], n3, n2)@matFAI(n2z[i], d2, k0)@mMATp(n2z[i], n1z[i], n2, n1)
    # p-偏光伝搬行列 T1

    r0[i] = -matT0[i, 1, 0] / matT0[i, 1, 1]    # p-偏光反射係数(誘電体なし)
    r1[i] = -matT1[i, 1, 0] / matT1[i, 1, 1]    # p-偏光反射係数(誘電体あり)

R0Abs = abs(r0)**2        # p-偏光反射率(誘電体なし)
R1Abs = abs(r1)**2        # p-偏光反射率(誘電体あり)

plt.figure(figsize=(8, 6))
plt.plot(t1DegOut, R1Abs, label="R1", linewidth=3.0, color='gray')
plt.plot(t1DegOut, R0Abs, label="R0", linewidth=3.0, color='black')
plt.xlabel(r"$\theta_1$ (deg.)", fontsize=20)
plt.ylabel(r"Reflectivity", fontsize=20)
plt.title("Surface Plasmon Resonance", fontsize=20)
plt.grid(True)
plt.legend(fontsize=20, loc='lower right')
plt.tick_params(labelsize=20)
plt.tight_layout()
plt.show()
