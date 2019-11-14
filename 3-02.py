import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
from scipy.special import spherical_jn, spherical_yn
from scipy.special import jv, jvp, hankel1, h1vp

import RI


def h1v(n, x):
    return hankel1(n, x)


def a(n, x, mA, mB):
    return (mB * jv(n, mB * x) * jvp(n, mA * x) - mA * jv(n, mA * x) * jvp(n, mB * x)) / \
        (mB * jv(n, mB * x) * h1vp(n, mA * x) -
         mA * h1v(n, mA * x) * jvp(n, mB * x))


def b(n, x, mA, mB):
    return (mA * jv(n, mB * x) * jvp(n, mA * x) - mB * jv(n, mA * x) * jvp(n, mB * x)) / \
        (mA * jv(n, mB * x) * h1vp(n, mA * x) -
         mB * h1v(n, mA * x) * jvp(n, mB * x))


def fn(n, phi):
    return (1 / k0) * (np.cos(n * phi) + 1j * np.sin(n * phi)) * (pow(-1j, n))


def matA_tm(n, m1, m2, m3, x1, x2):
    return np.array([[m1 * jv(n, m1 * x1), -m2 * jv(n, m2 * x1), m2 * h1v(n, m2 * x1), 0],
                     [m1**2 * jvp(n, m1 * x1), -m2**2 * jvp(n, m2 * x1), m2 **
                      2 * h1vp(n, m2 * x1), 0],
                     [0, m2 * jv(n, m2 * x2), -m2 *
                      h1v(n, m2 * x2), m3 * h1v(n, m3 * x2)],
                     [0, m2**2 * jvp(n, m2 * x2), -m2**2 * h1vp(n, m2 * x2), m3**2 * h1vp(n, m3 * x2)]])


def matF_tm(n, m3, x2):
    return np.array([[0],
                     [0],
                     [m3 * jv(n, m3 * x2)],
                     [m3**2 * jvp(n, m3 * x2)]])


def matX_tm(n, m1, m2, m3, x1, x2):
    return np.linalg.solve(matA_tm(n, m1, m2, m3, x1, x2), matF_tm(n, m3, x2))


def matA_te(n, m1, m2, m3, x1, x2):
    return np.array([[m1**2 * jv(n, m1 * x1), -m2**2 * jv(n, m2 * x1), m2**2 * h1v(n, m2 * x1), 0],
                     [m1 * jvp(n, m1 * x1), -m2 * jvp(n, m2 * x1), m2 *
                      h1vp(n, m2 * x1), 0],
                     [0, m2**2 * jv(n, m2 * x2), -m2**2 * h1v(n, m2 * x2),
                      m3**2 * h1v(n, m3 * x2)],
                     [0, m2 * jvp(n, m2 * x2), -m2 * h1vp(n, m2 * x2), m3 * h1vp(n, m3 * x2)]])


def matF_te(n, m3, x2):
    return np.array([[0],
                     [0],
                     [m3**2 * jv(n, m3 * x2)],
                     [m3 * jvp(n, m3 * x2)]])


def matX_te(n, m1, m2, m3, x1, x2):
    return np.linalg.solve(matA_te(n, m1, m2, m3, x1, x2), matF_te(n, m3, x2))


r1 = 50  # コアの半径
r2 = 55  # シェルの半径
qq = 5   # ベッセル関数の次数

Qsca_tm = np.zeros(RI.NumWLx, dtype=float)
Qext_tm = np.zeros(RI.NumWLx, dtype=float)
Qabs_tm = np.zeros(RI.NumWLx, dtype=float)
Qsca_te = np.zeros(RI.NumWLx, dtype=float)
Qext_te = np.zeros(RI.NumWLx, dtype=float)
Qabs_te = np.zeros(RI.NumWLx, dtype=float)

k0 = 2 * np.pi / RI.WLx  # 真空中の波数
m1 = 1.5            # シェルの屈折率
m2 = RI.RIAg           # 円柱コアの屈折率
m3 = 1.0            # 周辺媒質の屈折率
x1 = k0 * r1          # サイズ パラメータ(コア)
x2 = k0 * r2          # サイズ パラメータ(シェル)

for i in range(RI.NumWLx):
    for n in range(-qq, qq):
        Qsca_tm[i] = Qsca_tm[i] + \
            (2 / x2[i]) * abs(matX_tm(n, m1, m2[i], m3, x1[i], x2[i])[3, 0])**2
        Qext_tm[i] = Qext_tm[i] + \
            (2 / x2[i]) * np.real(matX_tm(n, m1, m2[i], m3, x1[i], x2[i])[3, 0])
        Qsca_te[i] = Qsca_te[i] + \
            (2 / x2[i]) * abs(matX_te(n, m1, m2[i], m3, x1[i], x2[i])[3, 0])**2
        Qext_te[i] = Qext_te[i] + \
            (2 / x2[i]) * np.real(matX_te(n, m1, m2[i], m3, x1[i], x2[i])[3, 0])

Qabs_tm = Qext_tm - Qsca_tm
Qabs_te = Qext_te - Qsca_te

plt.figure(figsize=(8, 6))        # figure size
# Scattering Cross-section
plt.plot(RI.WLx, Qsca_tm, label=r"$Q_{\rm sca}$", linewidth=3.0, color='black')
# Extinction Cross-section
plt.plot(RI.WLx, Qabs_tm, label=r"$Q_{\rm abs}$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.title("scat / absorb rate", fontsize=22)
plt.grid(True)                                 # Show Grid
plt.axis([300, 800, 0, 1])                   # Plot Range
# show legend and fontsize=16pt
plt.legend(fontsize=20, loc='lower right')
plt.tick_params(labelsize=20)  # Show Legend
plt.tight_layout()

plt.figure(figsize=(8, 6))        # figure size
# Scattering Cross-section
plt.plot(RI.WLx, Qsca_te, label=r"$Q_{\rm sca}$", linewidth=3.0, color='black')
# Extinction Cross-section
plt.plot(RI.WLx, Qabs_te, label=r"$Q_{\rm abs}$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.title("scat / absorb rate", fontsize=22)
plt.grid(True)                                 # Show Grid
plt.axis([300, 1000, 0, 10])                   # Plot Range
# show legend and fontsize=16pt
plt.legend(fontsize=20, loc='upper left')
plt.tick_params(labelsize=20)  # Show Legend
plt.tight_layout()
plt.show()
