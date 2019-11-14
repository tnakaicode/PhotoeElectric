import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
import math
from scipy.special import spherical_jn, spherical_yn
from scipy.special import jv, jvp, hankel1, h1vp

import RI


def kjo(k):
    return math.factorial(k)


def perpen(k, j, r0, ep1, ep2, ep3):
    return ((ep2 - ep1) * (ep1 - ep3) * k * kjo(k + j)) / ((ep2 + ep1) * ((k + 1) * ep1 + k * ep3) * kjo(k) * kjo(j) * (2 * r0)**(k + j + 1))


def parallel(k, j, r0, ep1, ep2, ep3):
    return ((ep2 - ep1) * (ep1 - ep3) * k * kjo(k + j)) / ((ep2 + ep1) * ((k + 1) * ep1 + k * ep3) * kjo(k + 1) * kjo(j - 1) * (2 * r0)**(k + j + 1))


def uhen(ep1, ep2, ep3):
    return (ep1 - ep3) / (2 * ep1 + ep3)


WLmin = RI.WLmin
WLmax = RI.WLmax
WLperiod = RI.WLperiod
WLx = RI.WLx
NumWLx = RI.NumWLx

k0 = 2 * np.pi / WLx    # 真空中の波数

qq = 15          # 計算する多重極の次数
r = 50           # 粒子の半径
gap = 1          # ギャップ長
d = gap + r        # D パラメータ
r0 = d / r         # r0 パラメータ

alpha_A = np.zeros(NumWLx, dtype=complex)    # A 係数初期化
alpha_B = np.zeros(NumWLx, dtype=complex)    # B 係数初期化

ep1 = np.zeros(NumWLx, dtype=complex)
ep2 = np.zeros(NumWLx, dtype=complex)
ep3 = np.zeros(NumWLx, dtype=complex)
al = np.zeros([NumWLx, qq, qq], dtype=complex)
bl = np.zeros([NumWLx, qq, qq], dtype=complex)
fl = np.zeros([NumWLx, qq], dtype=complex)
Xal = np.zeros([NumWLx, qq], dtype=complex)
Xbl = np.zeros([NumWLx, qq], dtype=complex)
a11l = np.zeros([NumWLx, qq], dtype=complex)
b11l = np.zeros([NumWLx, qq], dtype=complex)


for i in range(NumWLx):
    ep1[i] = 1                 # 周辺媒質の誘電率
    ep2[i] = RI.epAu[i]        # 球の誘電率
    ep3[i] = RI.epAu[i]        # 基板の誘電率
    for k in range(qq):      # 連立方程式の作成 垂直方向(A 係数
        for j in range(qq):
            if k == j:
                al[i, k, j] = 1 + \
                    perpen(k + 1, j + 1, r0, ep1[i], ep2[i], ep3[i])
            else:
                al[i, k, j] = perpen(k + 1, j + 1, r0, ep1[i], ep2[i], ep3[i])

    for k in range(qq):      # 連立方程式の作成 垂直方向(B 係数)
        for j in range(qq):
            if k == j:
                bl[i, k, j] = 1 + \
                    parallel(k + 1, j + 1, r0, ep1[i], ep2[i], ep3[i])
            else:
                bl[i, k, j] = parallel(
                    k + 1, j + 1, r0, ep1[i], ep2[i], ep3[i])

    for k in range(qq):     # 連立方程式の作成 右辺
        if k == 0:
            fl[i, k] = uhen(ep1[i], ep2[i], ep3[i])
        else:
            fl[i, k] = 0

    Xal[i] = np.linalg.solve(al[i], fl[i])    # 連立方程式を解く(A 係数)
    Xbl[i] = np.linalg.solve(bl[i], fl[i])    # 連立方程式を解く(B 係数)

    alpha_A[i] = -4 * np.pi * r**3 * ep1[i] * Xal[i, 0]    # 分極率を求める(A 係数)
    alpha_B[i] = -4 * np.pi * r**3 * ep1[i] * Xbl[i, 0]    # 分極率を求める(B 係数)

Csca_A = k0**4 / (6 * np.pi) * abs(alpha_A)**2      # 散乱断面積を求める(A 係数
Csca_B = k0**4 / (6 * np.pi) * abs(alpha_B)**2      # 散乱断面積を求める(B 係数)
Cabs_A = k0 * np.imag(alpha_A)                  # 吸収断面積を求める(A 係数)
Cabs_B = k0 * np.imag(alpha_B)                  # 吸収断面積を求める(B 係数

Qsca_A = Csca_A / ((r**2) * np.pi)
Qabs_A = Cabs_A / ((r**2) * np.pi)
Qsca_B = Csca_B / ((r**2) * np.pi)
Qabs_B = Cabs_B / ((r**2) * np.pi)

plt.figure(figsize=(8, 6))
plt.plot(WLx, Qsca_A, label=r"$Q_{\rm sca}$", linewidth=3.0, color='black')
plt.plot(WLx, Qabs_A, label=r"$Q_{\rm abs}$", linewidth=3.0, color='gray')
plt.axis([400, 700, 0, 12])
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='upper left')
plt.tight_layout()

plt.figure(figsize=(8, 6))
plt.plot(WLx, Qsca_B, label=r"$Q_{\rm sca}$", linewidth=3.0, color='black')
plt.plot(WLx, Qabs_B, label=r"$Q_{\rm abs}$", linewidth=3.0, color='gray')
plt.axis([400, 700, 0, 4])
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='upper right')
plt.tight_layout()
plt.show()
