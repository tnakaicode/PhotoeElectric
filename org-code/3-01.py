import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
from scipy.special import spherical_jn, spherical_yn
from scipy.special import jv, jvp, hankel1, h1vp

import RI


def h1v(n, x):               # ハンケル関数(短く定義)
    return hankel1(n, x)


def a(n, x, mA, mB):           # 散乱係数 an
    return (mB * jv(n, mB * x) * jvp(n, mA * x) - mA * jv(n, mA * x) * jvp(n, mB * x)) / \
        (mB * jv(n, mB * x) * h1vp(n, mA * x) -
         mA * h1v(n, mA * x) * jvp(n, mB * x))


def b(n, x, mA, mB):           # 散乱係数 bn
    return (mA * jv(n, mB * x) * jvp(n, mA * x) - mB * jv(n, mA * x) * jvp(n, mB * x)) / \
        (mA * jv(n, mB * x) * h1vp(n, mA * x) -
         mB * h1v(n, mA * x) * jvp(n, mB * x))


def fn(n, phi):              # 係数 fn
    return (1 / k0) * (np.cos(n * phi) + 1j * np.sin(n * phi)) * (pow(-1j, n))


rr = 25      # 円柱の半径
qq = 20      # ベッセル関数の次数

k0 = 2 * np.pi / RI.WLx  # 真空中の波数
m1 = RI.RIAg                    # 円柱の屈折率
m2 = 1.0                     # 周辺媒質の屈折率
x = k0 * m2 * rr  # サイズ パラメータ

Qsca_tm = np.zeros(RI.NumWLx, dtype=float)
Qext_tm = np.zeros(RI.NumWLx, dtype=float)
Qsca_te = np.zeros(RI.NumWLx, dtype=float)
Qext_te = np.zeros(RI.NumWLx, dtype=float)

for n in range(-qq, qq):
    Qsca_tm = Qsca_tm + (2 / x) * abs(b(n, x, m2, m1))**2   # TM 偏光散乱効率
    Qext_tm = Qext_tm + (2 / x) * np.real(b(n, x, m2, m1))     # TM 偏光消光効率
    Qsca_te = Qsca_te + (2 / x) * abs(a(n, x, m2, m1))**2   # TE 偏光散乱効率
    Qext_te = Qext_te + (2 / x) * np.real(a(n, x, m2, m1))     # TE 偏光消光効率

Qabs_tm = Qext_tm - Qsca_tm
Qabs_te = Qext_te - Qsca_te

plt.figure(figsize=(8, 6))
plt.plot(RI.WLx, Qsca_tm,
         label=r"$Q_{\rm sca}(\rm TM)$", linewidth=3.0, color='black')
plt.plot(RI.WLx, Qsca_te,
         label=r"$Q_{\rm sca}(\rm TE)$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.title("scat / absorb rate", fontsize=22)
plt.grid(True)
plt.axis([300, 800, 0, 5])
plt.legend(fontsize=20, loc='upper right')
plt.tick_params(labelsize=20)
plt.tight_layout()
plt.show()
