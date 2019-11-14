import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
from scipy.special import spherical_jn, spherical_yn
from scipy.special import jv, jvp, hankel1, h1vp

import RI

n1 = RI.RIAu    # 回転楕円体の屈折率
n2 = 1       # 周辺媒質の屈折率
a = b = 10   # 非回転軸の長さ〔nm〕
c = 50       # 回転軸の長さ〔nm〕
ee = np.sqrt(1 - (a / c)**2)    # 離心率
lz = (1 - ee**2) / ee**2 * (1 / (2 * ee) *
                            np.log((1 + ee) / (1 - ee)) - 1)   # z 方向 反電場係数
lx = (1 - lz) / 2            # x 方向 反電場係数
k = 2 * np.pi / RI.WLx  # 真空中の波数

alphax = 4 * np.pi * a * b * c * (n2**2) * ((n1**2) - (n2**2)) / \
    (3 * ((n2**2) + lx * ((n1**2) - (n2**2))))    # x 方向 分極率
alphaz = 4 * np.pi * a * b * c * (n2**2) * ((n1**2) - (n2**2)) / \
    (3 * ((n2**2) + lz * ((n1**2) - (n2**2))))    # z 方向 分極率

Csca_x = k**4 / (6 * np.pi) * abs(alphax)**2     # x 方向 散乱断面積
Cabs_x = k * np.imag(alphax)      # x 方向 吸収断面積
Qsca_x = Csca_x / (a * a * np.pi)     # x 方向 散乱効率
Qabs_x = Cabs_x / (a * a * np.pi)     # x 方向 吸収効率

Csca_z = k**4 / (6 * np.pi) * abs(alphaz)**2      # z 方向 散乱断面積
Cabs_z = k * np.imag(alphaz)      # z 方向 吸収断面積
Qsca_z = Csca_z / (a * c * np.pi)     # z 方向 散乱効率
Qabs_z = Cabs_z / (a * c * np.pi)     # z 方向 吸収効率

plt.figure(figsize=(8, 6))
plt.plot(RI.WLx, Qsca_x, label=r"$Q_{{\rm sca},a}$", linewidth=3.0, color='black')
plt.plot(RI.WLx, Qabs_x, label=r"$Q_{{\rm abs},a}$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.title("scat / absorb rate", fontsize=22)
plt.grid(True)                                 # Show Grid
plt.axis([300, 1000, 0, 1])                       # Plot Range
plt.tick_params(labelsize=20)
plt.legend(fontsize=20, loc='upper right')
plt.tight_layout()

plt.figure(figsize=(8, 6))
plt.plot(RI.WLx, Qsca_z, label=r"$Q_{{\rm sca},c}$", linewidth=3.0, color='black')
plt.plot(RI.WLx, Qabs_z, label=r"$Q_{{\rm abs},c}$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)
plt.ylabel("scat / absorb rate", fontsize=22)
plt.title("scat / absorb rate", fontsize=22)
plt.grid(True)                                 # Show Grid
plt.axis([300, 1000, 0, 50])                     # Plot Range
plt.tick_params(labelsize=20)
plt.legend(fontsize=20, loc='upper left')
plt.tight_layout()
plt.show()
