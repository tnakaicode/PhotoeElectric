import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

import RI

n1 = 1       # 周辺媒質の屈折率
n2 = RI.RIAg    # 粒子の屈折率
r = 25         # 粒子の半径
k = 2 * np.pi / RI.WLx    # 波数の配列
alpha = 4 * np.pi * (r**3) * (n1**2) * (n2**2 - n1**2) / \
    (n2**2 + 2 * n1**2)  # 分極率の計算
Csca = k**4 / (6 * np.pi) * abs(alpha)**2   # 散乱断面積
Cabs = k * np.imag(alpha)                   # 吸収断面積
Qsca = Csca / ((r**2) * np.pi)              # 散乱効率
Qabs = Cabs / ((r**2) * np.pi)              # 吸収効率

plt.figure(figsize=(8, 6))        # figure size
# real part of the refractive index of gold
plt.plot(RI.WLx, np.real(RI.RIAg), label="real", linewidth=3.0, color='black')
# real part of the refractive index of gold
plt.plot(RI.WLx, np.imag(RI.RIAg), label="imaginary", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)       # x-axis label
plt.ylabel("refractive", fontsize=22)          # y-axis label
plt.title("refractive of SIlver", fontsize=22)   # Title of the graph
plt.grid(True)                                 # Show Grid
plt.axis([300, 700, 0, 5])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='lower right')
plt.tight_layout()

plt.figure(figsize=(8, 6))        # figure size
# real part of the refractive index of gold
plt.plot(RI.WLx, np.real(RI.epAg), label="real", linewidth=3.0, color='black')
# real part of the refractive index of gold
plt.plot(RI.WLx, np.imag(RI.epAg), label="imaginary", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)       # x-axis label
plt.ylabel("dielective", fontsize=22)          # y-axis label
plt.title("dierective of SIlver", fontsize=22)   # Title of the graph
plt.grid(True)                                 # Show Grid
plt.axis([300, 700, -10, 5])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='lower right')
plt.tight_layout()

plt.figure(figsize=(8, 6))        # figure size
# real part of the refractive index of gold
plt.plot(RI.WLx, Csca / 10**4,
         label=r"$C_{{\rm sca}}$", linewidth=3.0, color='black')
# real part of the refractive index of gold
plt.plot(RI.WLx, Cabs / 10**4,
         label=r"$C_{{\rm abs}}$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)       # x-axis label
plt.ylabel("cross section(×10000 nm^2)", fontsize=22)          # y-axis label
plt.title("cross section", fontsize=22)   # Title of the graph
plt.grid(True)                                 # Show Grid
plt.axis([300, 700, 0, 5])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='lower right')
plt.tight_layout()

plt.figure(figsize=(8, 6))        # figure size
# real part of the refractive index of gold
plt.plot(RI.WLx, Qsca, label=r"$Q_{{\rm sca}}$", linewidth=3.0, color='black')
# real part of the refractive index of gold
plt.plot(RI.WLx, Qabs, label=r"$Q_{{\rm abs}}$", linewidth=3.0, color='gray')
plt.xlabel("wave (nm)", fontsize=22)       # x-axis label
plt.ylabel("scat/absorb rate", fontsize=22)          # y-axis label
plt.title("cross section", fontsize=22)   # Title of the graph
plt.title("scat/absorb rate (R=25 nm)", fontsize=22)   # Title of the graph
plt.grid(True)                                 # Show Grid
plt.axis([300, 700, 0, 20])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='lower right')
plt.tight_layout()
plt.show()
