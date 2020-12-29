import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
from scipy.special import spherical_jn, spherical_yn

from src.base import plot2d
from src.RI import RIAu, WLx

r3 = 25      # コアの半径
s = 1.5      # シェルの半径/コアの半径
r2 = r3 * s  # シェルの半径
n1 = 1       # 周辺媒質の屈折率
n2 = RIAu    # 微粒子コアの屈折率
n3 = 1.5     # シェルの屈折率
k = 2 * np.pi / WLx   # 波数の配列

delta = (n2**2) * (2 * (n1**2) * (1 + 2 * (s**3)) + (n3**2) *
                   (2 + s**3)) - 2 * ((n2**2)**2 + (n1**2) * (n3**2)) * (1 - s**3)
b11 = (s**3 / delta) * ((n2**2) * ((n1**2) * (1 + 2 * s**3) - (n3**2) *
                                   (2 + s**3)) + (2 * (n2**2)**2 - (n1**2) * (n3**2)) * (1 - s**3))

alpha = -4 * np.pi * r2**3 * (n1**2) * b11    # 分極率
Csca = k**4 / (6 * np.pi) * abs(alpha)**2     # 散乱断面積
Cabs = k * np.imag(alpha)          # 吸収断面積
Qsca = Csca / ((r2**2) * np.pi)    # 散乱効率
Qabs = Cabs / ((r2**2) * np.pi)   # 吸収効率

obj = plot2d(aspect="auto")
obj.axs.plot(WLx, Qsca, label=r"$Q_{{\rm sca}}$", color='black')
obj.axs.plot(WLx, Qabs, label=r"$Q_{{\rm abs}}$", color='gray')
obj.axs.set_xlabel("wave (nm)")
obj.axs.set_ylabel("scat / absorb rate")
obj.axs.axis([400, 1000, 0, 20])
obj.axs.tick_params()
obj.fig.tight_layout()
obj.SavePng()
