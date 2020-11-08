import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import matplotlib as mpl
import scipy.constants as cnt

from src.base import pol_plot
from material import mat_refract

freq = 170.0 * 10**(3 * 3)  # Hz
wave = cnt.c / freq  # m
knum = 2 * np.pi / wave
print(wave, knum)

# Refractive index
n1 = 1.0
n2 = mat_refract["Diamond"]

w1 = cnt.c / (n1 * freq)
w2 = cnt.c / (n2 * freq)
k1 = 2 * np.pi / w1
k2 = 2 * np.pi / w2
print(w1, k1)
print(w2, k2)

# Brewstar angle
# ref-rate = 0
# trs-rate = 1
theta_b = np.arctan(n2 / n1)

t1Deg = np.linspace(0, 90, 90)  # 入射角 t1 の配列の生成
t1 = t1Deg / 180 * np.pi          # 入射角をラジアンに直す
t2 = np.arcsin((n1 / n2) * np.sin(t1))   # 屈折角 t2 を求める

s1 = np.sin(t1)
c1 = np.cos(t1)
s2 = n1 / n2 * s1
c2 = np.sqrt(1 - s2**2)
n1z = n1 * c1
n2z = n2 * c2

# p-pol: TM: parallel to incident surface
# s-pol: TE: vertical to incident aurface
# tp: p-pol trs-coef: p-polarize transmission coefficient
# rp: p-pol ref-coef: p-polarize reflection coefficient
# ts: s-pol trs-coef: s-polarize transmission coefficient
# rs: s-pol ref-coef: s-polarize reflection coefficient
tp = 2 * (n1 * c1) / (n2 * c1 + n1 * np.cos(t2))
rp = 1 * (n2 * c1 - n1 * np.cos(t2)) / (n2 * c1 + n1 * np.cos(t2))
ts = 2 * (n1 * c1) / (n1 * c1 + n2 * np.cos(t2))
rs = 1 * (n1 * c1 - n2 * np.cos(t2)) / (n1 * c1 + n2 * np.cos(t2))

# Tp: p-pol trs-rate
# Rp: p-pol ref-rate
# Ts: s-pol trs-rate
# Rs: s-pol ref-rate
Rp = rp**2
Tp = tp**2 * (n2 * np.cos(t2)) / (n1 * c1)
Rs = rs**2
Ts = ts**2 * (n2 * np.cos(t2)) / (n1 * c1)

pol_plot(t1Deg, tp, rp, ts, rs)
pol_plot(t1Deg, Rp, Tp, Rs, Ts)
plt.show()
