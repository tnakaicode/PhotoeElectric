import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import pi, sin, cos, tan, arcsin, linspace, arange, sqrt

from base import pol_plot, pol_plot_ref
from material import mat_refract

n1 = mat_refract["Diamond"]
n2 = 1.0
ep1 = n1**2         # 媒質１の誘電率
ep2 = n2**2         # 媒質２の誘電率

t1Deg = linspace(0, 90, 90)
t1 = t1Deg / 180 * pi
s1 = sin(t1)
c1 = cos(t1)
s2 = n1 / n2 * s1
c2 = sqrt(1 - s2**2)
n1z = n1 * c1
n2z = n2 * c2

rs = (n1z - n2z) / (n1z + n2z)                  # s-偏光反射係数
rp = (ep2 * n1z - ep1 * n2z) / (ep2 * n1z + ep1 * n2z)  # p-偏光反射係数

RsAbs = abs(rs)**2  # s-偏光反射率
RpAbs = abs(rp)**2  # p-偏光反射率

pol_plot_ref(t1Deg, RpAbs, RsAbs)
plt.show()
