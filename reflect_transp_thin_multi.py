import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

n1 = 1.0            # 媒質１の屈折率
n2 = 1.5            # 媒質２の屈折率
n3 = 1.0            # 媒質３の屈折率
ep1 = n1**2         # 媒質１の誘電率
ep2 = n2**2         # 媒質２の誘電率
ep3 = n3**2         # 媒質３の誘電率
WL = 500            # 真空中の波長 WL〔nm〕
k0 = 2 * sp.pi / WL        # 真空中の波数

d2 = sp.linspace(0, 500, 500)      # 媒質 2 の厚さの配列  0 nm から 500 nm の 501 個

t1Deg = 0  # 入射角 t1
t1 = t1Deg / 180 * sp.pi  # 入射角をラジアンに直す
s1 = sp.sin(t1)          # sin(t1)
c1 = sp.cos(t1)          # cos(t1)
s2 = n1 / n2 * s1        # sin(t1)
c2 = sp.sqrt(1 - s2**2)  # cos(t2)
s3 = n1 / n3 * s1        # sin(t1)
c3 = sp.sqrt(1 - s3**2)  # cos(t3)

n1z = n1 * c1         # n1z=k1z/k0
n2z = n2 * c2         # n2z=k1z/k0
n3z = n3 * c3         # n2z=k1z/k0

rs12 = (n1z - n2z) / (n1z + n2z)         # s-偏光反射係数 rs12
rp12 = (ep2 * n1z - ep1 * n2z) / (ep2 * n1z + ep1 * n2z)  # p-偏光反射係数 rp12
rs23 = (n2z - n3z) / (n2z + n3z)         # rs-偏光反射係数 rs23
rp23 = (ep3 * n2z - ep2 * n3z) / (ep3 * n2z + ep2 * n3z)  # p-偏光反射係数 rp23

rs = (rs12 + rs23 * sp.exp(2 * 1j * n2z * k0 * d2)) / \
    (1 + rs23 * rs12 * sp.exp(2 * 1j * n2z * k0 * d2))
rp = (rp12 + rp23 * sp.exp(2 * 1j * n2z * k0 * d2)) / \
    (1 + rp23 * rp12 * sp.exp(2 * 1j * n2z * k0 * d2))

RsAbs = abs(rs)**2        # s-偏光反射率
RpAbs = abs(rp)**2        # p-偏光反射率

# p-偏光反射率のプロット
plt.plot(d2, RpAbs, label="$R_p$", linewidth=3.0, color='black')
plt.xlabel(r"tin (nm)", fontsize=20)
plt.ylabel("injection (deg.)", fontsize=20)
plt.title("reflect rate (depend on depth)", fontsize=20)
plt.grid(True)
plt.axis([0.0, 500, 0, 0.2])
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
