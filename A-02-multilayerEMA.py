import numpy as np
import matplotlib.pyplot as plt


def func_nAg(WLs):
    ep = 3.691 - 9.1522**2 / ((1240 / WLs)**2 + 1j * 0.021 * (1240 / WLs))
    index = np.sqrt(ep)
    return index


def func_nTiO2(WLs):
    ep = 5.193 + 0.244 / ((WLs / 1000)**2 - 0.0803)
    index = np.sqrt(ep)
    return index


def mMATs(n1z, n2z):
    # Mij np.matrix for s-polarization
    return (1 / (2 * n1z)) * np.matrix([[n1z + n2z, n1z - n2z], [n1z - n2z, n1z + n2z]])


def mMATp(n1z, n2z, n1, n2):
    return (1 / (2 * n1 * n2 * n1z)) * np.matrix([[n1**2 * n2z + n2**2 * n1z, n1**2 * n2z - n2**2 * n1z], [n1**2 * n2z - n2**2 * n1z, n1**2 * n2z + n2**2 * n1z]])
    # Mij np.matrix for p-polarization


def matFAI1(n1z, d1, k0):
    # fai np.matrix for s-polarization
    return np.matrix([[np.exp(1j * n1z * k0 * d1), 0], [0, np.exp(-1j * n1z * k0 * d1)]])


def matPI1(n1, n1z):
    return np.matrix([[n1z / n1, -n1z / n1, 0, 0], [n1, n1, 0, 0], [0, 0, 1, 1], [0, 0, n1z, -n1z]])


def matPI2(n2o, n2oz, n2ez, c2dash):
    return np.matrix([[c2dash, -c2dash, 0, 0], [n2o**2 / n2ez * c2dash, n2o**2 / n2ez * c2dash, 0, 0], [0, 0, 1, 1], [0, 0, n2oz, -n2oz]])


def matPI2inv(n2o, n2oz, n2ez, c2dash):
    return 0.5 * np.matrix([[1 / c2dash, n2ez / (n2o**2 * c2dash), 0, 0], [-1 / c2dash, n2ez / (n2o**2 * c2dash), 0, 0], [0, 0, 1, 1 / n2oz], [0, 0, 1, -1 / n2oz]])


def matPI3(n3, n3z):
    return np.matrix([[n3z / n3, -n3z / n3, 0, 0], [n3, n3, 0, 0], [0, 0, 1, 1], [0, 0, n3z, -n3z]])


def matPI3inv(n3, n3z):
    return 0.5 * np.matrix([[n3 / n3z, 1 / n3, 0, 0], [-n3 / n3z, 1 / n3, 0, 0], [0, 0, 1, 1 / n3z], [0, 0, 1, -1 / n3z]])


def matFAI2(k0, n2ez, n2oz, d2):
    return np.matrix([[np.exp(1j * n2ez * k0 * d2), 0, 0, 0], [0, np.exp(-1j * n2ez * k0 * d2), 0, 0], [0, 0, np.exp(1j * n2oz * k0 * d2), 0], [0, 0, 0, np.exp(-1j * n2oz * k0 * d2)]])

############# 初期化 ##############


WLmin = 300    # 始めの波長
WLmax = 1000   # 終わりの波長
WLperiod = 1   # 波長間隔
WLx = np.arange(WLmin, WLmax + 1, WLperiod)   # 波長の配列
NumWLx = int((WLmax - WLmin) / WLperiod) + 1   # 波長の配列数
k0 = 2 * np.pi / WLx   # 波数の配列数

t1Deg = 45  # 入射角
t1 = t1Deg / 180 * np.pi   # 入射角をラジアンに直す

############# 多層膜A モデルの計算 ##############

n1 = 1
nA = 1
n2 = np.zeros(NumWLx, dtype=complex)
n3 = np.zeros(NumWLx, dtype=complex)
n4 = np.zeros(NumWLx, dtype=complex)
n5 = np.zeros(NumWLx, dtype=complex)
n6 = np.zeros(NumWLx, dtype=complex)
n7 = np.zeros(NumWLx, dtype=complex)
n8 = np.zeros(NumWLx, dtype=complex)
n9 = np.zeros(NumWLx, dtype=complex)
d2 = d4 = d6 = d8 = 10
d3 = d5 = d7 = d9 = 10

for i in range(NumWLx):
    n2[i] = n4[i] = n6[i] = n8[i] = func_nAg(WLx[i])
    n3[i] = n5[i] = n7[i] = n9[i] = func_nTiO2(WLx[i])

s1 = np.sin(t1)
c1 = np.cos(t1)
s2, s3, s4, s5, s6, s7, s8, s9, sA = n1 / n2 * s1, n1 / n3 * s1, n1 / n4 * s1, n1 / \
    n5 * s1, n1 / n6 * s1, n1 / n7 * s1, n1 / n8 * s1, n1 / n9 * s1, n1 / nA * s1
c2, c3, c4, c5, c6, c7, c8, c9, cA = np.sqrt(1 - s2**2), np.sqrt(1 - s3**2), np.sqrt(1 - s4**2), np.sqrt(1 - s5**2), np.sqrt(1 - s6**2), np.sqrt(1 - s7**2), \
    np.sqrt(1 - s8**2), np.sqrt(1 - s9**2), np.sqrt(1 - sA**2),
n1z, n2z, n3z, n4z, n5z, n6z, n7z, n8z, n9z, nAz = n1 * c1, n2 * c2, n3 * \
    c3, n4 * c4, n5 * c5, n6 * c6, n7 * c7, n8 * c8, n9 * c9, nA * cA

matTs = np.zeros((NumWLx, 2, 2), dtype=complex)          # s-偏光伝搬行列の初期化
matTp = np.zeros((NumWLx, 2, 2), dtype=complex)          # p-偏光伝搬行列の初期化
rsML1 = np.zeros((NumWLx), dtype=complex)              # s-偏光反射係数の初期化
tsML1 = np.zeros((NumWLx), dtype=complex)              # s-偏光透過係数の初期化
rpML1 = np.zeros((NumWLx), dtype=complex)              # p-偏光反射係数の初期化
tpML1 = np.zeros((NumWLx), dtype=complex)              # p-偏光透過係数の初期化

for i in range(NumWLx):

    matTs[i] = mMATs(nAz, n9z[i])@matFAI1(n9z[i], d9, k0[i])@mMATs(n9z[i], n8z[i])@matFAI1(n8z[i], d8, k0[i])@mMATs(n8z[i], n7z[i]) \
        @matFAI1(n7z[i], d7, k0[i])@mMATs(n7z[i], n6z[i])@matFAI1(n6z[i], d6, k0[i])@mMATs(n6z[i], n5z[i])@matFAI1(n5z[i], d5, k0[i]) \
        @mMATs(n5z[i], n4z[i])@matFAI1(n4z[i], d4, k0[i])@mMATs(n4z[i], n3z[i])@matFAI1(n3z[i], d3, k0[i])@mMATs(n3z[i], n2z[i]) \
        @matFAI1(n2z[i], d2, k0[i])@mMATs(n2z[i], n1z)
    # s-偏光伝搬行列計算
    matTp[i] = mMATp(nAz, n9z[i], nA, n9[i])@matFAI1(n9z[i], d9, k0[i])@mMATp(n9z[i], n8z[i], n9[i], n8[i])@matFAI1(n8z[i], d8, k0[i]) \
        @mMATp(n8z[i], n7z[i], n8[i], n7[i])@matFAI1(n7z[i], d7, k0[i])@mMATp(n7z[i], n6z[i], n7[i], n6[i])@matFAI1(n6z[i], d6, k0[i]) \
        @mMATp(n6z[i], n5z[i], n6[i], n5[i])@matFAI1(n5z[i], d5, k0[i])@mMATp(n5z[i], n4z[i], n5[i], n4[i])@matFAI1(n4z[i], d4, k0[i]) \
        @mMATp(n4z[i], n3z[i], n4[i], n3[i])@matFAI1(n3z[i], d3, k0[i])@mMATp(n3z[i], n2z[i], n3[i], n2[i])@matFAI1(n2z[i], d2, k0[i]) \
        @mMATp(n2z[i], n1z, n2[i], n1)
    # p-偏光伝搬行列計算

    # s-偏光反射係数(多層膜 A)
    rsML1[i] = -matTs[i, 1, 0] / matTs[i, 1, 1]
    tsML1[i] = matTs[i, 0, 0] - matTs[i, 0, 1] * \
        matTs[i, 1, 0] / matTs[i, 1, 1]     # s-偏光透過係数(多層膜 A)
    # p-偏光反射係数(多層膜 A)
    rpML1[i] = -matTp[i, 1, 0] / matTp[i, 1, 1]
    tpML1[i] = matTp[i, 0, 0] - matTp[i, 0, 1] * \
        matTp[i, 1, 0] / matTp[i, 1, 1]     # p-偏光透過係数(多層膜 A)

RsML1 = abs(rsML1)**2        # s-偏光反射率(多層膜 A)
RpML1 = abs(rpML1)**2        # p-偏光反射率(多層膜 A)
TsML1 = abs(tsML1)**2        # s-偏光透過率(多層膜 A)
TpML1 = abs(tpML1)**2        # p-偏光透過率(多層膜 A)

########## 多層膜B モデルの計算 ############

n1 = 1
nA = 1
n2 = np.zeros(NumWLx, dtype=complex)
n3 = np.zeros(NumWLx, dtype=complex)
n4 = np.zeros(NumWLx, dtype=complex)
n5 = np.zeros(NumWLx, dtype=complex)
n6 = np.zeros(NumWLx, dtype=complex)
n7 = np.zeros(NumWLx, dtype=complex)
n8 = np.zeros(NumWLx, dtype=complex)
n9 = np.zeros(NumWLx, dtype=complex)
d2 = d4 = d6 = d8 = 10
d3 = d5 = d7 = d9 = 10

for i in range(NumWLx):
    n2[i] = n4[i] = n6[i] = n8[i] = func_nTiO2(WLx[i])
    n3[i] = n5[i] = n7[i] = n9[i] = func_nAg(WLx[i])

s1 = np.sin(t1)
c1 = np.cos(t1)
s2, s3, s4, s5, s6, s7, s8, s9, sA = n1 / n2 * s1, n1 / n3 * s1, n1 / n4 * s1, n1 / \
    n5 * s1, n1 / n6 * s1, n1 / n7 * s1, n1 / n8 * s1, n1 / n9 * s1, n1 / nA * s1
c2, c3, c4, c5, c6, c7, c8, c9, cA = np.sqrt(1 - s2**2), np.sqrt(1 - s3**2), np.sqrt(1 - s4**2), np.sqrt(1 - s5**2), np.sqrt(1 - s6**2), np.sqrt(1 - s7**2), \
    np.sqrt(1 - s8**2), np.sqrt(1 - s9**2), np.sqrt(1 - sA**2),
n1z, n2z, n3z, n4z, n5z, n6z, n7z, n8z, n9z, nAz = n1 * c1, n2 * c2, n3 * \
    c3, n4 * c4, n5 * c5, n6 * c6, n7 * c7, n8 * c8, n9 * c9, nA * cA

matTs = np.zeros((NumWLx, 2, 2), dtype=complex)          # s-偏光伝搬行列初期化
matTp = np.zeros((NumWLx, 2, 2), dtype=complex)          # p-偏光伝搬行列初期化
rsML2 = np.zeros((NumWLx), dtype=complex)              # rs 初期化
tsML2 = np.zeros((NumWLx), dtype=complex)              # ts 初期化
rpML2 = np.zeros((NumWLx), dtype=complex)              # rp 初期化
tpML2 = np.zeros((NumWLx), dtype=complex)              # tp 初期化

for i in range(NumWLx):

    matTs[i] = mMATs(nAz, n9z[i])@matFAI1(n9z[i], d9, k0[i])@mMATs(n9z[i], n8z[i])@matFAI1(n8z[i], d8, k0[i])@mMATs(n8z[i], n7z[i]) \
        @matFAI1(n7z[i], d7, k0[i])@mMATs(n7z[i], n6z[i])@matFAI1(n6z[i], d6, k0[i])@mMATs(n6z[i], n5z[i])@matFAI1(n5z[i], d5, k0[i]) \
        @mMATs(n5z[i], n4z[i])@matFAI1(n4z[i], d4, k0[i])@mMATs(n4z[i], n3z[i])@matFAI1(n3z[i], d3, k0[i])@mMATs(n3z[i], n2z[i]) \
        @matFAI1(n2z[i], d2, k0[i])@mMATs(n2z[i], n1z)
    # s-偏光伝搬行列の計算
    matTp[i] = mMATp(nAz, n9z[i], nA, n9[i])@matFAI1(n9z[i], d9, k0[i])@mMATp(n9z[i], n8z[i], n9[i], n8[i])@matFAI1(n8z[i], d8, k0[i]) \
        @mMATp(n8z[i], n7z[i], n8[i], n7[i])@matFAI1(n7z[i], d7, k0[i])@mMATp(n7z[i], n6z[i], n7[i], n6[i])@matFAI1(n6z[i], d6, k0[i]) \
        @mMATp(n6z[i], n5z[i], n6[i], n5[i])@matFAI1(n5z[i], d5, k0[i])@mMATp(n5z[i], n4z[i], n5[i], n4[i])@matFAI1(n4z[i], d4, k0[i]) \
        @mMATp(n4z[i], n3z[i], n4[i], n3[i])@matFAI1(n3z[i], d3, k0[i])@mMATp(n3z[i], n2z[i], n3[i], n2[i])@matFAI1(n2z[i], d2, k0[i]) \
        @mMATp(n2z[i], n1z, n2[i], n1)
    # p-偏光伝搬行列の計算

    rsML2[i] = -matTs[i, 1, 0] / matTs[i, 1, 1]    # s-偏光反射係数(多層膜モデル)
    tsML2[i] = matTs[i, 0, 0] - matTs[i, 0, 1] * \
        matTs[i, 1, 0] / matTs[i, 1, 1]     # s-偏光透過係数(多層膜モデル)
    rpML2[i] = -matTp[i, 1, 0] / matTp[i, 1, 1]    # p-偏光反射係数(多層膜モデル)
    tpML2[i] = matTp[i, 0, 0] - matTp[i, 0, 1] * \
        matTp[i, 1, 0] / matTp[i, 1, 1]     # s-偏光透過係数(多層膜モデル)

RsML2 = abs(rsML2)**2        # s-偏光反射率(多層膜 B)
RpML2 = abs(rpML2)**2        # p-偏光反射率(多層膜 B)
TsML2 = abs(tsML2)**2        # s-偏光透過率(多層膜 B)
TpML2 = abs(tpML2)**2        # p-偏光透過率(多層膜 B)

##### EMA モデルの計算(異方性をもつ薄膜の 3 層問題)######

n1 = 1
n3 = 1
d2 = 80

nTiO2 = np.zeros((NumWLx), dtype=complex)
nAg = np.zeros((NumWLx), dtype=complex)

for i in range(NumWLx):
    nTiO2[i] = func_nTiO2(WLx[i])
    nAg[i] = func_nAg(WLx[i])

epx = 0.5 * (nTiO2**2 + nAg**2)
epz = 2 * (nTiO2**2) * (nAg**2) / ((nTiO2**2) + (nAg**2))

no = np.sqrt(epx)
ne = np.sqrt(epz)

s1 = np.sin(t1)
c1 = np.cos(t1)
kappa = n1 * s1
s3 = n1 / n3 * s1
c3 = np.sqrt(1 - s3**2)

n1z = n1 * c1
n3z = n3 * c3

n2oz = np.sqrt(no**2 - kappa**2)
n2ez = (no / ne) * np.sqrt(ne**2 - kappa**2)
c2dash = (ne * np.sqrt(ne**2 - kappa**2)) / \
    np.sqrt(ne**4 + kappa**2 * (no**2 - ne**2))

matT = np.zeros((NumWLx, 4, 4), dtype=complex)             # s-偏光伝搬行列初期化
rsEMA = np.zeros((NumWLx), dtype=complex)                 # rs 初期化
tsEMA = np.zeros((NumWLx), dtype=complex)                 # ts 初期化
rpEMA = np.zeros((NumWLx), dtype=complex)                 # rp 初期化
tpEMA = np.zeros((NumWLx), dtype=complex)                 # tp 初期化

for i in range(NumWLx):
    matT[i] = matPI3inv(n3, n3z)@matPI2(no[i], n2oz[i], n2ez[i], c2dash[i])@matFAI2(
        k0[i], n2ez[i], n2oz[i], d2)@matPI2inv(no[i], n2oz[i], n2ez[i], c2dash[i])@matPI1(n1, n1z)
    rsEMA[i] = -matT[i, 3, 2] / matT[i, 3, 3]
    tsEMA[i] = matT[i, 2, 2] - matT[i, 2, 3] * matT[i, 3, 2] / matT[i, 3, 3]
    rpEMA[i] = -matT[i, 1, 0] / matT[i, 1, 1]
    tpEMA[i] = matT[i, 0, 0] - matT[i, 0, 1] * matT[i, 1, 0] / matT[i, 1, 1]

RsEMA = abs(rsEMA)**2        # s-偏光反射率(EMA モデル
RpEMA = abs(rpEMA)**2        # p-偏光反射率(EMA モデル)
TsEMA = abs(tsEMA)**2        # s-偏光透過率(EMA モデル)
TpEMA = abs(tpEMA)**2        # p-偏光透過率(EMA モデル)


############# PLOT ##############

plt.figure(figsize=(8, 6))
plt.plot(WLx, RsML1, label="RsML1", linewidth=3.0, color='black')
plt.plot(WLx, RsML2, label="RsML2", linewidth=3.0, color='gray')
plt.plot(WLx, RsEMA, label="RsEMA", linewidth=3.0,
         color='black', linestyle='dashed')
plt.xlabel(r"Wavelength(nm)", fontsize=22)
plt.ylabel(r"reflectivity", fontsize=22)
plt.title("", fontsize=22)
plt.grid(True)
plt.axis([300, 1000, 0, 1.1])
plt.legend(fontsize=20, loc='lower right')
plt.tick_params(labelsize=20)
plt.tight_layout()

plt.figure(figsize=(8, 6))
plt.plot(WLx, RpML1, label="RpML1", linewidth=3.0, color='black')
plt.plot(WLx, RpML2, label="RpML2", linewidth=3.0, color='gray')
plt.plot(WLx, RpEMA, label="RpEMA", linewidth=3.0,
         color='black', linestyle='dashed')
plt.xlabel(r"Wavelength(nm)", fontsize=22)
plt.ylabel(r"reflectivity", fontsize=22)
plt.title("", fontsize=22)
plt.grid(True)
plt.axis([300, 1000, 0, 1.1])
plt.legend(fontsize=20, loc='lower right')
plt.tick_params(labelsize=20)
plt.tight_layout()

plt.figure(figsize=(8, 6))
plt.plot(WLx, TsML1, label="TsML1", linewidth=3.0, color='black')
plt.plot(WLx, TsML2, label="TsML2", linewidth=3.0, color='gray')
plt.plot(WLx, TsEMA, label="TsEMA", linewidth=3.0,
         color='black', linestyle='dashed')
plt.xlabel(r"Wavelength(nm)", fontsize=22)
plt.ylabel(r"transmittance", fontsize=22)
plt.title("", fontsize=22)
plt.grid(True)
plt.axis([300, 1000, 0, 1.1])
plt.legend(fontsize=20, loc='lower right')
plt.tick_params(labelsize=20)
plt.tight_layout()

plt.figure(figsize=(8, 6))
plt.plot(WLx, TpML1, label="TpML1", linewidth=3.0, color='black')
plt.plot(WLx, TpML2, label="TpML2", linewidth=3.0, color='gray')
plt.plot(WLx, TpEMA, label="TpEMA", linewidth=3.0,
         color='black', linestyle='dashed')
plt.xlabel(r"Wavelength(nm)", fontsize=22)
plt.ylabel(r"transmittance", fontsize=22)
plt.title("", fontsize=22)
plt.grid(True)
plt.axis([300, 1000, 0, 1.1])
plt.legend(fontsize=20, loc='lower right')
plt.tick_params(labelsize=20)
plt.tight_layout()
plt.show()
