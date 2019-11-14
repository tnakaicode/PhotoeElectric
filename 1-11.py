import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import pi, arange, sqrt, zeros, array


def func_nAg(WLs):
    ep = 3.691 - 9.1522**2 / ((1240 / WLs)**2 + 1j * 0.021 * (1240 / WLs))
    index = sqrt(ep)
    return index        # 銀の誘電関数


def func_nTiO2(WLs):
    ep = 5.193 + 0.244 / ((WLs / 1000)**2 - 0.0803)
    index = sqrt(ep)
    return index       # TiO2の誘電関数


WLmin = 300    # 波長(短波長側) 〔nm〕
WLmax = 1000  # 波長(長波長側) 〔nm〕
WLperiod = 1  # 波長間隔  〔nm〕
WLx = arange(WLmin, WLmax + 1, WLperiod)  # 波長の配列
NumWLx = int((WLmax - WLmin) / WLperiod) + 1  # 波長の数
k0 = 2 * pi / WLx  # 各波長の波数

nTiO2 = zeros((NumWLx), dtype=complex)  # Ti02 屈折率の配列の初期化
nAg = zeros((NumWLx), dtype=complex)  # Ag 屈折率の配列の初期化

for i in range(NumWLx):
    nTiO2[i] = func_nTiO2(WLx[i])  # Ti02 屈折率の生成
    nAg[i] = func_nAg(WLx[i])      # Ag 屈折率の生成

epx = 0.5 * (nTiO2**2 + nAg**2)       # EMA による誘電率の計算 x 方向
epz = 2 * (nTiO2**2) * (nAg**2) / ((nTiO2**2) + (nAg**2))  # EMA による誘電率の計算 z 方向

plt.figure(figsize=(8, 6))
# x 方向有効媒質の誘電率(実部)のプロット
plt.plot(WLx, epx.real,
         label=r"Re$(\epsilon_{\rm \parallel})$", linewidth=3.0, color='black')
# x 方向有効媒質の誘電率(虚部)のプロット
plt.plot(WLx, epx.imag,
         label=r"Im$(\epsilon_{\rm \parallel})$", linewidth=3.0, color='gray')
plt.xlabel(r"Wavelength(nm)", fontsize=20)                 # x 軸のラベル
plt.ylabel(r"$ \epsilon_{\rm \parallel}$", fontsize=20)    # y 軸のラベル
plt.title("", fontsize=20)           # グラフのタイトル
plt.grid(True)                      # グリッドを表示
plt.axis([300, 1000, -30, 30])         # プロット範囲
plt.legend(fontsize=20, loc='upper right')   # 凡例の表示
plt.tick_params(labelsize=20)           # 軸の目盛の表示とフォントサイズの指定

plt.figure(figsize=(8, 6))
# z 方向有効媒質の誘電率(実部)のプロット
plt.plot(WLx, epz.real,
         label=r"Re$(\epsilon_{\rm z})$", linewidth=3.0, color='black')
# z 方向有効媒質の誘電率(虚部)のプロット
plt.plot(WLx, epz.imag,
         label=r"Im$(\epsilon_{\rm z})$", linewidth=3.0, color='gray')
plt.xlabel(r"Wavelength(nm)", fontsize=20)           # x 軸のラベル
plt.ylabel(r"$ \epsilon_{\rm z}$", fontsize=20)      # y 軸のラベル
plt.title("", fontsize=20)                           # グラフのタイトル
plt.grid(True)                                      # グリッドを表示
plt.axis([300, 1000, -100, 100])                       # プロット範囲
plt.legend(fontsize=20, loc='lower right')           # 凡例の表示
plt.tick_params(labelsize=20)                   # 軸の目盛表示とフォントサイズの指定
plt.show()                  # グラフを表示
