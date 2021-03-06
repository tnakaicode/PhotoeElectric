import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import pi,sin,cos,tan,arcsin,linspace,arange,sqrt
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid,axis,tight_layout

n1 = 1.5            # 媒質１の屈折率
n2 = 1.0            # 媒質２の屈折率
ep1 = n1**2         # 媒質１の誘電率
ep2 = n2**2         # 媒質２の誘電率

t1Deg = linspace(0, 90, 90)  # 入射角 t1 の配列の生成
t1 = t1Deg /180*pi   # 入射角をラジアンに直す
s1 = sin(t1)         # sin(t1) 
c1 = cos(t1)         # cos(t1) 
s2 = n1/n2*s1        # sin(t1) 
c2 = sqrt(1-s2**2)   # cos(t2) i
n1z = n1*c1         # n1z=k1z/k0  
n2z = n2*c2         # n2z=k1z/k0  

rs = (n1z-n2z)/(n1z+n2z)                  # s-偏光反射係数
rp = (ep2*n1z-ep1*n2z)/(ep2*n1z+ep1*n2z)  # p-偏光反射係数

RsAbs = abs(rs)**2  # s-偏光反射率
RpAbs = abs(rp)**2  # p-偏光反射率

plt.figure(figsize=(8,6))
plot(t1Deg,RpAbs, label=r"$R_{12}^{\rm{p}}$",linewidth = 3.0, color='black') # p-偏光反射率のプロット
plot(t1Deg,RsAbs, label=r"$R_{12}^{\rm{s}}$",linewidth = 3.0, color='gray') #  s-偏光反射率のプロット
xlabel(r"入射角 (deg.)",fontsize=20)   # x 軸のラベル
ylabel(r"反射・透過係数",fontsize=20)              # y 軸のラベル
title("反射率",fontsize=20)          # グラフタイトル
grid(True)                                 # グリッドを表示
axis([0.0,90,0,1.1])                       # プロット範囲
legend(fontsize=20,loc='lower right')     # 凡例の表示とフォントサイズ
plt.tick_params(labelsize=20)     # 軸の目盛表示とフォントサイズの指定
tight_layout()  # 枠に収まるようなグラフにするコマンド
show()                                     # グラフを表示
