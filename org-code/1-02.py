import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import pi,sin,cos,tan,arcsin,linspace,arange
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid,axis,tight_layout

n1 = 1            # 媒質１の屈折率
n2 = 1.5          # 媒質２の屈折率

t1Deg = linspace(0, 90, 90)  # 入射角 t1 の配列の生成
t1  = t1Deg /180*pi          # 入射角をラジアンに直す
t2 = arcsin((n1/n2)*sin(t1))   # 屈折角 t2 を求める

tp = 2*n1*cos(t1)/(n2*cos(t1)+n1*cos(t2))               # tp: p-偏光透過係数
rp = (n2*cos(t1)-n1*cos(t2))/(n2*cos(t1)+n1*cos(t2))    # rp: p-偏光反射係数
ts = 2*n1*cos(t1)/(n1*cos(t1)+n2*cos(t2))              # ts: s-偏光透過係数
rs = (n1*cos(t1)-n2*cos(t2))/(n1*cos(t1)+n2*cos(t2))   # rs: s-偏光反射係数

Rp = rp**2                                # p-偏光反射率
Tp = tp**2*(n2*cos(t2))/(n1*cos(t1))      # p-偏光透過率
Rs = rs**2                                # s-偏光反射率
Ts = ts**2*(n2*cos(t2))/(n1*cos(t1))      # s-偏光透過率

plt.figure(figsize=(8,6))    # figure size    
plot(t1Deg,Rp, label=r"$R_{12}^{\rm{p}}$",linewidth = 3.0, color='black', linestyle='dashed')       #  Rp をプロット
plot(t1Deg,Tp, label=r"$T_{12}^{\rm{p}}$",linewidth = 3.0, color='black')  # Tp をプロット
plot(t1Deg,Rs, label=r"$R_{12}^{\rm{s}}$",linewidth = 3.0, color='gray', linestyle='dashed')       # Rs をプロット
plot(t1Deg,Ts, label=r"$T_{12}^{\rm{s}}$",linewidth = 3.0, color='gray')   # Ts をプロット

xlabel(r"入射角 (deg.)",fontsize=20)   # x 軸のラベル
ylabel(r"反射率・透過率",fontsize=20)              # y 軸のラベル
title("反射率・透過率",fontsize=18)         # グラフタイトル
grid(True)                                # グリッドを表示
axis([0.0,90,0,1.1])                      # プロット範囲 
legend(fontsize=20,loc='lower left')      # 凡例を表示
plt.tick_params(labelsize=20) 
tight_layout()
show()                                    # グラフを表示

