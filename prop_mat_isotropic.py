import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
import glob
import imageio
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base import plot2d
from src.index import ref_index, ref_ratio, ref_ratio_complex


def mMATs(n1z, n2z):
    """
    S-Pol Mij Boundary Condition
    """
    mat = np.empty((2, 2), dtype=complex)
    mat[0, 0] = n1z + n2z
    mat[0, 1] = n1z - n2z
    mat[1, 0] = n1z - n2z
    mat[1, 1] = n1z + n2z
    return (1 / (2 * n1z)) * mat


def mMATp(n1z, n2z, n1, n2):
    """
    P-Pol Mij Boundary Condition
    """
    mat = np.empty((2, 2), dtype=complex)
    mat[0, 0] = n1**2 * n2z + n2**2 * n1z
    mat[0, 1] = n1**2 * n2z - n2**2 * n1z
    mat[1, 0] = n1**2 * n2z - n2**2 * n1z
    mat[1, 1] = n1**2 * n2z + n2**2 * n1z
    return (1 / (2 * n1 * n2 * n1z)) * mat


def matFAI(n1z, d1, k0):
    """
    Phi-Matrix i-index Phase Trasition
    """
    mat = np.empty((2, 2), dtype=complex)
    mat[0, 0] = np.exp(1j * n1z * k0 * d1)
    mat[0, 1] = 0
    mat[1, 0] = 0
    mat[1, 1] = np.exp(-1j * n1z * k0 * d1)
    return mat


n1 = 1.86           # 媒質 1(プリズム)の屈折率
n2 = np.sqrt(-10.8 + 1j * 1.47, dtype=complex)    # 媒質 2(金)の屈折率
n3 = 1.5            # 媒質 3(誘電体薄膜)の屈折率
n4 = 1.33           # 媒質 4(水)の屈折率
ep1 = n1**2         # 媒質 1 の誘電率
ep2 = n2**2         # 媒質 2 の誘電率
ep3 = n3**2         # 媒質 3 の誘電率
ep4 = n4**2         # 媒質 4 の誘電率
d2 = 47             # 媒質 2(金)の厚さ d2〔nm〕
d3 = 10             # 媒質 2(誘電体)の厚さ d2〔nm〕
WL = 633            # 真空中の波長 WL〔nm〕
k0 = 2 * np.pi / WL        # 真空中の波数

t1start = 40        # 始めの角度
t1end = 70          # 終わりの角度
t1points = 300      # プロット数

t1DegOut = np.linspace(t1start, t1end, t1points)
t1 = 0.25 * np.pi + (1 / n1) * np.arcsin((t1DegOut - 45) / 180 * np.pi)
s1 = np.sin(t1)
c1 = np.cos(t1)
s2 = n1 / n2 * s1
c2 = np.sqrt(1 - s2**2, dtype=complex)
s3 = n1 / n3 * s1
c3 = np.sqrt(1 - s3**2, dtype=complex)
s4 = n1 / n4 * s1
c4 = np.sqrt(1 - s4**2, dtype=complex)

n1z = n1 * c1
n2z = n2 * c2
n3z = n3 * c3
n4z = n4 * c4

matT0 = np.zeros((t1points, 2, 2), dtype=complex)
matT1 = np.zeros((t1points, 2, 2), dtype=complex)
r0 = np.zeros((t1points), dtype=complex)
r1 = np.zeros((t1points), dtype=complex)

for i in range(t1points):
    matT0[i] = mMATp(n4z[i], n2z[i], n4, n2) @ matFAI(n2z[i],
                                                      d2, k0) @ mMATp(n2z[i], n1z[i], n2, n1)
    matT1[i] = mMATp(n4z[i], n3z[i], n4, n3) @ matFAI(n3z[i], d3, k0) @ mMATp(
        n3z[i], n2z[i], n3, n2) @ matFAI(n2z[i], d2, k0) @ mMATp(n2z[i], n1z[i], n2, n1)
    r0[i] = -matT0[i, 1, 0] / matT0[i, 1, 1]
    r1[i] = -matT1[i, 1, 0] / matT1[i, 1, 1]

R0Abs = abs(r0)**2
R1Abs = abs(r1)**2

obj = plot2d(aspect="auto")
obj.axs.plot(t1DegOut, R1Abs, label="R1", linewidth=3.0, color='gray')
obj.axs.plot(t1DegOut, R0Abs, label="R0", linewidth=3.0, color='black')
obj.axs.set_xlabel(r"$\theta_1$ (deg.)", fontsize=20)
obj.axs.set_ylabel(r"Reflectivity", fontsize=20)
obj.axs.set_title("Surface Plasmon Resonance", fontsize=20)
obj.axs.legend(fontsize=20, loc='lower right')
obj.SavePng()
