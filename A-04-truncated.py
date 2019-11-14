import RI
import numpy as np
import scipy as sp
import scipy.special
import math
import cmath
import matplotlib.pyplot as plt

from scipy import pi, sin, cos, tan, arcsin, exp, linspace, arange, sqrt, zeros, array, matrix, asmatrix, real, imag, interpolate, integrate
from scipy.special import spherical_jn, spherical_yn, factorial, lpmv, eval_legendre


def kjo(k):
    return math.factorial(k)


def perpen(k, j, r0, ep1, ep2, ep3):
    return ((ep2 - ep1) * (ep1 - ep3) * k * kjo(k + j)) / ((ep2 + ep1) * ((k + 1) * ep1 + k * ep3) * kjo(k) * kjo(j) * (2 * r0)**(k + j + 1))


def parallel(k, j, r0, ep1, ep2, ep3):
    return ((ep2 - ep1) * (ep1 - ep3) * k * kjo(k + j)) / ((ep2 + ep1) * ((k + 1) * ep1 + k * ep3) * kjo(k + 1) * kjo(j - 1) * (2 * r0)**(k + j + 1))


def uhen(ep1, ep2, ep3):
    return (ep1 - ep3) / (2 * ep1 + ep3)


def funcIMG(r, r0, t):
    return r * r - 4 * r * r0 * t + 4 * r0 * r0


def funcIMG2(r0, t):
    return 1 + 4 * r0 * r0 - 4 * r0 * t


def funcV(m, j, r, t, r0):
    return funcIMG(r, r0, t)**(-(j + 1) / 2) * lpmv(m, j, (r * t - 2 * r0) * funcIMG(r, r0, t)**(-1 / 2))


def funcV2(m, j, t, r0):
    return funcIMG2(r0, t)**(-0.5 * (3 + j)) * (-(1 + j) * funcIMG2(r0, t) * lpmv(m, j, (t - 2 * r0) / sqrt(funcIMG2(r0, t))) - 2 * (1 + j - m) * r0 * sqrt(funcIMG2(r0, t)) * lpmv(m, j + 1, (t - 2 * r0) / sqrt(funcIMG2(r0, t))))


def funcW(m, j, r, t, r0):
    return funcIMG(r, r0, t)**(j / 2) * lpmv(m, j, (r * t - 2 * r0) * funcIMG(r, r0, t)**(-1 / 2))


def funcW2(m, j, t, r0):
    return funcIMG2(r0, t)**(j / 2 - 1) * ((j - 4 * j * r0 * r0 + 2 * r0 * (t - 2 * r0)) * lpmv(m, j, (t - 2 * r0) / sqrt(funcIMG2(r0, t))) - 2 * (1 + j - m) * r0 * sqrt(funcIMG2(r0, t)) * lpmv(m, j + 1, (t - 2 * r0) / sqrt(funcIMG2(r0, t))))


qq = 11
idosh = 90   # capping angle
r0 = cos(idosh * pi / 180)
rr = 25   # radius

WLmin = RI.WLmin
WLmax = RI.WLmax
WLperiod = RI.WLperiod
WLx = RI.WLx
NumWLx = RI.NumWLx
k0 = 2 * pi / WLx

ep1 = zeros(NumWLx, dtype=complex)
ep2 = zeros(NumWLx, dtype=complex)
ep3 = zeros(NumWLx, dtype=complex)
ep4 = zeros(NumWLx, dtype=complex)

matrixCinteg = zeros([qq, qq], dtype=float)
matrixDinteg = zeros([qq, qq], dtype=float)
matrixEinteg = zeros([qq], dtype=float)
matrixFinteg = zeros([qq, qq], dtype=float)
matrixGinteg = zeros([qq, qq], dtype=float)
matrixJinteg = zeros([qq, qq], dtype=float)
matrixKinteg = zeros([qq, qq], dtype=float)
matrixMinteg = zeros([qq, qq], dtype=float)
matrixNinteg = zeros([qq, qq], dtype=float)
matrixPinteg = zeros([qq], dtype=float)

matrixClist = zeros([NumWLx, qq, qq], dtype=complex)
matrixDlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixElist = zeros([NumWLx, qq], dtype=complex)
matrixFlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixGlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixHlist = zeros([NumWLx, qq], dtype=complex)
matrixJlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixKlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixLlist = zeros([NumWLx, qq], dtype=complex)
matrixMlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixNlist = zeros([NumWLx, qq, qq], dtype=complex)
matrixPlist = zeros([NumWLx, qq], dtype=complex)

matrixRA = zeros([NumWLx, 2 * qq, 2 * qq], dtype=complex)
matrixRB = zeros([NumWLx, 2 * qq, 2 * qq], dtype=complex)
matrixQA = zeros([NumWLx, 2 * qq], dtype=complex)
matrixQB = zeros([NumWLx, 2 * qq], dtype=complex)
vectorA = zeros([NumWLx, 2 * qq], dtype=complex)
vectorB = zeros([NumWLx, 2 * qq], dtype=complex)
alpha_A = zeros([NumWLx], dtype=complex)
alpha_B = zeros([NumWLx], dtype=complex)

ep1 = 1
ep2 = ep4 = 1.5

for i in range(NumWLx):
    ep3[i] = RI.epAu[i]

for k in range(qq):
    matrixEinteg[k], dummy = integrate.quad(
        lambda t: lpmv(0, k + 1, t) * (t - r0), -1, r0)
    matrixPinteg[k], dummy = integrate.quad(
        lambda t: lpmv(1, k + 1, t) * lpmv(1, 1, t), -1, r0)

    for j in range(qq):
        matrixCinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            0, k + 1, t) * (lpmv(0, j + 1, t) - (-1)**(j + 1) * funcV(0, j + 1, 1, t, r0)), -1, r0)
        matrixDinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            0, k + 1, t) * (lpmv(0, j + 1, t) - (-1)**(j + 1) * funcW(0, j + 1, 1, t, r0)), -1, r0)
        matrixFinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            0, k + 1, t) * ((j + 2) * lpmv(0, j + 1, t) - (-1)**(j + 1) * funcV2(0, j + 1, t, r0)), -1, r0)
        matrixGinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            0, k + 1, t) * ((j + 1) * lpmv(0, j + 1, t) + (-1)**(j + 1) * funcW2(0, j + 1, t, r0)), -1, r0)

        matrixJinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            1, k + 1, t) * (lpmv(1, j + 1, t) + (-1)**(j + 1) * funcV(1, j + 1, 1, t, r0)), -1, r0)
        matrixKinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            1, k + 1, t) * (lpmv(1, j + 1, t) + (-1)**(j + 1) * funcW(1, j + 1, 1, t, r0)), -1, r0)
        matrixMinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            1, k + 1, t) * ((j + 2) * lpmv(1, j + 1, t) + (-1)**(j + 1) * funcV2(1, j + 1, t, r0)), -1, r0)
        matrixNinteg[k, j], dummy = integrate.quad(lambda t: lpmv(
            1, k + 1, t) * ((j + 1) * lpmv(1, j + 1, t) - (-1)**(j + 1) * funcW2(1, j + 1, t, r0)), -1, r0)

for i in range(NumWLx):
    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixClist[i, k, j] = 4 * ep1 / \
                    ((ep1 + ep2) * (2 * (k + 1) + 1)) - (ep1 - ep2) / \
                    (ep1 + ep2) * matrixCinteg[k, j]
            else:
                matrixClist[i, k, j] = -(ep1 - ep2) / \
                    (ep1 + ep2) * matrixCinteg[k, j]

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixDlist[i, k, j] = -4 * ep3[i] / ((ep3[i] + ep4) * (2 * (k + 1) + 1)) + (
                    ep3[i] - ep4) / (ep3[i] + ep4) * matrixDinteg[k, j]
            else:
                matrixDlist[i, k, j] = (ep3[i] - ep4) / \
                    (ep3[i] + ep4) * matrixDinteg[k, j]

    for k in range(qq):
        if k == 0:
            matrixElist[i, k] = -2 * ep1 / \
                (3 * ep2) - (1 - ep1 / ep2) * matrixEinteg[k]
        else:
            matrixElist[i, k] = -(1 - ep1 / ep2) * matrixEinteg[k]

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixFlist[i, k, j] = -4 * ep1 * ep2 * (k + 2) / ((ep1 + ep2) * (2 * (k + 1) + 1)) - (
                    ep1 * (ep1 - ep2)) / (ep1 + ep2) * matrixFinteg[k, j]
            else:
                matrixFlist[i, k, j] = - \
                    (ep1 * (ep1 - ep2)) / (ep1 + ep2) * matrixFinteg[k, j]

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixGlist[i, k, j] = -4 * ep3[i] * ep4 * (k + 1) / ((ep3[i] + ep4) * (2 * (k + 1) + 1)) - (
                    ep3[i] * (ep3[i] - ep4)) / (ep3[i] + ep4) * matrixGinteg[k, j]
            else:
                matrixGlist[i, k, j] = - \
                    (ep3[i] * (ep3[i] - ep4)) / \
                    (ep3[i] + ep4) * matrixGinteg[k, j]

    for k in range(qq):
        if k == 0:
            matrixHlist[i, k] = -2 * ep1 / 3
        else:
            matrixHlist[i, k] = 0

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixJlist[i, k, j] = 4 * ep1 * (k + 1) * (k + 2) / ((ep1 + ep2) * (
                    2 * (k + 1) + 1)) - (ep1 - ep2) / (ep1 + ep2) * matrixJinteg[k, j]
            else:
                matrixJlist[i, k, j] = -(ep1 - ep2) / \
                    (ep1 + ep2) * matrixJinteg[k, j]

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixKlist[i, k, j] = -4 * ep3[i] * (k + 1) * (k + 2) / ((ep3[i] + ep4) * (
                    2 * (k + 1) + 1)) + (ep3[i] - ep4) / (ep3[i] + ep4) * matrixKinteg[k, j]
            else:
                matrixKlist[i, k, j] = (ep3[i] - ep4) / \
                    (ep3[i] + ep4) * matrixKinteg[k, j]

    for k in range(qq):
        if k == 0:
            matrixLlist[i, k] = -4 / 3
        else:
            matrixLlist[i, k] = 0

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixMlist[i, k, j] = -4 * ep1 * ep2 * (k + 1) * (k + 2)**2 / (
                    (ep1 + ep2) * (2 * (k + 1) + 1)) - (ep1 * (ep1 - ep2)) / (ep1 + ep2) * matrixMinteg[k, j]
            else:
                matrixMlist[i, k, j] = - \
                    (ep1 * (ep1 - ep2)) / (ep1 + ep2) * matrixMinteg[k, j]

    for k in range(qq):
        for j in range(qq):
            if k == j:
                matrixNlist[i, k, j] = -4 * ep3[i] * ep4 * (k + 1)**2 * (k + 2) / ((ep3[i] + ep4) * (
                    2 * (k + 1) + 1)) - (ep3[i] * (ep3[i] - ep4)) / (ep3[i] + ep4) * matrixNinteg[k, j]
            else:
                matrixNlist[i, k, j] = - \
                    (ep3[i] * (ep3[i] - ep4)) / \
                    (ep3[i] + ep4) * matrixNinteg[k, j]

    for k in range(qq):
        if k == 0:
            matrixPlist[i, k] = -4 * ep2 / 3 - (ep1 - ep2) * matrixPinteg[k]
        else:
            matrixPlist[i, k] = -(ep1 - ep2) * matrixPinteg[k]

for i in range(NumWLx):
    for k in range(qq):
        for j in range(qq):
            matrixRA[i, k, j] = matrixClist[i, k, j]
            matrixRA[i, k, j + qq] = matrixDlist[i, k, j]
            matrixRA[i, k + qq, j] = matrixFlist[i, k, j]
            matrixRA[i, k + qq, j + qq] = matrixGlist[i, k, j]

for i in range(NumWLx):
    for k in range(qq):
        matrixQA[i, k] = matrixElist[i, k]
        matrixQA[i, k + qq] = matrixHlist[i, k]

for i in range(NumWLx):
    for k in range(qq):
        for j in range(qq):
            matrixRB[i, k, j] = matrixJlist[i, k, j]
            matrixRB[i, k, j + qq] = matrixKlist[i, k, j]
            matrixRB[i, k + qq, j] = matrixMlist[i, k, j]
            matrixRB[i, k + qq, j + qq] = matrixNlist[i, k, j]

for i in range(NumWLx):
    for k in range(qq):
        matrixQB[i, k] = matrixLlist[i, k]
        matrixQB[i, k + qq] = matrixPlist[i, k]

for i in range(NumWLx):
    vectorA[i] = np.linalg.solve(matrixRA[i], matrixQA[i])
    vectorB[i] = np.linalg.solve(matrixRB[i], matrixQB[i])

    alpha_A[i] = -4 * pi * rr**3 * ep1 * vectorA[i, 0]
    alpha_B[i] = -4 * pi * rr**3 * ep1 * vectorB[i, 0]

alpha_A_Re = real(alpha_A)
alpha_B_Re = real(alpha_B)
alpha_A_Im = imag(alpha_A)
alpha_B_Im = imag(alpha_B)
alpha_A_Abs = abs(alpha_A)
alpha_B_Abs = abs(alpha_B)

Csca_A = k0**4 / (6 * pi) * alpha_A_Abs**2
Csca_B = k0**4 / (6 * pi) * alpha_B_Abs**2
Cabs_A = k0 * alpha_A_Im
Cabs_B = k0 * alpha_B_Im

Qsca_A = Csca_A / ((rr**2) * pi)
Qabs_A = Cabs_A / ((rr**2) * pi)
Qsca_B = Csca_B / ((rr**2) * pi)
Qabs_B = Cabs_B / ((rr**2) * pi)

plt.figure(figsize=(8, 6))
plt.plot(WLx, Qabs_A, label=r"$Q_{\rm abs}$", linewidth=3.0, color='black')
plt.plot(WLx, Qabs_B, label=r"$Q_{\rm abs}$", linewidth=3.0, color='gray')
plt.axis([400, 700, 0, 5])
plt.xlabel("wavelength (nm)", fontsize=22)
plt.ylabel("efficiency", fontsize=22)
plt.tick_params(labelsize=20)  # scale fontsize=18pt
plt.legend(fontsize=20, loc='upper left')
plt.tight_layout()
plt.show()
