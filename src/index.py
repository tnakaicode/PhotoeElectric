import numpy as np


def mMATs(n1z, n2z):
    """
    S-Pol M-Matrix Boundary Condition
    """
    mat = np.empty((2, 2))
    mat[0, 0] = n1z + n2z
    mat[0, 1] = n1z - n2z
    mat[1, 0] = n1z - n2z
    mat[1, 1] = n1z + n2z
    return (1 / (2 * n1z)) * mat


def mMATp(n1z, n2z, n1, n2):
    """
    P-Pol M-Matrix Boundary Condition
    """
    mat = np.empty((2, 2))
    mat[0, 0] = n1**2 * n2z + n2**2 * n1z
    mat[0, 1] = n1**2 * n2z - n2**2 * n1z
    mat[1, 0] = n1**2 * n2z - n2**2 * n1z
    mat[1, 1] = n1**2 * n2z + n2**2 * n1z
    return (1 / (2 * n1 * n2 * n1z)) * mat


def matFAI(n1z, d1, k0):
    """
    Phi-Matrix Phase Trasition
    """
    mat = np.empty((2, 2))
    mat[0, 0] = np.exp(1j * n1z * k0 * d1)
    mat[0, 1] = 0
    mat[1, 0] = 0
    mat[1, 1] = np.exp(-1j * n1z * k0 * d1)
    return mat


def brewster_angle(n1=1.0, n2=1.0):
    return np.arctan(n2 / n1)


def ref_index(n1=1.0, n2=1.0, t0=0.0):
    t1 = np.deg2rad(t0)
    t2 = np.arcsin((n1 / n2) * np.sin(t1))
    cs = n1 * np.cos(t1) + n2 * np.cos(t2)
    ct = n2 * np.cos(t1) + n1 * np.cos(t2)
    rs = (n1 * np.cos(t1) - n2 * np.cos(t2)) / cs
    ts = 2 * n1 * np.cos(t1) / cs
    rp = (n2 * np.cos(t1) - n1 * np.cos(t2)) / ct
    tp = 2 * n1 * np.cos(t1) / ct
    return tp, rp, ts, rs


def ref_ratio(n1=1.0, n2=1.0, t0=0.0):
    t1 = np.deg2rad(t0)
    t2 = np.arcsin((n1 / n2) * np.sin(t1))
    tp, rp, ts, rs = ref_index(n1, n2, t0)
    R_p = rp**2
    T_p = tp**2 * (n2 * np.cos(t2) / (n1 * np.cos(t1)))
    R_s = rs**2
    T_s = ts**2 * (n2 * np.cos(t2) / (n1 * np.cos(t1)))
    return [tp, rp, ts, rs], [T_p, R_p, T_s, R_s]


def ref_ratio_complex(n1=1.0, n2=1.0, t0=0.0):
    t1 = np.deg2rad(t0)
    s1 = np.sin(t1)
    c1 = np.cos(t1)
    e1 = n1**2
    s2 = n1 / n2 * s1
    c2 = np.sqrt(1 - s2**2, dtype=complex)
    e2 = n2**2
    cs = n1 * c1 + n2 * c2
    ct = n2 * c1 + n1 * c2

    rs = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
    rp = (e2 * n1 * c1 - e1 * n2 * c2) / (e2 * n1 * c1 + e1 * n2 * c2)
    ts = 2 * n1 * c1 / cs
    tp = 2 * n1 * c1 / ct
    R_p = np.abs(rp)**2
    T_p = np.abs(tp)**2 * (n2 * c2 / (n1 * c1))
    R_s = np.abs(rs)**2
    T_s = np.abs(ts)**2 * (n2 * c2 / (n1 * c1))
    return [tp, rp, ts, rs], [T_p, R_p, T_s, R_s]
