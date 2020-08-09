import numpy as np


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

def ref_ratio (n1=1.0, n2=1.0, t0=0.0):
    t1 = np.deg2rad(t0)
    t2 = np.arcsin((n1 / n2) * np.sin(t1))
    tp, rp, ts, rs = ref_index(n1, n2, t0)
    R_p = rp**2
    T_p = tp**2*(n2*np.cos(t2)/ (n1*np.cos(t1)))
    R_s = rs**2
    T_s = ts**2*(n2*np.cos(t2)/ (n1*np.cos(t1)))
    return [tp, rp, ts, rs], [T_p, R_p, T_s, R_s]

