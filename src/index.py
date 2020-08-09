import numpy as np


def ref_index(n1, n2, t0):
    t1 = np.deg2rad(t0)
    t2 = np.arcsin((n1 / n2) * np.sin(t1))
    rs = (n1 * np.cos(t1) - n2 * np.cos(t2)) / \
        (n1 * np.cos(t1) + n2 * np.cos(t2))
    ts = 2 * n1 * np.cos(t1) / (n1 * np.cos(t1) + n2 * np.cos(t2))
    rp = (n2 * np.cos(t1) - n1 * np.cos(t2)) / \
        (n2 * np.cos(t1) + n1 * np.cos(t2))
    tp = 2 * n1 * np.cos(t1) / (n2 * np.cos(t1) + n1 * np.cos(t2))
    return tp, rp, ts, rs
