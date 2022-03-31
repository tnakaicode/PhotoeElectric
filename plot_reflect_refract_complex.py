import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
import glob
import imageio
import argparse

sys.path.append(os.path.join("./"))
from src.base import plot2d
from src.index import ref_index, ref_ratio, ref_ratio_complex

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--n12", dest="n12",
                      default=[1.5, 1.0], type=float, nargs=2)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = plot2d(aspect="auto")
    n1, n2 = opt.n12
    pt = np.linspace(0, 1, 100) * 90
    [tp, rp, ts, rs], [T_p, R_p, T_s, R_s] = ref_ratio_complex(n1=n1, n2=n2, t0=pt)

    title = "n1={:.2f} -> n2={:.2f}".format(n1, n2)
    obj.new_2Dfig(aspect="auto")
    obj.axs.set_title("Reflection and Transmition Index\n" + title)
    obj.axs.plot(pt, tp, label="P-pol Transmit")
    obj.axs.plot(pt, rp, label="P-pol Reflect")
    obj.axs.plot(pt, ts, label="S-pol Transmit")
    obj.axs.plot(pt, rs, label="S-pol Reflect")
    obj.axs.legend()
    obj.axs.set_xlim(0, 90)
    obj.axs.set_ylim(-1.1, 5.0)
    obj.SavePng(obj.tempname + "-index.png")

    obj.new_2Dfig(aspect="auto")
    obj.axs.set_title("Reflection and Transmition Ratio\n" + title)
    obj.axs.plot(pt, T_p, label="P-pol Transmit")
    obj.axs.plot(pt, R_p, label="P-pol Reflect")
    obj.axs.plot(pt, T_s, label="S-pol Transmit")
    obj.axs.plot(pt, R_s, label="S-pol Reflect")
    obj.axs.legend()
    obj.axs.set_xlim(0, 90)
    obj.axs.set_ylim(-0.1, 1.1)
    obj.SavePng(obj.tempname + "-ratio.png")
