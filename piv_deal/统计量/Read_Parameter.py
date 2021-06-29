import os
import numpy as np


def read_p(para_metre, nu):
    mul = np.zeros(nu)
    f = open(para_metre, 'r')
    for i in range(nu):
        line1 = f.readline()
        sdy = line1.split("=")
        mul[i] = float(sdy[1])
    f.close()
    return mul
