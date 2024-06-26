import numpy as np
from addlib import *
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from coordrots import ijkfrompqw_els, kep2cart
from scipy.optimize import root
from scipy.integrate import solve_ivp
import datetime as dt

def readFile(fname):
    # Reads in file that is just newline-
    # delimited
    arr = []
    with open(fname+".txt", 'r') as f:
        for line in f:
            arr.append(float(line))
    return np.array(arr)

time = readFile("./times_hes")
analytic = readFile("./hes477")
finite_diff = readFile("./hes477_fd")
difference = readFile("./hes477_diff")

plt.plot(time,abs(analytic), label='analytic')
plt.plot(time,abs(finite_diff), label='findiff')
plt.plot(time,abs(difference), label='difference')
plt.yscale("log")
plt.legend()
plt.show()
