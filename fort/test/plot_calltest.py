import numpy as np
from addlib import *
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from coordrots import ijkfrompqw_els, kep2cart
from scipy.optimize import root
import datetime as dt

# ROUGH RULE OF THUMB: 100 degrees per month for moon

def readFile(fname):
    # Reads in file that is just newline-
    # delimited
    arr = []
    with open(fname+".txt", 'r') as f:
        for line in f:
            arr.append(float(line))
    return np.array(arr)

body = "Test Spacecraft"

xpos_spice = readFile("spice_resamp_x")
ypos_spice = readFile("spice_resamp_y")
zpos_spice = readFile("spice_resamp_z")
xpos_tbtest = readFile("base_sol_x")
ypos_tbtest = readFile("base_sol_y")
zpos_tbtest = readFile("base_sol_z")
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.scatter(xpos_spice, ypos_spice, zpos_spice)
ax.scatter(xpos_tbtest, ypos_tbtest, zpos_tbtest)
ax.set_aspect("equal")
fig2, ax2 = plt.subplots()
ax2.plot(xpos_tbtest-xpos_spice)
ax2.plot(ypos_tbtest-ypos_spice)
ax2.plot(zpos_tbtest-zpos_spice)
plt.show()
