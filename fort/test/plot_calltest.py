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

body = "Test Spacecraft"
times = []
posx = []
with open("testx.txt", 'r') as f:
    for line in f:
        posx.append(float(line))
posy = []
with open("testy.txt", 'r') as f:
    for line in f:
        posy.append(float(line))
posz = []
with open("testz.txt", 'r') as f:
    for line in f:
        posz.append(float(line))

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.scatter(posx,posy,posz)
ax.set_aspect("equal")
plt.show()
