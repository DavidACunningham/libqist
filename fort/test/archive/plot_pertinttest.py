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
with open("twobody_integrated_t.txt", 'r') as f:
    for line in f:
        times.append(float(line))
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
intposx = []
with open("twobody_integrated_x.txt", 'r') as f:
    for line in f:
        intposx.append(float(line))
intposy = []
with open("twobody_integrated_y.txt", 'r') as f:
    for line in f:
        intposy.append(float(line))
intposz = []
with open("twobody_integrated_z.txt", 'r') as f:
    for line in f:
        intposz.append(float(line))
posx = np.array(posx)
posy = np.array(posy)
posz = np.array(posz)
intposx = np.array(intposx)
intposy = np.array(intposy)
intposz = np.array(intposz)
diffs = [posx - intposx, posy - intposy, posz - intposz]
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.scatter(posx,posy,posz)
# ax.scatter(intposx,intposy,intposz,label="integrated")
ax.legend()
ax.set_aspect("equal")
fig2, ax2 = plt.subplots()
for item in diffs:
    ax2.plot(times, item)
plt.show()
