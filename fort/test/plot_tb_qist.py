import numpy as np
from addlib import *
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from coordrots import ijkfrompqw_els, kep2cart, cart2kep
from scipy.optimize import root
import datetime as dt

# ROUGH RULE OF THUMB: 100 degrees per month for moon
def plt_sphere(ax,c, r):
  # draw sphere
  u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
  x = r*np.cos(u)*np.sin(v)
  y = r*np.sin(u)*np.sin(v)
  z = r*np.cos(v)
  ax.plot_surface(x-c[0], y-c[1], z-c[2], color=np.random.choice(['g','b']), alpha=0.5*np.random.random()+0.5)
def kepler_prop_from_periapse(elements,dt, mu=1, fast_variable="true_anomaly",
                output_mode="els"):
    elements_out = elements.copy()
    n = np.sqrt(mu/elements[0]**3)
    M = n*dt
    kepfunc = lambda E: E - elements[1]*np.sin(E) - M
    dkepfunc = lambda E: 1 - elements[1]*np.cos(E)
    E = root(kepfunc,M,jac=dkepfunc,tol=1e-16)
    E = E.x
    factor = np.sqrt((1+elements[1])/(1-elements[1]))
    nu_f = 2*np.arctan2(factor*np.tan(E/2),1)
    elements_out[5] = nu_f
    if output_mode=="els":
        return elements_out
    elif output_mode=="pqw":
        return kep2cart(elements_out,mu=mu)
    elif output_mode=="ijk":
        return ijkfrompqw_els(elements_out)@kep2cart(elements_out,mu=mu)

posx = []
with open("twobody_qist_x.txt", 'r') as f:
    for line in f:
        posx.append(float(line))
posy = []
with open("twobody_qist_y.txt", 'r') as f:
    for line in f:
        posy.append(float(line))
posz = []
with open("twobody_qist_z.txt", 'r') as f:
    for line in f:
        posz.append(float(line))
time = []
with open("twobody_qist_t.txt", 'r') as f:
    for line in f:
        time.append(float(line))
xdot = []
with open("twobody_qist_xdot.txt", 'r') as f:
    for line in f:
        xdot.append(float(line))
ydot = []
with open("twobody_qist_ydot.txt", 'r') as f:
    for line in f:
        ydot.append(float(line))
zdot = []
with open("twobody_qist_zdot.txt", 'r') as f:
    for line in f:
        zdot.append(float(line))
<<<<<<< HEAD


spicex = []
with open("spice_resamp_x.txt", 'r') as f:
    for line in f:
        spicex.append(float(line))
spicey = []
with open("spice_resamp_y.txt", 'r') as f:
    for line in f:
        spicey.append(float(line))
spicez = []
with open("spice_resamp_z.txt", 'r') as f:
    for line in f:
        spicez.append(float(line))
ax.plot(spicex,spicey,spicez)
    if ind == 4:
        plota = []
        for item in (els[ind] - kepels[ind]):
            if item>np.pi:
                plota.append(item-2*np.pi)
            elif item < -np.pi:
                plota.append(item + 2*np.pi)
            else:
                plota.append(item)
        a.plot(time, plota)
    else:
        a.plot(time,els[ind]-kepels[ind])

plt.show()
breakpoint()
