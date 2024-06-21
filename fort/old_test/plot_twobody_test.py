import numpy as np
from addlib import *
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from coordrots import ijkfrompqw_els, kep2cart
from scipy.optimize import root
import datetime as dt

# ROUGH RULE OF THUMB: 100 degrees per month for moon

subprocess.run(["gfortran", "-g", "./twobody_test.f90", "../src/cheby.f90", "/home/david/libf/spicelib.a"])
subprocess.run(["./a.out"])
body = "Test Spacecraft"
times = []
with open("times.txt", 'r') as f:
    for line in f:
        times.append(float(line))
spiceposx = []
with open("spiceposx.txt", 'r') as f:
    for line in f:
        spiceposx.append(float(line))
spiceposy = []
with open("spiceposy.txt", 'r') as f:
    for line in f:
        spiceposy.append(float(line))
spiceposz = []
with open("spiceposz.txt", 'r') as f:
    for line in f:
        spiceposz.append(float(line))

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
mu = 3.986e5
inc = 70*np.pi/180
x0 = np.array([16000., 5.e-1, inc, np.pi/2, 0, 0]) 
x = []
for t in times:
    x.append(kepler_prop_from_periapse(x0,t,mu=mu,output_mode='ijk'))
truth = np.array(x)
fig,ax = plt.subplots(3,1,sharex='col')
ax[0].plot(times,truth.T[0],label="pos")
ax[0].plot(times,np.array(spiceposx),label="pos")
ax[0].set_title("x (km)")
ax[1].plot(times,truth.T[1],label="vel")
ax[1].plot(times,np.array(spiceposy),label="vel")
ax[1].set_title("y (km)")
ax[2].plot(times,np.array(spiceposz),label="acc")
ax[2].plot(times,truth.T[2],label="acc")
ax[2].set_title("z (km)")
fig2,ax2 = plt.subplots()
spicepos = np.array([spiceposx, spiceposy, spiceposz])
ax2.plot(times,la.norm(truth.T[:3]- spicepos,axis=0))
ax[2].set_xlabel("time (s)")
ax2.set_xlabel("time (s)")
ax2.set_ylabel("absolute error (km)")
ax2.set_yscale("log")
fig.suptitle(body+" J2000 coords, relative to Earth")
plt.show()
