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

# ROUGH RULE OF THUMB: 100 degrees per month for moon
def kepler_prop_from_periapse(elements,dt, mu=1, fast_variable="true_anomaly",
                output_mode="els"):
    elements_out = elements.copy()
    a = elements[0]
    ecc = elements[1]
    n = np.sqrt(mu/a**3)
    M = n*dt
    kepfunc = lambda x: x - ecc*np.sin(x) - M
    dkepfunc = lambda x: 1 - ecc*np.cos(x)
    E = root(kepfunc,M,jac=dkepfunc,tol=1e-16)
    E = E.x[0]
    nu = 2.*np.arctan2(np.sqrt((1+ecc)/(1-ecc))*np.tan(E/2.),1.)
    elements_out[5] = nu
    if output_mode=="els":
        return elements_out
    elif output_mode=="pqw":
        return kep2cart(elements_out,mu=mu)
    elif output_mode=="ijk":
        return ijkfrompqw_els(elements)@kep2cart(elements_out,mu=mu)
mu = 398600.5
inc = 70.*np.pi/180.
x0 = np.array([16000., 5.e-1, inc, np.pi/2., 0., 0.]) 
print("Initial State:")
x0_cart = ijkfrompqw_els(x0)@kep2cart(x0,mu=mu)
for xx in x0_cart:
    print(xx)
x = []

def readFile(fname):
    # Reads in file that is just newline-
    # delimited
    arr = []
    with open(fname+".txt", 'r') as f:
        for line in f:
            arr.append(float(line))
    return np.array(arr)

body = "Test Spacecraft"

time = readFile("./times_resamp")
xpos_spice = readFile("spice_resamp_x")
ypos_spice = readFile("spice_resamp_y")
zpos_spice = readFile("spice_resamp_z")
xpos_preresamp = readFile("integrate_pre_resamp_x")
ypos_preresamp = readFile("integrate_pre_resamp_y")
zpos_preresamp = readFile("integrate_pre_resamp_z")
for t in time:
    x.append(kepler_prop_from_periapse(x0,t,mu=mu,output_mode='ijk'))
x = np.array(x)
x_unperturbed = x.T[0]
y_unperturbed = x.T[1]
z_unperturbed = x.T[2]
def xdot_kep(t,y):
    y1  =  mu/np.sqrt(sum(y[:3]**2))**3.
    res = np.r_[y[3:6], -y1*y[:3]]
    return res
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.scatter(xpos_spice, ypos_spice, zpos_spice)
ax.scatter(x_unperturbed, y_unperturbed, z_unperturbed)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_aspect("equal")
fig2, ax2 = plt.subplots(nrows=4, ncols=2)
cut = 1000
ax2[0][0].plot(time[:cut],xpos_spice[:cut],label="spice")
ax2[1][0].plot(time[:cut],ypos_spice[:cut],label="spice")
ax2[2][0].plot(time[:cut],zpos_spice[:cut],label="spice")
ax2[0][1].plot(time[:cut],xpos_preresamp[:cut]-x_unperturbed[:cut],label="pre-resamp")
ax2[1][1].plot(time[:cut],ypos_preresamp[:cut]-y_unperturbed[:cut],label="pre-resamp")
ax2[2][1].plot(time[:cut],zpos_preresamp[:cut]-z_unperturbed[:cut],label="pre-resamp")
ax2[0][1].plot(time[:cut],xpos_spice[:cut]-x_unperturbed[:cut],label="spice")
ax2[1][1].plot(time[:cut],ypos_spice[:cut]-y_unperturbed[:cut],label="spice")
ax2[2][1].plot(time[:cut],zpos_spice[:cut]-z_unperturbed[:cut],label="spice")
# for ax in ax2.T[1]:
#     ax.set_yscale("log")
for axrow in ax2[:-1]:
    for ax in axrow:
        ax.legend()
plt.show()
