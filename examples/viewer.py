import matplotlib.pyplot as plt
import numpy as np
import spiceypy as spy
from itraj import Itraj
def makesphere(ax,x0=0,y0=0,z0=0,r=1.,color=None):
    # make data
    rx = 8.04
    ry = 5.89
    rz = 5.11
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = rx*np.outer(np.cos(u), np.sin(v)) + x0
    y = ry*np.outer(np.sin(u), np.sin(v)) + y0
    z = rz*np.outer(np.ones(np.size(u)), np.cos(v)) + z0

    # plot the surface
    if color is None:
        color='gray'
    ax.plot_surface(x, y, z, color=color)
mu_d = 0.0000985
t0_utc = "2026 Nov 26 12:00:00.00"
tf_utc = "2026 Nov 26 17:45:00.00"
spy.furnsh('../kernels/mk_example_with_traj.tf')

t0 = spy.utc2et(t0_utc)
tf = spy.utc2et(tf_utc)

times = np.linspace(t0,tf,1000)

states = np.array([spy.spkgeo(-31415, tt, "MDROTBAR", 402)[0] for tt in times])

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(*states.T[:3])
makesphere(ax)
ax.set_aspect('equal')

it = Itraj("./curve_deimos_config_namelist_2026Nov26120000002026Nov2617450000.nml")
dx0 = np.array([0, 0.1, 0., 0.1*np.sqrt(mu_d/10**3)/2, 0., 0.])
dxbs = np.array([it.prop(t0, tt, np.r_[dx0,0.,0]) for tt in times])
fig2 = plt.figure()
ax2 = fig2.add_subplot(111,projection='3d')
ax2.plot(*dxbs.T[:3])
ax2.set_aspect('equal')

plt.show()
