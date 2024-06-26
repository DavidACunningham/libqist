import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
moon_mu = 0.490279980693169e4
moon_r = 1738.
energy = lambda x: np.sqrt(sum(x[3:6]**2))/2. - moon_mu/np.sqrt(sum(x[:3]**2))
efroma = lambda a: -moon_mu/(2.*a)
P = lambda a: 2.*np.pi*np.sqrt(a**3/moon_mu)
orbital_speed = lambda a,r: np.sqrt(moon_mu*(2./r - 1./a))
def rv(rmag,vmag):
    rhat = np.random.rand(3)
    vxvy = np.random.rand(2)
    rhat = rhat/np.sqrt(sum(rhat**2))
    r = rhat*rmag
    vz = (- (rhat[0]*vxvy[0] + rhat[1]*vxvy[1]))/rhat[2]
    vhat = np.array([vxvy[0], vxvy[1], vz])
    vhat = vhat/np.sqrt(sum(vhat**2))
    v = vhat*vmag
    print("rdotv: ", r@v)
    return (r,v)

rp = 1.1*moon_r
e = 0.3
a = rp/(1-e)
vp = orbital_speed(a,rp)
rpvec, vpvec = rv(rp,vp)
x0 = np.r_[rpvec, vpvec]
x0 = np.zeros(6)
x0[0] =  8.059849800040121E+02
x0[1] =  1.732107143792893E+03
x0[2] =  7.191866537665294E+01
x0[3] =  1.586030952010121E-01
x0[4] =  1.724388489487553E-03
x0[5] = -1.818978945873842E+00
np.set_printoptions(precision=15)
print("rpvec = ", rpvec)
print("vpvec = ", vpvec)
period = P(a)
print(period)

def eoms(t,x):
    return np.r_[x[3:6], -moon_mu*x[:3]/np.sqrt(sum(x[:3]**2))**3]

orbit = solve_ivp(eoms, [0, period], x0, rtol=1.e-12, atol=1.e-14)

fig = plt.figure()
ax  = fig.add_subplot(projection="3d")
ax.plot(*orbit.y[:3])
ax.set_aspect("equal")
ax.scatter(0,0,0)
plt.show()

