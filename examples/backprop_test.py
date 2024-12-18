import matplotlib.pyplot as plt
import numpy as np
import spiceypy as spy
from itraj import Itraj
mu_d = 0.0000985
t0_utc = "2026 Nov 26 12:00:00.00"
tf_utc = "2026 Nov 26 17:45:00.00"
spy.furnsh('../kernels/mk_example_with_traj.tf')

t0 = spy.utc2et(t0_utc)
tf = spy.utc2et(tf_utc)

times = np.linspace(t0,tf,1000)

states = np.array([spy.spkgeo(-31415, tt, "MDROTBAR", 402)[0] for tt in times])

it = Itraj("./curve_deimos_config_namelist_2026Nov26120000002026Nov2617450000.nml")
dxa = np.array([0, 0.1, 0., 0.1*np.sqrt(mu_d/10**3)/2, 0., 0.,0.,0.])
forward_prop = np.array([it.prop(t0,t,dxa) for t in times])
back_prop_from_end = np.array([it.prop_back(tf, t, forward_prop[-1]) for t in times])
back_prop_by_step = [forward_prop[-1]]
for idx in range(len(times)-1,0,-1):
    back_prop_by_step.append(it.prop_back(times[idx],times[idx-1],back_prop_by_step[-1]))

back_prop_by_step = np.array(back_prop_by_step[::-1])
re_prop_error_from_end = np.array([it.prop(t,tf,bp) - forward_prop[-1] for t, bp in zip(times,back_prop_from_end)])
re_prop_error_by_step = np.array([it.prop(t,tf,bp) - forward_prop[-1] for t, bp in zip(times,back_prop_by_step)])

forward_prop = forward_prop.T
back_prop_from_end=back_prop_from_end.T
back_prop_by_step = back_prop_by_step.T
re_prop_error_from_end = re_prop_error_from_end.T
re_prop_error_by_step = re_prop_error_by_step.T
fig,ax = plt.subplots(6,1,sharex=True)
for a,b1,b2 in zip(ax, re_prop_error_from_end,re_prop_error_by_step):
    a.plot(times,abs(b1),ls=":")
    a.plot(times,abs(b2),ls="--")
    a.set_yscale('log')
# for a,f,b1,b2 in zip(ax,forward_prop, back_prop_from_end,prop_back_by_step):
#     a.plot(times,f)
#     a.plot(times,b1,ls=":")
#     a.plot(times,b2,ls="--")
np.set_printoptions(precision=8, linewidth=np.nan)
plt.show()
