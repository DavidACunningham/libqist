from itraj import Itraj_dense, Itraj_dtime, denseloads, Itraj, rfrloads
import matplotlib.pyplot as plt
import numpy as np
from dac3B import mu, LU, TU
from plotfuncs import gridplot, plotNRHO

qistd = Itraj_dense(**denseloads[1.5])
qistt = Itraj_dtime(**denseloads[0.0])
qistc = Itraj(**rfrloads[1.5])
qists = [ qistd, qistt]#, qistc]
qlabels = ["1.5", "0", "1.5_cheb"]
fig,ax = plt.subplots(constrained_layout=True)
for qist,n in zip(qists,qlabels):
    ts = np.linspace(qist.tau0, qist.tauf,1000)
    stmtestels = [qist.rstate(t)[0] for t in ts]
    ax.plot(ts/ts[-1],stmtestels,
            label=''.join(["$\\alpha=",n,"$"]))
ax.set_ylabel("$\phi^{23}$ (unitless)")
ax.set_xlabel("Fraction of Gateway Period")
ax.grid()

plt.legend()

# qist = Itraj(**rfrloads[1.5])
# ts = np.linspace(qist.tau0, qist.tauf, 1000)
# states = np.array([qist.rstate(t) for t in ts]).T

# plotNRHO(*states[:3])

plt.show()
