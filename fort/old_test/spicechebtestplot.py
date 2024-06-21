import numpy as np
import subprocess
import matplotlib.pyplot as plt

subprocess.run(["gfortran", "-g", "./spicechebytest.f90", "../src/cheby.f90", "/home/david/libf/spicelib.a"])
subprocess.run(["./a.out", "1", "50000", "6678000"])
testpoints = []
with open("testpoints.txt", 'r') as f:
    for line in f:
        testpoints.append(float(line))
truth = []
with open("truth.txt", 'r') as f:
    for line in f:
        truth.append(float(line))
dtruth = []
with open("dtruth.txt", 'r') as f:
    for line in f:
        dtruth.append(float(line))
cheb = []
with open("cheb.txt", 'r') as f:
    for line in f:
        cheb.append(float(line))
dcheb = []
with open("dcheb.txt", 'r') as f:
    for line in f:
        dcheb.append(float(line))

fig,ax = plt.subplots()
ax.plot(testpoints,np.array(truth),label="sin(x)")
ax.plot(testpoints,np.array(dtruth),label="cos(x)")
fig2,ax2 = plt.subplots()
ax2.plot(testpoints,np.abs(np.array(truth) - np.array(cheb)),label="sin error")
ax2.plot(testpoints,np.abs(np.array(dtruth) - np.array(dcheb)),label="deriv error")
ax.legend()
ax2.legend()
ax.set_xlabel("x")
ax2.set_xlabel("x")
ax2.set_ylabel("absolute error")
ax2.set_yscale("log")
plt.show()
