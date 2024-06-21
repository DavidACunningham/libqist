import numpy as np
import subprocess
import matplotlib.pyplot as plt

subprocess.run(["gfortran", "chebytest.f90", "../src/cheby.f90"])
subprocess.run(["./a.out"])
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
ddtruth = []
with open("ddtruth.txt", 'r') as f:
    for line in f:
        ddtruth.append(float(line))
cheb = []
with open("cheb.txt", 'r') as f:
    for line in f:
        cheb.append(float(line))
dcheb = []
with open("dcheb.txt", 'r') as f:
    for line in f:
        dcheb.append(float(line))
ddcheb = []
with open("ddcheb.txt", 'r') as f:
    for line in f:
        ddcheb.append(float(line))
coeffs = []
with open("coeffs.txt", 'r') as f:
    for line in f:
        coeffs.append(float(line))
dcoeffs = []
with open("dcoeffs.txt", 'r') as f:
    for line in f:
        dcoeffs.append(float(line))

ncoeffs = np.array(range(len(coeffs)))
ndcoeffs = np.array(range(len(dcoeffs)))
coeffs = np.array(coeffs)
dcoeffs = np.array(dcoeffs)

fig,ax = plt.subplots()
ax.plot(testpoints,np.array(truth),label="sin(x)")
ax.plot(testpoints,np.array(dtruth),label="cos(x)")
fig2,ax2 = plt.subplots()
ax2.plot(testpoints,np.abs(np.array(truth) - np.array(cheb)),label="sin error")
ax2.plot(testpoints,np.abs(np.array(dtruth) - np.array(dcheb)),label="deriv error")
ax2.plot(testpoints,np.abs(np.array(ddtruth) - np.array(ddcheb)),label="2nd deriv error")
ax.legend()
ax2.legend()
ax.set_xlabel("x")
ax2.set_xlabel("x")
ax2.set_ylabel("absolute error")
ax2.set_yscale("log")
fig3,ax3 = plt.subplots()
ax3.scatter(ncoeffs[coeffs>1.e-14],coeffs[coeffs>1.e-14],label="Sin Coeffs",marker="x")
ax3.scatter(ndcoeffs[dcoeffs>1.e-14],dcoeffs[dcoeffs>1.e-14],label="Cos Coeffs",marker="x")
ax3.set_yscale("log")
ax3.set_xticks(ncoeffs[::2])
ax3.legend()
plt.show()
