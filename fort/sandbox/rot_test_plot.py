import numpy as np
import spiceypy as spy
import numpy.linalg as la
import matplotlib.pyplot as plt

def readFile(fname):
    # Reads in file that is just newline-
    # delimited
    arr = []
    with open(fname+".txt", 'r') as f:
        for line in f:
            arr.append(float(line))
    return np.array(arr)

q0 = readFile("q0")
q1 = readFile("q1")
q2 = readFile("q2")
q3 = readFile("q3")
qdot0 = readFile("qdot0")
qdot1 = readFile("qdot1")
qdot2 = readFile("qdot2")
qdot3 = readFile("qdot3")
qdot_fd0 = readFile("qdotfd0")
qdot_fd1 = readFile("qdotfd1")
qdot_fd2 = readFile("qdotfd2")
qdot_fd3 = readFile("qdotfd3")
times = readFile("testtimes")

fig,ax = plt.subplots(4,1)
fig2,ax2 = plt.subplots(4,1)
fig3,ax3 = plt.subplots(4,1)
fig4,ax4 = plt.subplots(4,1)

ax[0].plot(times,qdot0)
ax[1].plot(times,qdot1)
ax[2].plot(times,qdot2)
ax[3].plot(times,qdot3)
ax2[0].plot(times,qdot_fd0)
ax2[1].plot(times,qdot_fd1)
ax2[2].plot(times,qdot_fd2)
ax2[3].plot(times,qdot_fd3)
ax3[0].plot(times,q0)
ax3[1].plot(times,q1)
ax3[2].plot(times,q2)
ax3[3].plot(times,q3)
ax4[0].plot(times,qdot0 - qdot_fd0)
ax4[1].plot(times,qdot1 - qdot_fd1)
ax4[2].plot(times,qdot2 - qdot_fd2)
ax4[3].plot(times,qdot3 - qdot_fd3)
fig.suptitle("Analytic")
fig2.suptitle("FD")
fig3.suptitle("q")
fig4.suptitle("deriv diff")
plt.show()
