import numpy as np
import subprocess
import matplotlib.pyplot as plt
import datetime as dt

bod_dict = {"Mercury" : "1",
            "Venus"   : "2",
            "Moon"    : "301",
            "Mars"    : "4",
            "Jupiter" : "5",
            "Saturn"  : "6",
            "Uranus"  : "7",
            "Neptune" : "8",
            "Pluto"   : "9",
            }


epoch = dt.date(2000,1,1)
start_date = dt.date(2023,10,10)
stop_date = dt.date(2024,4,10)
start_sec = (start_date - epoch).total_seconds()
stop_sec = (stop_date - epoch).total_seconds()
body = "Moon"
# ROUGH RULE OF THUMB: 100 degrees per month for moon
degree = 200

subprocess.run(["gfortran", "-g", "./resamptest.f90", "../src/cheby.f90", "/home/david/libf/spicelib.a"])
subprocess.run(["./a.out", bod_dict[body], start_sec.__repr__(), stop_sec.__repr__(), degree.__repr__()])
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
ddcheb = []
with open("ddcheb.txt", 'r') as f:
    for line in f:
        ddcheb.append(float(line))

fig,ax = plt.subplots(3,1,sharex='col')
ax[0].plot(testpoints,np.array(truth),label="pos")
ax[0].plot(testpoints,np.array(cheb),label="pos")
ax[0].set_title("pos (km)")
ax[1].plot(testpoints,np.array(dtruth),label="vel")
ax[1].plot(testpoints,np.array(dcheb),label="vel")
ax[1].set_title("vel (km/s)")
ax[2].plot(testpoints,np.array(ddcheb),label="acc")
ax[2].set_title("acc (km/s^2)")
fig2,ax2 = plt.subplots()
ax2.plot(testpoints,np.abs(np.array(truth) - np.array(cheb)),label="pos error")
ax2.plot(testpoints,np.abs(np.array(dtruth) - np.array(dcheb)),label="vel error")
ax2.legend()
ax[2].set_xlabel("time (s)")
ax2.set_xlabel("time (s)")
ax2.set_ylabel("absolute error (km or km/s)")
ax2.set_yscale("log")
fig.suptitle(body+" J2000 x-coord, relative to Earth")
plt.show()
