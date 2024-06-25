
import numpy as np
import matplotlib.pyplot as plt
def readFile(fname):
    # Reads in file that is just newline-
    # delimited
    arr = []
    with open(fname+".txt", 'r') as f:
        for line in f:
            arr.append(float(line))
    return np.array(arr)
