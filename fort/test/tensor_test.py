import numpy as np
# these lines implement the operation
# ia, abc, jb -> ijc with RIU, HU, RIU
mat = np.zeros((3,3))
tens = np.zeros((3,3,3))
mat[0,:] = np.array([1., 2., 3.])
mat[1,:] = np.array([4., 5., 6.])
mat[2,:] = np.array([7., 8., 9.])

tens[0,0,:] = np.array([10., 19., 28.])
tens[0,1,:] = np.array([11., 20., 29.])
tens[0,2,:] = np.array([12., 21., 30.])
tens[1,0,:] = np.array([13., 22., 31.])
tens[1,1,:] = np.array([14., 23., 32.])
tens[1,2,:] = np.array([15., 24., 33.])
tens[2,0,:] = np.array([16., 25., 34.])
tens[2,1,:] = np.array([17., 26., 35.])
tens[2,2,:] = np.array([18., 27., 36.])

vec = np.array([37., 38., 39.])
output = np.einsum('ia, abc, jb, c -> ij',mat, tens, mat, vec)

for row in output:
    print(row)

