import numpy as np
import numpy.linalg as la
np.set_printoptions(linewidth=np.nan)
# these lines implement the operation
# ia, abc, jb -> ijc with RIU, HU, RIU

def fort_print_output(varname,out):
    if len(out.shape) == 1:
        print(f"{varname} = ", np.array2string(out, separator=f", "))
    if len(out.shape) == 2:
        for i, row in enumerate(out):
            print(f"{varname}({i+1},:) = ", np.array2string(row, separator=f", "))
    if len(out.shape) == 3:
        for i, page in enumerate(out):
            for j, row in enumerate(page):
                print(f"{varname}({i+1},{j+1},:) = ", np.array2string(row, separator=", "))
    print("")

mat = np.round(100*np.random.rand(3,3))
tens = np.round(100*np.random.rand(3,3,3))
vec = np.round(100*np.random.rand(3))
stm6 = np.round(100*np.random.rand(6,6))
stm8 = np.round(100*np.random.rand(8,8))
stm1 = np.round(100*np.random.rand(8,8))
stm2 = np.round(100*np.random.rand(8,8))
stt1 = np.round(100*np.random.rand(8,8,8))
stt2 = np.round(100*np.random.rand(8,8,8))

zmat = np.zeros((6,6))
zmat[:3,3:] = np.eye(3)
zmat[3:,:3] = -np.eye(3)
print("! vec")
fort_print_output("vec",vec)
print("! mat")
fort_print_output("mat",mat)
print("! tens")
fort_print_output("tens",tens)
print("!stm6")
fort_print_output("stm6",stm6)
print("!stm8")
fort_print_output("stm8",stm8)
print("!stm1")
fort_print_output("stm1",stm1)
print("!stm2")
fort_print_output("stm2",stm2)
print("!stt1")
fort_print_output("stt1",stt1)
print("!stt2")
fort_print_output("stt2",stt2)

print("! mattens")
output = np.einsum('ia, ajk -> ijk', mat, tens)
fort_print_output("mt_truth",output)

print("! quad")
output = np.einsum('iab, aj, bk -> ijk', tens, mat, mat)
fort_print_output("quad_truth",output)

print("! vectens1")
output = np.einsum('ajk, a -> jk', tens, vec)
fort_print_output("vectens1_truth",output)
print("! vectens2")
output = np.einsum('iak, a -> ik', tens, vec)
fort_print_output("vectens2_truth",output)
print("! vectens3")
output = np.einsum('ija, a -> ij', tens, vec)
fort_print_output("vectens3_truth",output)

print("! vectensquad")
output = np.einsum('iab, a, b -> i', tens, vec, vec)
fort_print_output("vectensquad_truth",output)

print("! STMinvert6")
output = -zmat@stm6.T@zmat
fort_print_output("stminvert6_truth",output)

print("! STMinvert8")
output = la.inv(stm8)
fort_print_output("stminvert8_truth",output)


