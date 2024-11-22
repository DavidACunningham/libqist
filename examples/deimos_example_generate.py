import spiceypy as spy
import subprocess
from time import time
from math import factorial as ff
from math import sqrt
kd = lambda m: {0:1}.get(m,0)
def knorm(l,m):
    return sqrt(
        ff(l+m)/ff(l-m) / ((2-kd(m))*(2*l+1))
            )
def get_shdegree(c,s):
    deg = lambda x: x[0]
    order = lambda x: x[1]
    maxes = [
             max(c,key=deg)[0],
             max(c,key=order)[1],
             max(s,key=deg)[0],
             max(s,key=order)[1],
            ]
    return str(max(maxes))


itime = time()
t0_utc = "2026 Nov 26 12:00:00.00"
tf_utc = "2026 Nov 26 17:45:00.00"
rtol_qist = 1.e-12
atol_qist = 1.e-20

rtol_kernel = 1.e-16
atol_kernel = 1.e-20

kernel_nodes = 400

trajectory_name = "curve_deimos"
trajectory_naif_id = -31415
absolute_qist_path = "/home/david/wrk/nstgro/libqist/"
kernel_path = "kernels/"
datadir = "datadir/"
relative_data_path = "examples/"
central_body_naif_id = 402
central_body_mu = 9.85000000E-05
central_body_reference_radius = 6.200000E+00

regularize = True
resample_interpolation_degree = 200
rotation_interpolation_degree = 200
rotating_frame = "MDROTBAR"
inertial_frame = "J2000"

perturbing_body_naif_ids =  [
                499, # Mars
                401, # Phobos
                10, # Sun
                5, # Jupiter Barycenter
                  ]
perturbing_body_mus   =  [ 
                    42828.3, # Mars
                    0.0007087, # Phobos
                    132712440041.939400, # Sun
                    126712764.800000, # Jupiter Barycenter
                   ]
init_state_array = [
        0.,
        3.42020143,
        9.39692621,
        0.003138471,
        0.,
        0.,
        ] 
cdict = {
    (0,0) : 1.0000000000000000E+00,
    (1,0) : 0.0000000000000000E+00, 
    (1,1) : 0.0000000000000000E+00, 
    (2,0) :-1.0790000000000000E-01,
    (2,1) :-9.3900000000000000E-04,
    (2,2) : 3.0810000000000000E-02,
    (3,0) : 2.5660000000000000E-02,
    (3,1) : 1.5450000000000000E-02,
    (3,2) :-3.7900000000000000E-03,
    (3,3) :-1.6800000000000000E-04,
           }
sdict = {
    (1,0) : 0.0000000000000000E+00, 
    (1,1) : 0.0000000000000000E+00, 
    (2,1) : 4.6800000000000000E-03, 
    (2,2) : 7.9100000000000000E-04, 
    (3,1) : 5.6200000000000000E-04,
    (3,2) :-1.7000000000000000E-03,
    (3,3) : 1.0100000000000000E-03,
        }
# ~~~~Edit between these lines

cbardict = {k:knorm(*k)*v for k,v in cdict.items()}
sbardict = {k:knorm(*k)*v for k,v in sdict.items()}


kernelpath = absolute_qist_path+kernel_path
path = absolute_qist_path + relative_data_path + datadir
spy.furnsh(absolute_qist_path+kernel_path+"mk_example.tf")
t0_utc_stripped = t0_utc.replace(" ","").replace(":","").replace(".","")
tf_utc_stripped = tf_utc.replace(" ","").replace(":","").replace(".","")
t0tf = t0_utc_stripped + tf_utc_stripped
qistname = path+trajectory_name + t0tf + ".qist"
kernelname = absolute_qist_path + kernel_path + "example_binary_kernel.bsp"
metakernel_filename_no_trajectory = kernelpath + "mk_example.tf"
metakernel_filename_with_trajectory = kernelpath +  "mk_example_with_traj.tf"
resample_filename_no_trajectory = path + "example_resample_no_traj_"+t0tf+".resamp"
resample_filename_with_trajectory = path + "example_resample_with_traj_"+t0tf+".resamp"
rot_filename = path + trajectory_name + rotating_frame +"_"+ t0tf + ".rot"
if regularize:
    regname = path + trajectory_name + t0tf + "_kvtau.odesolution"
else:
    regname = ""
t0_tdb = spy.utc2et(t0_utc)
tf_tdb = spy.utc2et(tf_utc)

qs = '"'
replace_dict = {
                "DATES"                   : t0tf,
                "QIST_T0_TDB"             : f"{t0_tdb:.15E}",
                "QIST_T0_UTC"             : t0_utc,
                "QIST_TF_TDB"             : f"{tf_tdb:.15E}",
                "QIST_TF_UTC"             : tf_utc,
                "RESAMP_T0_TDB"           : f"{t0_tdb:.15E}",
                "RESAMP_T0_UTC"           : t0_utc,
                "RESAMP_TF_TDB"           : f"{tf_tdb:.15E}",
                "RESAMP_TF_UTC"           : tf_utc,
                "QISTFILE"                : qs+qistname+qs,
                "SH_DEGREE"               : get_shdegree(cbardict,sbardict),
                "REGFILE"                 : qs+regname+qs,
                "KERNEL_FILE"             : qs+kernelname+qs,
                "MK_WITH_TRAJECTORY"      : qs+metakernel_filename_with_trajectory+qs,
                "MK_NO_TRAJECTORY"        : qs+metakernel_filename_no_trajectory+qs,
                "RESAMP_WITH_TRAJECTORY"  : qs+resample_filename_with_trajectory+qs,
                "RESAMP_NO_TRAJECTORY"    : qs+resample_filename_no_trajectory+qs,
                "RESAMP_DEG"              : str(resample_interpolation_degree),
                "ROT_FILE"                : qs+rot_filename+qs,
                "ROT_DEG"                 : str(rotation_interpolation_degree),
                "INERTIAL_FRAME"          : qs+inertial_frame+qs,
                "ROTATING_FRAME"          : qs+rotating_frame+qs,
                "TRAJNAME"                : trajectory_name,
                "RTOL_QIST"               : f"{rtol_qist:.15E}",
                "ATOL_QIST"               : f"{atol_qist:.15E}",
                "REGULARIZE"              : f".{regularize}.".lower(),
                "CENTRAL_BODY_ID"         : str(central_body_naif_id),
                "CENTRAL_BODY_MU"         : f"{central_body_mu:.15E}",
                "CENTRAL_BODY_REF_RADIUS" : f"{central_body_reference_radius}",
                "TRAJ_ID"                 : str(trajectory_naif_id),
                "NUM_BODIES"              : str(len(perturbing_body_naif_ids)),
                "KERNEL_NODES"            : str(kernel_nodes),
                "RTOL_KERNEL"             : f"{rtol_kernel:.15E}",
                "ATOL_KERNEL"             : f"{atol_kernel:.15E}",
               }

newname = absolute_qist_path \
        + relative_data_path \
        + f"{trajectory_name}_config_namelist_{t0tf}.nml"
mu_list ={i+1:v for i,v in enumerate(perturbing_body_mus)}
body_list ={i+1:v for i,v in enumerate(perturbing_body_naif_ids)}
init_state = {i+1:v for i,v in enumerate(init_state_array)}
array_replace_list = { "BODY_LIST"  : ("body_list",body_list),
                       "MU_LIST"    : ("mu_list", mu_list), 
                       "INIT_STATE" : ("x0", init_state), 
                       "CBAR"       : ("cbar", cbardict), 
                       "SBAR"       : ("sbar", sbardict),
                      }

def write_array(handle, arrayname, arraydict, startind=1):
    for k,v in arraydict.items():
        if type(v) is int:
            handle.write(f"     {arrayname}({k})            = {v}\n")
        else:
            if type(k) is tuple:
                handle.write(f"     {arrayname}{k}            = {v:.15E}\n")
            else:
                handle.write(f"     {arrayname}({k})            = {v:.15E}\n")
def parseline(handle, line):
    replaceline = False
    for k,v in array_replace_list.items():
        if k in line:
            replaceline = True
            write_array(new,*v)
    if not replaceline:
        newline = line
        for k,v in replace_dict.items():
            newline = newline.replace(k,v)
        new.write(newline)
print("Writing " + trajectory_name  +"  config file from " + t0_utc + " to " + tf_utc)
with open("./config_namelist_template.nml",'r') as old, open(newname,"w") as new:
    for line in old:
        parseline(new,line)
print("GENERATING NEW KERNEL")
subprocess.run(["../fort/exe/make_new_kernel", newname, "True", "True"])
print("GENERATING NEW QIST MODEL")
subprocess.run(["../fort/exe/make_qist_existing_kernel", newname, "True", "True"])
ftime = time()
print(f"Done. . . {ftime-itime} s elapsed")
