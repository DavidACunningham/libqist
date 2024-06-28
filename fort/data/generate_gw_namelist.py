import spiceypy as spy
t0_utc = "2024 Jun 12 00:00:00"
tf_utc = "2024 Jun 14 00:00:00"
t0_utc_stripped = t0_utc.replace(" ","").replace(":","")
tf_utc_stripped = tf_utc.replace(" ","").replace(":","")
t0tf = t0_utc_stripped + tf_utc_stripped
spy.furnsh("../../../kernels/mk_gw.tf")
t0_tdb = spy.utc2et(t0_utc)
tf_tdb = spy.utc2et(tf_utc)
absolute_qist_path = "/home/david/wrk/nstgro/qist"

replace_dict = {
                "DATES"          : t0tf,
                "INIT_TIME_TDB"  : str(t0_tdb),
                "INIT_TIME_UTC"  : t0_utc,
                "FINAL_TIME_TDB" : str(tf_tdb),
                "FINAL_TIME_UTC" : tf_utc,
                "QISTPATH"       : absolute_qist_path,
               }
                

print("Writing Gateway config file from " + t0_utc + " to " + tf_utc)
with open("./gw_config_namelists_TEMPLATE.nml",'r') as f:
    template_data = f.read()

formatted_data = template_data

for k,v in replace_dict.items():
    formatted_data= formatted_data.replace(k,v)

with open(f"./gw_config_namelists_{t0tf}.nml","w") as f:
    f.write(formatted_data)
