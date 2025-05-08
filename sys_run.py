import os
import re
import gc

import ROOT
ROOT.gROOT.SetBatch(True)

def clear_vars():
    for name in dir():
        if not name.startswith('_'):
            del globals()[name]
    gc.collect()

#make configuration file
def make_track_config(file_name,pars_dict):

    os.system("cp Configurations.py Sys_config/"+file_name+".py -rf")

    for key in pars_dict:
        key_re = re.compile(rf"{chr(34)}{key}{chr(34)}\s*:\s*.*?,")

        config_file = "Sys_config/" + file_name + ".py"
        config = open(config_file,"r+")
        config_content = config.read()

        try:
            key_value = key_re.search(config_content).group()[:-1]

        except:
            print(f"Key {key} not found in the configuration file")
            continue

        new_content = config_content.replace(key_value,f'"{key}": {pars_dict[key]}')

        config.seek(0)
        config.write(new_content)
        config.truncate()
        config.close()

            

if __name__ == "__main__":

    # make_track_config("config_sys_track_0",{"Ana_name": "\"SysAna_track_0\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # make_track_config("config_sys_track_1",{"Ana_name": "\"SysAna_track_1\"","Min_cls_ITS": 1, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # make_track_config("config_sys_track_2",{"Ana_name": "\"SysAna_track_2\"","Min_cls_ITS": 2, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # make_track_config("config_sys_track_3",{"Ana_name": "\"SysAna_track_3\"","Min_cls_ITS": 3, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # make_track_config("config_sys_track_4",{"Ana_name": "\"SysAna_track_4\"","Min_cls_ITS": 4, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # make_track_config("config_sys_track_5",{"Ana_name": "\"SysAna_track_5\"","Min_cls_ITS": 0, "Min_cls_TPC": 1, "Min_eta_track": 0})
    # make_track_config("config_sys_track_6",{"Ana_name": "\"SysAna_track_6\"","Min_cls_ITS": 0, "Min_cls_TPC": 2, "Min_eta_track": 0})
    # make_track_config("config_sys_track_7",{"Ana_name": "\"SysAna_track_7\"","Min_cls_ITS": 0, "Min_cls_TPC": 3, "Min_eta_track": 0})
    # make_track_config("config_sys_track_8",{"Ana_name": "\"SysAna_track_8\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 1})
    # make_track_config("config_sys_track_9",{"Ana_name": "\"SysAna_track_9\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 2})
    # make_track_config("config_sys_track_10",{"Ana_name": "\"SysAna_track_10\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 3})

    outstd_dir = "/home/mingze/work/dstar/Dstar_Spin_Alignment/Output/SysUncer/run_outstd"

    run_track_sys = True
    # run_track_sys = False

    # run_fit_sys = True
    run_fit_sys = False

    # run_cut_variation_sys = True
    run_cut_variation_sys = False

    # run_bkg_cut_sys = True
    run_bkg_cut_sys = False

    # run_cos_range_sys = True
    run_cos_range_sys = False
    
    if run_track_sys == True:
        os.system(f"python3 Analysis.py -c config_sys_track_0 2>&1 | tee {outstd_dir}/output-track_sys-0.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_1 2>&1 | tee {outstd_dir}/output-track_sys-1.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_2 2>&1 | tee {outstd_dir}/output-track_sys-2.log")
        clear_vars()   

        os.system(f"python3 Analysis.py -c config_sys_track_3 2>&1 | tee {outstd_dir}/output-track_sys-3.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_4 2>&1 | tee {outstd_dir}/output-track_sys-4.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_5 2>&1 | tee {outstd_dir}/output-track_sys-5.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_6 2>&1 | tee {outstd_dir}/output-track_sys-6.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_7 2>&1 | tee {outstd_dir}/output-track_sys-7.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_8 2>&1 | tee {outstd_dir}/output-track_sys-8.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_9 2>&1 | tee {outstd_dir}/output-track_sys-9.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_track_10 2>&1 | tee {outstd_dir}/output-track_sys-10.log")
        clear_vars()

    if run_fit_sys == True:

        os.system(f"python3 Analysis.py -c config_sys_fit_1 2>&1 | tee {outstd_dir}/output-fit_sys-1.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_2 2>&1 | tee {outstd_dir}/output-fit_sys-2.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_3 2>&1 | tee {outstd_dir}/output-fit_sys-3.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_4 2>&1 | tee {outstd_dir}/output-fit_sys-4.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_5 2>&1 | tee {outstd_dir}/output-fit_sys-5.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_6 2>&1 | tee {outstd_dir}/output-fit_sys-6.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_7 2>&1 | tee {outstd_dir}/output-fit_sys-7.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_8 2>&1 | tee {outstd_dir}/output-fit_sys-8.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_9 2>&1 | tee {outstd_dir}/output-fit_sys-9.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_fit_10 2>&1 | tee {outstd_dir}/output-fit_sys-10.log")
        clear_vars()

    if run_cut_variation_sys == True:

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_1 2>&1 | tee {outstd_dir}/output-cut_variation_sys-1.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_2 2>&1 | tee {outstd_dir}/output-cut_variation_sys-2.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_3 2>&1 | tee {outstd_dir}/output-cut_variation_sys-3.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_4 2>&1 | tee {outstd_dir}/output-cut_variation_sys-4.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_5 2>&1 | tee {outstd_dir}/output-cut_variation_sys-5.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_6 2>&1 | tee {outstd_dir}/output-cut_variation_sys-6.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_7 2>&1 | tee {outstd_dir}/output-cut_variation_sys-7.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cut_variation_8 2>&1 | tee {outstd_dir}/output-cut_variation_sys-8.log")
        clear_vars()

    if run_bkg_cut_sys == True:

        os.system(f"python3 Analysis.py -c config_sys_bkg_cut_1 2>&1 | tee {outstd_dir}/output-bkg_cut_sys-1.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_bkg_cut_2 2>&1 | tee {outstd_dir}/output-bkg_cut_sys-2.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_bkg_cut_3 2>&1 | tee {outstd_dir}/output-bkg_cut_sys-3.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_bkg_cut_4 2>&1 | tee {outstd_dir}/output-bkg_cut_sys-4.log")
        clear_vars()

    if run_cos_range_sys == True:

        os.system(f"python3 Analysis.py -c config_sys_cos_range_1 2>&1 | tee {outstd_dir}/output-cos_range_sys-1.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cos_range_2 2>&1 | tee {outstd_dir}/output-cos_range_sys-2.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cos_range_3 2>&1 | tee {outstd_dir}/output-cos_range_sys-3.log")
        clear_vars()

        os.system(f"python3 Analysis.py -c config_sys_cos_range_4 2>&1 | tee {outstd_dir}/output-cos_range_sys-4.log")
        clear_vars()
