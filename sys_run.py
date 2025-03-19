import os
import sys
import argparse
import numpy as np
import ctypes
from datetime import datetime
import re

import ROOT
ROOT.gROOT.SetBatch(True)

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

    make_track_config("config_sys_track_0",{"Ana_name": "\"SysAna_track_0_wo_bkg\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 0})
    os.system(f"python3 Analysis.py -c config_sys_track_0 2>&1 | tee output.log")

    # make_track_config("config_sys_track_1",{"Ana_name": "\"SysAna_track_1\"","Min_cls_ITS": 1, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_1 2>&1 | tee output.log")

    # make_track_config("config_sys_track_2",{"Ana_name": "\"SysAna_track_2\"","Min_cls_ITS": 2, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_2 2>&1 | tee output.log")

    # make_track_config("config_sys_track_3",{"Ana_name": "\"SysAna_track_3\"","Min_cls_ITS": 3, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_3 2>&1 | tee output.log")

    # make_track_config("config_sys_track_4",{"Ana_name": "\"SysAna_track_4\"","Min_cls_ITS": 4, "Min_cls_TPC": 0, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_4 2>&1 | tee output.log")

    # make_track_config("config_sys_track_5",{"Ana_name": "\"SysAna_track_5\"","Min_cls_ITS": 0, "Min_cls_TPC": 1, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_5 2>&1 | tee output.log")

    # make_track_config("config_sys_track_6",{"Ana_name": "\"SysAna_track_6\"","Min_cls_ITS": 0, "Min_cls_TPC": 2, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_6 2>&1 | tee output.log")

    # make_track_config("config_sys_track_7",{"Ana_name": "\"SysAna_track_7\"","Min_cls_ITS": 0, "Min_cls_TPC": 3, "Min_eta_track": 0})
    # os.system(f"python3 Analysis.py -c config_sys_track_7 2>&1 | tee output.log")

    # make_track_config("config_sys_track_8",{"Ana_name": "\"SysAna_track_8\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 1})
    # os.system(f"python3 Analysis.py -c config_sys_track_8 2>&1 | tee output.log")

    # make_track_config("config_sys_track_9",{"Ana_name": "\"SysAna_track_9\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 2})
    # os.system(f"python3 Analysis.py -c config_sys_track_9 2>&1 | tee output.log")

    # make_track_config("config_sys_track_10",{"Ana_name": "\"SysAna_track_10\"","Min_cls_ITS": 0, "Min_cls_TPC": 0, "Min_eta_track": 3})
    # os.system(f"python3 Analysis.py -c config_sys_track_10 2>&1 | tee output.log")


