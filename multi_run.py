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
# def make_config(file_name,pars_dict):

#     os.system("cp Configurations_PbPb.py Sys_config/PbPb/"+file_name+".py -rf")

#     for key in pars_dict:
#         key_re = re.compile(rf"{chr(34)}{key}{chr(34)}\s*:\s*[^\n]*")

#         config_file = "Sys_config/PbPb/" + file_name + ".py"
#         config = open(config_file,"r+")
#         config_content = config.read()

#         try:
#             key_value = key_re.search(config_content).group()[:-1]
#         except:
#             print(f"Key {key} not found in the configuration file")
#             continue

#         new_content = config_content.replace(key_value,f'"{key}": {pars_dict[key]}')

#         config.seek(0)
#         config.write(new_content)
#         config.truncate()
#         config.close()

def make_config(file_name, pars_dict, con_dir):

    os.system(f"cp Configurations_PbPb.py Sys_config/{con_dir}/{file_name}.py -rf")

    for key in pars_dict:
        # key_re = re.compile(rf"{chr(34)}{key}{chr(34)}\s*:\s*[^\n]*")

        # config_file = "Sys_config/PbPb/" + file_name + ".py"
        config_file = "./Configurations_PbPb.py"
        config = open(config_file,"r+")
        config_content = config.read()

        # for key in pars_dict:

        # pattern = rf'("{key}"\s*:\s*)([^#\n]*?)(?=,)'
        pattern = rf'("{key}"\s*:\s*)(.+)(,)'

        replacement = rf'\g<1>{repr(pars_dict[key])}\g<3>'
        # replacement = rf'\1{repr(pars_dict[key])}'

        # pattern = rf'("{key}"\s*:\s*)(.*?)(,?\s*(#.*)?$)'
        # replacement = rf'\1{repr(pars_dict[key])}\3'

        new_content, n = re.subn(pattern, replacement, config_content)

        if n == 0:
            print(f"Key {key} not found in the configuration file")
        else:
            print(f"Replaced {n} occurrence(s) of {key}")

        # try:
        #     key_value = key_re.search(config_content).group()[:-1]
        # except:
        #     print(f"Key {key} not found in the configuration file")
        #     continue

        # new_content = config_content.replace(key_value,f'"{key}": {pars_dict[key]}')

        config.seek(0)
        config.write(new_content)
        config.truncate()
        config.close()

if __name__ == "__main__":

    rapidity_cut = [
        [0, 0.8],
        [0, 0.3],
        [0.3, 0.8]
    ]
    frac_cuts = [
        "max-0.1",
        "min+0.1",
        "min+0.2",
        "S-gaussian_B-expopowext",
        "S-gaussian_B-expopow",
        "step_0.04"
    ]
    frac_name = [
        "max-0p1",
        "min+0p1",
        "min+0p2",
        "S-gaussian_B-expopowext",
        "S-gaussian_B-expopow",
        "step_0p04"
    ]
    bkg_func = [
        ["expopowext"],
        ["expopow"]
    ]
    signal_func = [
        ["gaussian"]
    ]
    mass_range_list = [
        [0.1396, 0.170]
    ]

    for rap_set in rapidity_cut:
        # for mass_range in mass_range_list:
        # for frac_set in frac_cuts:
            # for signal_func_set in signal_func:
        variation_name = f"rapidity_{int(rap_set[0]*10)}-{int(rap_set[1]*10)}_Data"
        Ana_name_control = {
            "Name_fraction" :   f"PbPb_30-50_0-0.8_FT0C",
            "Name_writing" :    f"PbPb_30-50_{rap_set[0]}-{rap_set[1]}_FT0C_Mix_Cands_LHC23",
            "Name_massfit" :    f"PbPb_30-50_{rap_set[0]}-{rap_set[1]}_FT0C_Mix_Cands_LHC23",
            "Name_rhoextract" : f"PbPb_30-50_{rap_set[0]}-{rap_set[1]}_FT0C_Mix_Cands_LHC23"
            }
        running_control = {
            "Cut_Variations":       False,
            "Rewriting_for_frac":   False,
            "Data_And_Efficiency":  True,
            "Signal_Extraction":    True,
            "Rho_Extraction":       True
            }
        pt_control = {
            "1-2":     False,
            "2-3":     True,
            "3-4":     True,
            "4-6":     True,
            "6-8":     True,
            "8-10":    True,
            "10-15":   True,
            "15-30":   False,
            "30-50":   False
        }
        variation_config = {
            "rapidity_cut" : rap_set
            }

        con_dir = "Mix_Cands"
        print(f"config_{variation_name}")
        # make_config(f"config_{variation_name}", Ana_name_control | running_control | variation_config | pt_control, con_dir)
        os.system(f"python3 Analysis.py -c config_{variation_name} 2>&1 | tee /home/mingze/work/dstar/Dstar_Spin_Alignment/Output/PbPb_Analysis/Conf_And_Logs/Run_log_{variation_name}.log")
        clear_vars()
