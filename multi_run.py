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
def make_config(file_name,pars_dict):

    os.system("cp Configurations_PbPb.py Sys_config/PbPb/"+file_name+".py -rf")

    for key in pars_dict:
        key_re = re.compile(rf"{chr(34)}{key}{chr(34)}\s*:\s*[^\n]*")

        config_file = "Sys_config/PbPb/" + file_name + ".py"
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

    rapidity_cut = [
        [0.0,0.2],
        [0.2,0.4],
        [0.4,0.8]
    ]
    for set in rapidity_cut:

        variation_name = f"rapidity_{int(set[0]*10)}_{int(set[1]*10)}"
        Ana_name_control = {
            "Name_fraction" :   f"\"PbPb_rapidity_0.0_0.8\"",
            "Name_writing" :    f"\"PbPb_rapidity_{set[0]}_{set[1]}\"",
            "Name_massfit" :    f"\"PbPb_rapidity_{set[0]}_{set[1]}\"",
            "Name_rhoextract" : f"\"PbPb_rapidity_{set[0]}_{set[1]}\""
            }
        running_control = {
            "Cut_Variations":       False,
            "Rewriting_for_frac":   False,
            "Data_And_Efficiency":  False,
            "Signal_Extraction":    False,
            "Rho_Extraction":       True
            }
        variation_config = {
            "rapidity_cut" : set
            }

        make_config(f"config_{variation_name}", Ana_name_control | running_control | variation_config)
        os.system(f"python3 Analysis.py -c config_{variation_name} 2>&1 | tee /home/mingze/work/dstar/Dstar_Spin_Alignment/Output/PbPb_Analysis/Conf_And_Logs/Run_log_{variation_name}.log")
        clear_vars()











        
# Do_sys_uncer = {
#     "Track_sys": {
#         "do": False,
#         "Min_cls_ITS": [0, 1, 2, 3, 4], #{1:4, 2:5, 3:6, 4:7}
#         "Min_cls_TPC": [0, 1, 2, 3], #{1:80, 2:100, 3:120}
#         "Min_eta_track": [0, 1, 2, 3], #{1:0, 2:0.1, 3:0.2}
#     },
#     "Fit_sys": {
#         "do": False
#     },
#     "Cut_variation_sys": {
#         "do": False
#     },
#     "Bkg_cut_sys": {
#         "do": False,
#         "bkg_cut": {
#             "3-5": [0.02, 0.04, 0.06, 0.08],
#             "5-7": [0.02, 0.04, 0.06, 0.08],
#             "7-10": [0.02, 0.04, 0.06, 0.08],
#             "10-20": [0.02, 0.04, 0.06, 0.08],
#             "20-30": [0.02, 0.04, 0.06, 0.08],
#             "30-50": [0.02, 0.04, 0.06, 0.08],
#             "50-100": [0.02, 0.04, 0.06, 0.08]
#         }
#     },
#     "Cos_range_sys": {
#         "do": False
#     }
# }