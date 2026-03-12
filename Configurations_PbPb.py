import os
import numpy as np

working_dir = "/home/mingze/work/dstar/Dstar_Spin_Alignment"

Analysis = {
    "Task_Name": "PbPb_Analysis",
    "Name_fraction":    "PbPb_test",
    "Name_writing":     "PbPb_test",
    "Name_massfit":     "PbPb_test",
    "Name_rhoextract":  "PbPb_test",
    # "Framework": ["Helicity", "EP"]
    "Framework": ["EP"]
    # "Framework": ["Helicity", "Production"]
}

Doing = {
    "Cut_Variations":       True,
    "Rewriting_for_frac":   True,
    "Data_And_Efficiency":  True,
    "Signal_Extraction":    False,
    "Rho_Extraction":       False,
    "1-2":     False,
    "2-3":     True,
    "3-5":     True,
    "5-7":     True,
    "7-10":    True,
    "10-20":   True,
    "20-30":   True,
    "30-50":   False
}

Files = {
    "data": [
             "AnalysisResults.root"
    ],
    "mc": [
           "AnalysisResults_MC_PbPb_LHC25a5_547031.root"
    ],
    "mc_fraction": [1.0], # if multiple mc files, please provide the fraction of each file
    "simulation": [
        "AnalysisResults_MC_PbPb_LHC25a5_547031.root"
    ],
}

BinSet = {
    "rapidity_cut": [0.0, 0.8],
    "Mc_reweight": False,
    "Min_cls_ITS": 0, #{1:4, 2:5, 3:6, 4:7}
    "Min_cls_TPC": 0, #{1:80, 2:100, 3:120}
    "Min_eta_track": 0, #{1:0, 2:0.1, 3:0.2}
    "pt_bin_num": 11,
    "pt_bin_edges": [1, 2, 3, 5, 7, 10, 20, 30, 50],
    "using_np_mc" :[True, True, True, True, True, True, True, True, True, True, True],
    # "using_np_mc" :[False, False, False, False, False, False, False, False, False, False, False], # --- IGNORE ---
    "corr_with_v2" :[True, True, True, True, True, True, True, True, True, True, True],
    "pars_fix_from_mc": [True, True, True, True, True, True, True, True, True, True, True],
    "pt_bin_set": {
        "1-2": {
            "doing": Doing["1-2"],
            "min": 1,
            "max": 2,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.004,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 30,
            "frac_remove_bin": []
        },
        "2-3": {
            "doing": Doing["2-3"],
            "min": 2,
            "max": 3,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.004,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 40,
            "frac_remove_bin": []
        },
        "3-5": {
            "doing": Doing["3-5"],
            "min": 3,
            "max": 5,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.01,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 40,
            "frac_remove_bin": []
        },
        "5-7": {
            "doing": Doing["5-7"],
            "min": 5,
            "max": 7,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.03,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 38,
            "frac_remove_bin": []
        },
        "7-10": {
            "doing": Doing["7-10"],
            "min": 7,
            "max": 10,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.15,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 38,
            "frac_remove_bin": []
        },
        "10-20": {
            "doing": Doing["10-20"],
            "min": 10,
            "max": 20,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.8,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 38,
            "frac_remove_bin": []
        },
        "20-30": {
            "doing": Doing["20-30"],
            "min": 20,
            "max": 30,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 1,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.165],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 30,
            "frac_remove_bin": []
        },
        "30-50": {
            "doing": Doing["30-50"],
            "min": 30,
            "max": 50,
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 1,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 38,
            "frac_remove_bin": []
        }
    }
}

Directories = {
    # "WorkDir": os.getcwd(),
    "WorkDir": working_dir,
    "InputDir": os.path.join(working_dir, "Input"),
    "OutputDir": os.path.join(working_dir, "Output"),
}

Weights = {
    "Do_reweight": False,
    "Hard_pt_weights": "ptweights_dstar_LHC24d5_apass4.root",
    "Normal_pt_weights": "ptweights_dstar_LHC24d3_apass6.root",
    "Mult_weights": "MultWeights_LHC24d5_trackTuner_208424.root",
    "B_Reco_axis": 12,
    "C_Reco_axis": 1,
    "Mult_Reco_axis": 2,
    "B_Gen_axis": 4,
    "C_Gen_axis": 0,
    "Mult_Gen_axis": 1
}

Template_fit = False
Plot_np_rho = False
Task_name = "hf-task-charm-polarisation"
Data_keep_frame = ["EP"]
Data_keep_axis = np.array([0, 4, 5, 6, 7, 8, 9, 10, 3], dtype=np.int32) # mass, mD0, cosstar, mlbkg, mlfd, mintrack, Nits, Ntpc, rapidity
Qc_check_axes = {
    "Reduced" : [3, 4, 8],
    "Reco": [1, 3, 4, 6],
    "Gen": [0, 2]
}
Force_Reducing = False
Cre_fit_root= True
EP_Resolution = {
    "doing": True,
    "Redoing": False,
    "Det_A": "FT0M",
    "Det_B": "TPCpos",
    "Det_C": "TPCneg"
}