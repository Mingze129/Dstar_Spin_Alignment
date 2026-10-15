import os
import numpy as np

working_dir = "/home/mingze/work/dstar/Dstar_Spin_Alignment"

Analysis = {
    "Task_Name":        "PbPb_Analysis",
    "Name_fraction":    'PbPb_30-50_0-0.8_FT0C',
    "Name_writing":     'PbPb_30-50_0-0.8_FT0C_Data',
    "Name_massfit":     'PbPb_30-50_0-0.8_FT0C_Data',
    "Name_rhoextract":  'PbPb_30-50_0-0.8_FT0C_Data',
    "Name_Mix_Cands":   'PbPb_30-50_0-0.8_FT0C_Mix_Cands_LHC23',
    # "Framework": ["Helicity", "EP"]
    "Framework": ["EP"]
    # "Framework": ["Helicity", "Production"]
}

Doing = {
    "Cut_Variations":       False,
    "Rewriting_for_frac":   False,
    "Data_And_Efficiency":  False,
    "Signal_Extraction":    False,
    "Rho_Extraction":       True,
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

Files = {
    "Primary_data": [
        # "AnalysisResults_LHC23_apass4_PbPb_20-50_634121_part_0.root",
        # "AnalysisResults_LHC23_apass4_PbPb_20-50_634121_part_1.root",
        # "AnalysisResults_LHC23_apass4_PbPb_20-50_634121_part_2.root",
        # "AnalysisResults_LHC23_apass4_PbPb_20-50_634121_part_3.root",
        # "AnalysisResults_LHC23_apass4_PbPb_20-50_634121_part_4.root",
        # "AnalysisResults_LHC23_apass4_PbPb_20-50_634121_part_5.root",
        # "AnalysisResults_LHC24_apass3_PbPb_20-50_634120_part_0.root",
        # "AnalysisResults_LHC24_apass3_PbPb_20-50_634120_part_1.root",
        # "AnalysisResults_LHC24_apass3_PbPb_20-50_634120_part_2.root",
        # "AnalysisResults_LHC24_apass3_PbPb_20-50_634120_part_3.root",
        # "AnalysisResults_LHC24_apass3_PbPb_20-50_634120_part_4.root"
        "AnalysisResults_LHC24_apass3_PbPb_20-50_669042_FT0C.root",
        "AnalysisResults_LHC23_apass4_PbPb_20-50_669043_FT0C.root"
    ],
    "Reduced_data": [
        # "AnalysisResults_LHC23_apass4_PbPb_30-50_634121_reduced.root",
        # "AnalysisResults_LHC24_apass3_PbPb_30-50_634120_reduced.root"
        # "AnalysisResults_LHC23_apass4_PbPb_30-50_669043_FT0C_reduced.root",
        # "AnalysisResults_LHC24_apass3_PbPb_30-50_669042_FT0C_reduced.root"
        # "Mix_Cands_reduced_20-25.root",
        # "Mix_Cands_reduced_26-29.root",
        "Mix_Cands_reduced_30-50.root"
    ],
    "mc": [
        #    "AnalysisResults_mc_LHC25a5_ptsmear-1p5_vsPhi_20-50_634283.root"
           "AnalysisResults_MC_LHC25a5_1p5_vsPhi_FT0C_669044.root"
    ],
    "mc_fraction": [1.0], # if multiple mc files, please provide the fraction of each file
    "simulation": [
        "npDstar_PbPb_rho_v2_0p1_thsparse.root"
    ],
}

BinSet = {
    "centrality_cut": [20, 50],
    "rapidity_cut": [0, 0.8],
    "Mc_reweight": False,
    "Min_cls_ITS": 1, 
    #{1:4, 2:5, 3:6, 4:7}
    "Min_cls_TPC": 1, 
    #{1:80, 2:100, 3:120}
    "Min_eta_track": 1,
    #{1:0, 2:0.1, 3:0.2}
    # "pt_bin_edges": [1, 2, 3, 5, 7, 10, 20, 30, 50],
    "pt_bin_edges": [1, 2, 3, 4, 6, 8, 10, 15, 30, 50],
    "using_np_mc" :[True, True, True, True, True, True, True, True, True, True, True],
    # "using_np_mc" :[False, False, False, False, False, False, False, False, False, False, False], # --- IGNORE ---
    "corr_with_v2" :[True, True, True, True, True, True, True, True, True, True, True],
    "pars_fix_from_mc": [True, True, True, True, True, True, True, True, True, True, True],
    "pt_bin_set": {
        "1-2": {
            "doing": Doing["1-2"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.02,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 30,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "2-3": {
            "doing": Doing["2-3"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.02,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 48,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "3-4": {
            "doing": Doing["3-4"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.04,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 49,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "4-6": {
            "doing": Doing["4-6"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.2,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 49,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "6-8": {
            "doing": Doing["6-8"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.5,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 15,
            "frac_max_bin": 45,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "8-10": {
            "doing": Doing["8-10"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.5,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 45,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "10-15": {
            "doing": Doing["10-15"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.5,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 45,
            "frac_remove_bin": [],
            "Mix_corr": True
        },
        "15-30": {
            "doing": Doing["15-30"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.5,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 42,
            "frac_remove_bin": [],
            "Mix_corr": False
        },
        "30-50": {
            "doing": Doing["30-50"],
            "D0Mass": [1.74, 1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.9],
            "Bkg_cut": 0.5,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [False, 0.1396, 0.160],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.02),
            "frac_min_bin": 10,
            "frac_max_bin": 38,
            "frac_remove_bin": [],
            "Mix_corr": False
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

# Data_keep_axis = np.array([0, 4, 5, 6, 7, 8, 9, 10, 3], dtype=np.int32) # mass, mD0, cosstar, mlbkg, mlfd, mintrack, Nits, Ntpc, rapidity
# Data_keep_axis = np.array([0, 5, 7, 8, 9, 10, 3], dtype=np.int32) # mass, cosstar, mlfd, mintrack, Nits, Ntpc, rapidity
Reducing_data = {
    "Force_Reducing": False,
    "Data_keep_frame": ["EP"],
    "Data_keep_axis": np.array([0, 5, 7, 8, 9, 10, 3], dtype=np.int32),
    "Mass": 0,
    "Cost": 1,
    "Np_scores": 2,
    "Min_track_eta": 3,
    "N_ITS": 4,
    "N_TPC": 5,
    "Rapidity": 6
}
Qc_check_axes = {
    "Reduced" : [2, 6],
    "Reco": [1, 3, 4, 6],
    "Gen": [0, 2]
}
EP_Resolution = {
    "doing": True,
    "Redoing": False,
    "Det_A": "FT0C",
    "Det_B": "FT0A",
    "Det_C": "TPCtot"
}