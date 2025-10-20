import os
import numpy as np

working_dir = "/home/mingze/work/dstar/Dstar_Spin_Alignment"

Analysis = {
    "Task_Name": "PbPb_Analysis",
    "Ana_name": "pp_test",
    # "Framework": ["Helicity", "EP"]
    "Framework": ["Helicity"]
}

Doing = {
    "Cut_Variations": False,
    "Data_And_Efficiency": True,
    "Signal_Extraction": True,
    "Rho_Extraction": False,
    "3-5":     False,
    "5-7":     False,
    "7-10":    False,
    "10-20":   False,
    "20-30":   True,
    "30-50":   False,
    "50-100":  False
}

Files = {
    "data": [
            #  "AnalysisResults_LHC22o_apass7.root",
            #  "AnalysisResults_LHC23_apass4_part1.root",
            #  "AnalysisResults_LHC23_apass4_part2.root",
            #  "AnalysisResults_LHC23_apass4_part3.root",
            #  "AnalysisResults_LHC23_apass4_part4.root"
             "AnalysisResults_LHC23_apass4_PbPb_489048.root"
    ],
    "mc": [
            #  "AnalysisResults_mc_merged_pt-30-cut_vsPhi.root"
           "AnalysisResults_mc_LHC25a2_ptsmear-1p5_vsPhi_PbPb.root"     
    ],
    "mc_fraction": [1.0], # if multiple mc files, please provide the fraction of each file
    "simulation": [
        "Pythia_EvtGen_pTHardBins.root"
    ],
}

BinSet = {
    "Mc_reweight": False,
    "Min_cls_ITS": 0, #{1:4, 2:5, 3:6, 4:7}
    "Min_cls_TPC": 0, #{1:80, 2:100, 3:120}
    "Min_eta_track": 0, #{1:0, 2:0.1, 3:0.2}
    "pt_bin_num": 5,
    "pt_bin_edges": [3, 5, 7, 10, 20, 30],
    "pt_bin_set": {
        "3-5": { 
            "doing": Doing["3-5"], 
            "min": 3,
            "max": 5,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1,-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.2,0.6,0.8,1.0],
            "Bkg_cut": 0.06,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [False],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.025),
            "frac_min_bin": 10,
            "frac_max_bin": 38,
            "frac_remove_bin": []
        },
        "5-7": {   
            "doing": Doing["5-7"],      
            "min": 5,
            "max": 7,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.2,0.4,0.8,1.0],
            "Bkg_cut": 0.06,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [False],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.160],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.03),
            "frac_min_bin": 8,
            "frac_max_bin": 30,
            "frac_remove_bin": []
        },
        "7-10": {      
            "doing": Doing["7-10"], 
            "min": 7,
            "max": 10,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.15,0.4,0.8,1.0],
            "Bkg_cut": 0.06,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [False],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.170],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.03),
            "frac_min_bin": 8,
            "frac_max_bin": 30,
            "frac_remove_bin": []
        },
        "10-20": {  
            "doing": Doing["10-20"],       
            "min": 10,
            "max": 20,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,1.0],
            "Bkg_cut": 0.06,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.170],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "Rebin": 1,
            "corr_bkg": [False],
            "var_fd_range": np.arange(0.0, 1.0001, 0.1),
            "frac_min_bin": 1,
            "frac_max_bin": 9,
            "frac_remove_bin": []
        },
        "20-30": {   
            "doing": Doing["20-30"],   
            "min": 20,
            "max": 30,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0, 1.0],
            "Bkg_cut": 0.4,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopow"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.175],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": [0.0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.18,0.22,0.26,0.30,0.34,0.38,0.42,0.46,0.50,0.54,0.58,0.62,0.66,0.70,0.74,0.78,0.82,0.86,0.90,0.94,0.98,1.0],
            "frac_min_bin": 4,
            "frac_max_bin": 22,
            "frac_remove_bin": []
        },
        "30-50": {  
            "doing": Doing["30-50"],     
            "min": 30,
            "max": 50,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,1.0],
            "Bkg_cut": 0.4,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopow"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "with_bkg": False,
            "Mass_range": [0.1396, 0.180],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": [0.0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.18,0.22,0.26,0.30,0.34,0.38,0.42,0.46,0.50,0.54,0.58,0.62,0.66,0.70,0.74,0.78,0.82,0.86,0.90,0.94,0.98,1.0],
            "frac_min_bin": 4,
            "frac_max_bin": 26,
            "frac_remove_bin": []
        },
        "50-100": {    
            "doing": Doing["50-100"],  
            "min": 50,
            "max": 100,
            "D0Mass": [1.74,1.99],
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,1.0],
            "Bkg_cut": 0.15,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopow"],
            "chi2_loss": False,
            "fix_pars": [True,"nl","nr","alphal","alphar"],
            "Mass_range": [0.1396, 0.180],
            "bin_counting": [True, 0.1396, 0.165],
            "threshold": [False],
            "with_bkg": False,
            "Rebin": 3,
            "corr_bkg": [False],
            "var_fd_range": [0.0,0.005,0.01,0.015,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.18,0.22,0.26,0.30,0.34,0.38,0.42,0.46,0.50,0.54,0.58,0.62,0.66,0.70,0.74,0.78,0.82,0.86,0.90],
            "frac_min_bin": 4,
            "frac_max_bin": 25,
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

Do_sys_uncer = {
    "Track_sys": {
        "do": False,
        "Min_cls_ITS": [0, 1, 2, 3, 4], #{1:4, 2:5, 3:6, 4:7}
        "Min_cls_TPC": [0, 1, 2, 3], #{1:80, 2:100, 3:120}
        "Min_eta_track": [0, 1, 2, 3], #{1:0, 2:0.1, 3:0.2}
    },
    "Fit_sys": {
        "do": False
    },
    "Cut_variation_sys": {
        "do": False
    },
    "Bkg_cut_sys": {
        "do": False,
        "bkg_cut": {
            "3-5": [0.02, 0.04, 0.06, 0.08],
            "5-7": [0.02, 0.04, 0.06, 0.08],
            "7-10": [0.02, 0.04, 0.06, 0.08],
            "10-20": [0.02, 0.04, 0.06, 0.08],
            "20-30": [0.02, 0.04, 0.06, 0.08],
            "30-50": [0.02, 0.04, 0.06, 0.08],
            "50-100": [0.02, 0.04, 0.06, 0.08]
        }
    },
    "Cos_range_sys": {
        "do": False
    }
}

Task_name = "hf-task-charm-polarisation"
# Task_name = "task-polarisation-charm-hadrons"
Data_keep_frame = ["Helicity", "EP"]
Data_keep_axis = np.array([0, 4, 5, 6, 7, 8, 9, 10], dtype=np.int32) # mass, mD0, cosstar, mlbkg, mlfd, mintrack, Nits, Ntpc
Force_Reducing = False
Cre_fit_root= True