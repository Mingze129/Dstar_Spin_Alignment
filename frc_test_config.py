import os
import sys
import numpy as np

Doing = {
    "3-5":     True,
    "5-7":     False,
    "7-10":    False,
    "10-20":   False,
    "20-30":   False,
    "30-50":   False,
    "50-100":  False
}

Analysis = {
    # "Framework": ["Helicity","Production","Beam"],
    "Framework": ["Helicity"],
    "Task_Name": "frc_write_test",
    "Mc_reweight": False,
    "Cos_Bin_Bkg": False,
    "MC_only": False,
    "Eff_Correction": False,
    "Cre_fit_root": True
}

Target = {
    "data": ["AnalysisResults_LHC23_apass7.root"],
    "mc": ["AnalysisResults_LHC24g5_apass7_ptsmear-1p5_wosreDCA-1p1.root"],
    "mc_factor": [1],
    "simulation": ["Pythia_EvtGen_pTHardBins.root"],
    "sim_factor": [1]
}

Directories = {
    # "WorkDir": os.getcwd(),
    "WorkDir": "/home/mingze/work/dstar/Dstar_Spin_Alignment",
    "LogDir": os.path.join(os.getcwd(), "Logs"),
    "InputDir": os.path.join("/home/mingze/work/dstar/Dstar_Spin_Alignment", "Input"),
    "OutputDir": os.path.join(os.getcwd(), "Output"),
}

Weights = {
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

BinSet = {
    "pt_bin_num": 7,
    "pt_bin_edges": [3, 5, 7, 10, 20, 30, 50, 100],
    "pt_bin_set": {
        "3-5": { 
            "doing": Doing["3-5"], 
            "min": 3,
            "max": 5,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.2,0.6,0.8,1.0],
            "Bkg_cut": 0.1,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.160],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.5)
        },
        "5-7": {   
            "doing": Doing["5-7"],      
            "min": 5,
            "max": 7,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.2,0.4,0.8,1.0],
            "Bkg_cut": 0.2,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.160],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.04)
        },
        "7-10": {      
            "doing": Doing["7-10"], 
            "min": 7,
            "max": 10,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.15,0.4,0.8,1.0],
            "Bkg_cut": 0.2,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.170],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.04)
        },
        "10-20": {  
            "doing": Doing["10-20"],       
            "min": 10,
            "max": 20,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.15,0.4,0.8,1.0],
            "Bkg_cut": 0.4,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.170],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.05)
        },
        "20-30": {   
            "doing": Doing["20-30"],   
            "min": 20,
            "max": 30,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,0.15,0.4,0.8,1.0],
            "Bkg_cut": 0.4,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.175],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.05)
        },
        "30-50": {  
            "doing": Doing["30-50"],     
            "min": 30,
            "max": 50,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,1.0],
            "Bkg_cut": 0.4,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.180],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.06)
        },
        "50-100": {    
            "doing": Doing["50-100"],  
            "min": 50,
            "max": 100,
            "cos_bin_edges": [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1],
            "fd_edges": [0.0,1.0],
            "Bkg_cut": 0.4,
            "Signal_func": ["doublecb"],
            "Bkg_func": ["expopowext"],
            "fix_sigma": False,
            "rot_bkg_used": False,
            "Mass_range": [0.1396, 0.180],
            "Rebin": 1,
            "var_fd_range": np.arange(0.0, 1.0001, 0.08)
        }
    }
    
}

Data_keep_axis = np.array([0, 4, 5, 6, 7], dtype=np.int32)
Force_Reducing = False

