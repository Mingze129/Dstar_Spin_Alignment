import os
import sys
import argparse
import numpy as np
import ctypes
from datetime import datetime

now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

import ROOT
ROOT.gROOT.SetBatch(True)

from utils import DataOps, FitOps, SpinOps, FracOps, logger_config
import Configurations as config

args = argparse.ArgumentParser(description="Run the analysis")
args.add_argument("-c", "--config", help="Configuration file", default="Configurations")
args = args.parse_args()

sys.path.append(os.path.join(os.getcwd(),"Sys_config"))

if args.config:
    config = __import__(args.config)
    
out_dir = os.path.join(config.Directories["OutputDir"],config.Analysis["Task_Name"])
os.makedirs(out_dir, exist_ok=True)

logger_file = os.path.join(out_dir, f"Run_{datetime.now().strftime('%Y-%m-%d')}.log")
logger = logger_config(log_path = logger_file, log_name = "Running Control")

logger.info("")
logger.info("-"*50)
logger.info("This is where the analyszer started!")
logger.info(f"Task Name: {config.Analysis['Task_Name']}, Analysis: {config.Analysis['Ana_name']}")
logger.info("-"*50)
logger.info("Creating work dir and copy configuration file to it....")

os.system(f"cp {config.__file__} {out_dir}/Configurations_{datetime.now().strftime('%Y-%m-%d')}.py -f")

if config.Doing["Cut_Variations"]:
    logger.info("Initiating cut-variation operator...")
    frac_ops = FracOps(config,logger_config(log_path = logger_file, log_name = "Cut-Variation Operation"))
    logger.info("Get input for cut-variation...")
    frac_ops.get_yield_input()
    frac_ops.get_eff_input()
    logger.info("Get fraction by cut-variation method...")
    frac_ops.get_fraction()

if config.Doing["Data_And_Efficiency"]:
    logger.info("Initiating data and efficiency operator...")
    data_ops = DataOps(config,logger_config(log_path = logger_file, log_name = "Data Operation"))
    logger.info("Writing data into analysis file...")
    data, mc = data_ops.load_data()
    data_ops.write_data(data)
    data_ops.write_mc(mc)
    data_ops.write_corr_fit() # Writing the correlated background template

if config.Doing["Signal_Extraction"]:
    logger.info("Initiating fit operator...")
    fit_ops = FitOps(config,logger_config(log_path = logger_file, log_name = "Fitting Operation"))
    logger.info("Fitting raw-yield...")
    fit_ops.get_raw_yield()

if config.Doing["Rho_Extraction"]:
    logger.info("Initiating rho extraction operator...")
    spin_ops = SpinOps(config,logger_config(log_path = logger_file, log_name = "Analysing Operation"))
    logger.info("Get corrected yield...")
    spin_ops.read_frac()
    logger.info("Extracted rho by fit the corrected yield...")
    spin_ops.get_rho()
    spin_ops.extro_rho()
    spin_ops.plot_rho()