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

args = argparse.ArgumentParser(description="Run the analysis")
args.add_argument("-c", "--config", help="Configuration file", default="Configurations")
args = args.parse_args()

sys.path.append(os.path.join(os.getcwd(),"Sys_config"))

if args.config:
    config = __import__(args.config)
    

out_dir = os.path.join(os.getcwd(),"Output","Mc_pars")
os.makedirs(out_dir, exist_ok=True)

logger_file = os.path.join(out_dir, f"Run_{datetime.now().strftime('%Y-%m-%d')}.log")
logger = logger_config(log_path = logger_file, log_name = "Running Control")

logger.info("")
logger.info("-"*50)
logger.info("This is where the analyszer started!")
logger.info("-"*50)
logger.info("Making work dir and copy configuration file to it....")

os.system(f"cp {config.__file__} {out_dir} -f")

data_ops = DataOps(config,logger_config(log_path = logger_file, log_name = "Data Operation"))
fit_ops = FitOps(config,logger_config(log_path = logger_file, log_name = "Fitting Operation"))

# logger.info("Writing mc into analysis file...")
# data, mc = data_ops.load_data()
# data_ops.write_mc_pars(mc)

fit_ops.get_mc_pars()

logger.info("All done!")