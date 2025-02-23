import os
import sys
import argparse
import numpy as np
import ctypes

import ROOT
ROOT.gROOT.SetBatch(True)

from utils import DataOps, FitOps, SpinOps, FracOps, logger_config
import Configurations as config

outdir = out_dir = os.path.join(config.Directories["OutputDir"],config.Analysis["Task_Name"])
os.makedirs(out_dir, exist_ok=True)
logger = logger_config(log_path = os.path.join(outdir, "Run.log"), log_name = "Running Control")

logger.info("")
logger.info("-"*50)
logger.info("This is where the analyszer started!")
logger.info("-"*50)
logger.info("Making work dir and copy configuration file to it....")

os.system(f"cp Configurations.py {outdir} -f")

data_ops = DataOps(config,logger_config(log_path = os.path.join(outdir, "Run.log"), log_name = "Data Operation"))
fit_ops = FitOps(config,logger_config(log_path = os.path.join(outdir, "Run.log"), log_name = "Fitting Operation"))
spin_ops = SpinOps(config,logger_config(log_path = os.path.join(outdir, "Run.log"), log_name = "Analysing Operation"))
frac_ops = FracOps(config,logger_config(log_path = os.path.join(outdir, "Run.log"), log_name = "Cut-Variation Operation"))



logger.info("Get efficiency and raw-yield for cut-variation")
# frac_ops.get_input()

# frac_ops.get_fraction()

logger.info("Writing data into analysis file...")
data, mc = data_ops.load_data()
# data_ops.write_data(data)
# data_ops.write_mc(mc)


# fit_ops.get_raw_yield()

data_ops.read_frac()

spin_ops.get_rho()
# spin_ops.write_simu_rho()



