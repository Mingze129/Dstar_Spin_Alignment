import numpy as np
import sys
sys.path.append("..")
import Configurations

gen_set = Configurations
bin_set = Configurations.BinSet["pt_bin_set"]
task_name = Configurations.Analysis["Task_Name"]

pt_min = 7
pt_max = 10
start_bin = 3
end_bin = 24
remove_bin = []

fd_edges = bin_set[f"{pt_min}-{pt_max}"]["var_fd_range"]
icut = np.arange(0,len(fd_edges),1)

raw_yield_list = [f"{i}_raw_yield_fd-cut_{fd_min:.3f}.root" for (i, fd_min,fd_max) in zip(icut[start_bin:end_bin],fd_edges[start_bin:end_bin-1],fd_edges[start_bin+1:end_bin])]
eff_list = [ f"{i}_efficiency_fd-cut_{fd_min:.3f}.root" for (i, fd_min,fd_max) in zip(icut[start_bin:end_bin],fd_edges[start_bin:end_bin-1],fd_edges[start_bin+1:end_bin])]

for nbin in remove_bin:
  raw_yield_list.remove(f"{nbin}_raw_yield_fd-cut_{fd_edges[nbin]:.3f}.root"),
  eff_list.remove(f"{nbin}_efficiency_fd-cut_{fd_edges[nbin]:.3f}.root")

cutvar_set = {
    "rawyields": {
      "inputdir": f"../Output/{task_name}/Cut-variation/pt_{pt_min}_{pt_max}/raw_yield",
      "inputfiles": raw_yield_list,
      "histoname": "hraw_yield"
    },
    "efficiencies": {
      "inputdir": f"../Output/{task_name}/Cut-variation/pt_{pt_min}_{pt_max}/efficiency",
      "inputfiles": eff_list,
      "histonames": {
        "prompt": "heff_prompt",
        "nonprompt": "heff_nonprompt"
      }
    },
    "minimisation": {
        "correlated": True
    },
    "central_efficiency": {
        "computerawfrac": True,
        "inputdir": f"../Output/{task_name}/Cut-variation/pt_{pt_min}_{pt_max}/efficiency",
        "inputfile": f"0_efficiency_fd-cut_{fd_edges[0]:.3f}.root",
        "histonames": {
            "prompt": "heff_prompt",
            "nonprompt": "heff_nonprompt"
        }
    },
    "output": {
        "directory": f"../Output/{task_name}/Cut-variation/pt_{pt_min}_{pt_max}/fraction",
        "file": f"CutVar_{task_name}.root"
    }
}