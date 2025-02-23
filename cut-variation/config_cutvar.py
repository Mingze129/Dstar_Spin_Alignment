import numpy as np
import sys
sys.path.append("..")
import Configurations

gen_set = Configurations
bin_set = Configurations.BinSet["pt_bin_set"]
task_name = Configurations.Analysis["Task_Name"]

pt_min = 30
pt_max = 50
start = 1
end = 12

fd_edges = bin_set[f"{pt_min}-{pt_max}"]["var_fd_range"]
icut = np.arange(0,len(fd_edges),1)

raw_yield_list = [f"{i}_raw_yield_{fd_min:.2f}_{fd_max:.2f}.root" for (i, fd_min,fd_max) in zip(icut[start:end],fd_edges[start:end-1],fd_edges[start+1:end])]
eff_list = [ f"{i}_eff_{fd_min:.2f}_{fd_max:.2f}.root" for (i, fd_min,fd_max) in zip(icut[start:end],fd_edges[start:end-1],fd_edges[start+1:end])]

raw_yield_list.remove(f"6_raw_yield_{fd_edges[6]:.2f}_{fd_edges[7]:.2f}.root"),
eff_list.remove(f"6_eff_{fd_edges[6]:.2f}_{fd_edges[7]:.2f}.root")

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
        "inputfile": f"0_eff_0.00_{fd_edges[1]:.2f}.root",
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