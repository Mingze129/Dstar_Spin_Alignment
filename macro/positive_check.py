import ROOT
import os
import numpy as np
from tool import plot_ratio

pt_bins = [3,5,7,10,20,30,50,100]

steps = [0.025,0.04,0.04, 0.05, 0.05, 0.06, 0.08]



for pt_min, pt_max ,step in zip(pt_bins[:-1],pt_bins[1:],steps):   

    fd_list = np.arange(0.0, 1.0001, step)


    out_name = f"positive-check_pt_{pt_min}_{pt_max}.root"

    yield_dir = f"../Output/Dev_test/Cut-variation/pt_{pt_min}_{pt_max}/raw_yield"
    eff_dir = f"../Output/Dev_test/Cut-variation/pt_{pt_min}_{pt_max}/efficiency"

    bin_num = len(fd_list)

    
    hist_check= ROOT.TH1F(f"raw_yield_pt_{pt_min}_{pt_max}", f"raw_yield",bin_num, -0.5, bin_num-0.5)
    eff_p_check=ROOT.TH1F(f"p_eff_pt_{pt_min}_{pt_max}", f"eff_p", bin_num, -0.5, bin_num-0.5)
    eff_np_check=ROOT.TH1F(f"np_eff_pt_{pt_min}_{pt_max}", f"eff_np", bin_num, -0.5, bin_num-0.5)

    for i in range(0,bin_num-1):

        try:
            yield_file =  ROOT.TFile(f"{yield_dir}/{i}_raw_yield_{fd_list[i]:.2f}_{fd_list[i+1]:.2f}.root", "read")
            data_yield = yield_file.Get("hraw_yield")
        except:
            pass

        eff_file = ROOT.TFile(f"{eff_dir}/{i}_eff_{fd_list[i]:.2f}_{fd_list[i+1]:.2f}.root", "read")
        
        
        prompt_eff = eff_file.Get("heff_prompt")
        nonprompt_eff = eff_file.Get("heff_nonprompt")

  
        try:
            hist_check.SetBinContent(i+1, data_yield.GetBinContent(1))
        except:
            hist_check.SetBinContent(i+1, 0)
    
        eff_p_check.SetBinContent(i+1, prompt_eff.GetBinContent(1))
        eff_p_check.SetBinError(i+1, prompt_eff.GetBinError(1))
        eff_np_check.SetBinContent(i+1, nonprompt_eff.GetBinContent(1))
        eff_np_check.SetBinError(i+1, nonprompt_eff.GetBinError(1))

        yield_file.Close()
        eff_file.Close()

    posi_check = ROOT.TFile(out_name, "recreate")
    posi_check.cd()


    hist_check.Write()
    eff_p_check.Write()
    eff_np_check.Write()

    posi_check.Close()
