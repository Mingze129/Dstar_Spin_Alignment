import os
import sys
import numpy as np
from tqdm import tqdm
import ctypes
import ROOT
ROOT.gROOT.SetBatch(True)


class McOps(object):

    def __init__(self,config):
        self.config = config
        self.weights_dir = os.path.join(self.config.Directories["InputDir"], "Weight")
        self.fraction = config.Target["mc_factor"]

    def get_mc_sparse(self, mc_list, frame):
        
        Reco_prompt_mc = list()
        Reco_nonprompt_mc = list()
        Gen_prompt_mc  = list()
        Gen_nonprompt_mc  = list()
        for (mc_file, frac) in zip(mc_list,self.fraction):

            (f"\nLoading MC samples from {mc_file.split('/')[-1]} ...")
            infile = ROOT.TFile(mc_file, "READ")
            if not infile:
                print(f"MC file {mc_file} not found")
                sys.exit(1)

            Reco_prompt = infile.Get(f"task-polarisation-charm-hadrons/hRecoPrompt{frame}")
            Reco_nonprompt = infile.Get(f"task-polarisation-charm-hadrons/hRecoNonPrompt{frame}")
            Gen_prompt = infile.Get(f"task-polarisation-charm-hadrons/hGenPrompt{frame}")
            Gen_nonprompt = infile.Get(f"task-polarisation-charm-hadrons/hGenNonPrompt{frame}")
            if self.config.Analysis["Mc_reweight"] == False:
                Reco_prompt_mc.append(Reco_prompt)
                Reco_nonprompt_mc.append(Reco_nonprompt)
                Gen_prompt_mc.append(Gen_prompt)
                Gen_nonprompt_mc.append(Gen_nonprompt)
            else:
                Reco_prompt = self.do_reweight(Reco_prompt, frac, mc_file,self.config.Weights["C_Reco_axis"], self.config.Weights["Mult_Reco_axis"])
                Reco_nonprompt = self.do_reweight(Reco_nonprompt, frac, mc_file,self.config.Weights["B_Reco_axis"], self.config.Weights["Mult_Reco_axis"])
                Gen_prompt = self.do_reweight(Gen_prompt, frac, mc_file,self.config.Weights["C_Gen_axis"], self.config.Weights["Mult_Gen_axis"])
                Gen_nonprompt =self.do_reweight(Gen_nonprompt, frac, mc_file,self.config.Weights["B_Gen_axis"], self.config.Weights["Mult_Gen_axis"])

                Reco_prompt_mc.append(Reco_prompt)
                Reco_nonprompt_mc.append(Reco_nonprompt)
                Gen_prompt_mc.append(Gen_prompt)
                Gen_nonprompt_mc.append(Gen_nonprompt)

        
        Reco_prompt_all, Reco_nonprompt_all, Gen_prompt_all, Gen_nonprompt_all = self.merge(Reco_prompt_mc, Reco_nonprompt_mc, Gen_prompt_mc, Gen_nonprompt_mc)
        return Reco_prompt_all, Reco_nonprompt_all, Gen_prompt_all, Gen_nonprompt_all
      
    def do_reweight(self, sparse, frac, file_dir, pt_axis, multi_axis):

        print(f"Reweighting {sparse.GetName()}...")

        if "hard" in file_dir:
            pt_weight_file = ROOT.TFile(self.weights_dir + '/' + self.config.Weights["Hard_pt_weights"], 'read')
            if 'NonPrompt' in sparse.GetName():
                pt_weight_hist = pt_weight_file.Get("h_weights_bmesons")
            else:
                pt_weight_hist = pt_weight_file.Get("h_weights_prompt")
           
        else:
            pt_weight_file = ROOT.TFile(self.weights_dir + '/' + self.config.Weights["Normal_pt_weights"], 'read')
            if 'NonPrompt' in sparse.GetName():
                pt_weight_hist = pt_weight_file.Get("h_weights_bmesons")
            else:
                pt_weight_hist = pt_weight_file.Get("h_weights_prompt")

        multi_file = ROOT.TFile(self.weights_dir + '/' + self.config.Weights["Mult_weights"], 'read')
        if not multi_file or not pt_weight_file:
            print(f"Weight file not found")
            sys.exit(1)

        multi_weight_hist = multi_file.Get("h_weights_candinmass")

        coord =  [0] * sparse.GetNdimensions()
        coord = np.array(coord, dtype=np.int32)

        for i in tqdm(range(sparse.GetNbins())):

            content = sparse.GetBinContent(i, coord)
            error = sparse.GetBinError(i)
            pt_weight = pt_weight_hist.GetBinContent(int(coord[pt_axis]))
            pt_weight_error = pt_weight_hist.GetBinError(int(coord[pt_axis]))

            multi_weight = multi_weight_hist.GetBinContent(int(coord[multi_axis]))
            multi_weight_error = multi_weight_hist.GetBinError(int(coord[multi_axis]))

            sparse.SetBinContent(i, frac*content*pt_weight*multi_weight)
            sparse.SetBinError(i, frac*error*pt_weight*multi_weight)

        return sparse

    def merge(self, Reco_prompt, Reco_nonprompt, Gen_prompt, Gen_nonprompt):
        
        print("\nMerging MC samples ...")
        for i in range(len(Reco_prompt)):
            if i == 0:
                continue
            Reco_prompt[0].Add(Reco_prompt[i])
            Reco_nonprompt[0].Add(Reco_nonprompt[i])
            Gen_prompt[0].Add(Gen_prompt[i])
            Gen_nonprompt[0].Add(Gen_nonprompt[i])

        Reco_coord = np.arange(0,self.config.Weights["B_Reco_axis"],1)
        Reco_coord = np.array(Reco_coord, dtype=np.int32)
        Reco_THP = Reco_prompt[0].ProjectionND(self.config.Weights["B_Reco_axis"], Reco_coord, "EO")
        Reco_THNP = Reco_nonprompt[0].ProjectionND(self.config.Weights["B_Reco_axis"], Reco_coord, "EO")

        if Reco_THP.GetNdimensions() != Reco_THNP.GetNdimensions():
            print("Dimensions mismatch when merging promtp and non-prompt")
            sys.exit(1)
        
        Gen_coord = np.arange(0,self.config.Weights["B_Gen_axis"],1)
        Gen_coord = np.array(Gen_coord, dtype=np.int32)
        Gen_THP = Gen_prompt[0].ProjectionND(self.config.Weights["B_Gen_axis"], Gen_coord, "EO")
        Gen_THNP = Gen_nonprompt[0].ProjectionND(self.config.Weights["B_Gen_axis"], Gen_coord, "EO")

        if Gen_THP.GetNdimensions() != Gen_THNP.GetNdimensions():
            print("Dimensions mismatch when merging promtp and non-prompt")
            sys.exit(1)

        return Reco_THP, Reco_THNP, Gen_THP, Gen_THNP
    
if __name__ == '__main__':
    
    mc_list = ["/home/mingze/work/dstar/Dstar_Spin_Alignment/Input/MC/AnalysisResults_LHC24g5_apass7_ptsmear-1p5_wosreDCA-1p1.root"]

    sys.path.append("../")
    import Configurations as config
    
    mc_ops = McOps(config)
    mc_ops.get_mc_sparse(mc_list, "Beam")

