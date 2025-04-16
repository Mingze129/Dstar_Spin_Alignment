import os
import sys
import numpy as np
from tqdm import tqdm
import ctypes
import ROOT
ROOT.gROOT.SetBatch(True)

from .Formula import constant_func
from .McOps import McOps

class DataOps(object):

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.out_dir = os.path.join(self.config.Directories["OutputDir"],self.config.Analysis["Task_Name"])
        self.logger.info("Locating work dir...")
        os.makedirs(self.out_dir + "/Analysis-root", exist_ok=True)
        os.makedirs(self.out_dir + "/Cut-variation", exist_ok=True)
        os.makedirs(self.out_dir + "/Figure", exist_ok=True)
        os.makedirs(self.out_dir + "/Mass-Fit", exist_ok=True)
        self.mc_ops = McOps(config)
    

    def load_data(self):
        self.logger.info("Loading data...")
        D_file_list = self.config.Target["data"]
        M_file_list = self.config.Target["mc"]
        self.logger.info(f"target data files: {D_file_list}")
        self.logger.info(f"target mc files: {M_file_list}")
        self.logger.info(f"Finding data in directory:  {self.config.Directories['InputDir']}")

        if not os.path.exists(self.config.Directories["InputDir"]):
            self.logger.error("Input directory does not exist.")
            sys.exit(1)
        try:
            D_file_finded = os.listdir(self.config.Directories["InputDir"]+"/Data/Primary")
            D_file_reduced = os.listdir(self.config.Directories["InputDir"]+"/Data/Reduced")
            M_file_finded = os.listdir(self.config.Directories["InputDir"]+"/MC")
        except:
            self.logger.error("Error when Locate files, please check the directory.")
            sys.exit(1)

        MC_list = []
        for M_file in M_file_list:
            if M_file not in M_file_finded:
                self.logger.error(f"MC file {M_file} not found.")
                sys.exit(1)
            MC_list.append(os.path.join(self.config.Directories["InputDir"],"MC",M_file))

        reduce_key = []
        for D_file in D_file_list:
            reduced_file = D_file.replace(".root", "_reduced.root")
            if reduced_file not in D_file_reduced:
                self.logger.warning(f"Reduced data file {reduced_file} not found!\nTry to reduce data file {D_file}...")
                if D_file not in D_file_finded:
                    self.logger.error(f"Data file {D_file} not found.")
                    sys.exit(1)
                reduce_key.append(True)
            elif self.config.Force_Reducing:
                reduce_key.append(True)
            else:
                reduce_key.append(False)
  
        Reduced_file = []
        for i, (D_file, key) in enumerate(zip(D_file_list, reduce_key)):
            if key:
                self.logger.info(f"Reducing data file {D_file}...")
                Reduced_file.append(self.Reduceing(D_file))
            else:
                Reduced_file.append(os.path.join(self.config.Directories["InputDir"],"Data/Reduced",D_file.replace(".root", "_reduced.root")))

        return Reduced_file, MC_list

    def write_data(self, data):

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        outfile_name = os.path.join(self.out_dir, "Analysis-root", self.config.Analysis["Ana_name"]+".root")

        outfile = ROOT.TFile(outfile_name, "UPDATE")
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            self.logger.info(f"Writing data of {frame} framework...")
            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue
                self.logger.info(f"    Working in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                
                Tdata = self.get_sparse(data, f"h{frame}_pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}_data")
                Tbkg = self.get_sparse(data,f"h{frame}_pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}_rotbkg")
                
                D0_min = Tdata.GetAxis(1).FindBin(pt_bin_set["D0Mass"][0]+1e-6)
                D0_max = Tdata.GetAxis(1).FindBin(pt_bin_set["D0Mass"][1]-1e-6)
                Tdata.GetAxis(1).SetRange(D0_min,D0_max)
                Tbkg.GetAxis(1).SetRange(D0_min,D0_max)

                Tdata.GetAxis(5).SetRange(self.config.Analysis["Min_eta_track"],100)
                Tbkg.GetAxis(5).SetRange(self.config.Analysis["Min_eta_track"],100)
                Tdata.GetAxis(6).SetRange(self.config.Analysis["Min_cls_ITS"],100)
                Tbkg.GetAxis(6).SetRange(self.config.Analysis["Min_cls_ITS"],100)
                Tdata.GetAxis(7).SetRange(self.config.Analysis["Min_cls_TPC"],100)
                Tbkg.GetAxis(7).SetRange(self.config.Analysis["Min_cls_TPC"],100)

                pt_bin_dir.cd()

                bkg_max = Tdata.GetAxis(3).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
                Tdata.GetAxis(3).SetRange(0, bkg_max)
                Tbkg.GetAxis(3).SetRange(0, bkg_max)
                
                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]

                fd_bin = np.array(fd_edges)
              
                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()
        
                    cos_max = Tdata.GetAxis(2).FindBin(1-1e-6)
                    Tdata.GetAxis(2).SetRange(0,cos_max)
                    Tbkg.GetAxis(2).SetRange(0,cos_max)

                    fd_min = Tdata.GetAxis(4).FindBin(fd_min_edge+1e-6)
                    fd_max = Tdata.GetAxis(4).FindBin(fd_max_edge-1e-6)
                    Tdata.GetAxis(4).SetRange(fd_min,fd_max)
                    Tbkg.GetAxis(4).SetRange(fd_min,fd_max)

                    hmass = Tdata.Projection(0).Clone(f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}")
                    hmass.GetXaxis().SetTitle(f"M(K#pi#pi) - M(K#pi) GeV/#it{{c}}^{{2}}")
                    hmass.GetYaxis().SetTitle("Entries")
                    hmass.Write("",ROOT.TObject.kOverwrite)

                    hbkg = Tbkg.Projection(0).Clone(f"hrotbkg_{frame}_pt_{pt_min_edge}_{pt_max_edge}")
                    hbkg.GetXaxis().SetTitle(f"M(K#pi#pi) - M(K#pi) GeV/#it{{c}}^{{2}}")
                    hbkg.GetYaxis().SetTitle("Entries")
                    hbkg.Write("",ROOT.TObject.kOverwrite)

                    plot = self.Compare(hbkg,hmass)
                    plot.Write("",ROOT.TObject.kOverwrite)
              
                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        cos_dir = fd_dir.mkdir(f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", "", ROOT.kTRUE)
                        cos_dir.cd()

                        cos_min = Tdata.GetAxis(2).FindBin(cos_min_edge+1e-6)
                        cos_max = Tdata.GetAxis(2).FindBin(cos_max_edge-1e-6)
                        Tdata.GetAxis(2).SetRange(cos_min,cos_max)
                        Tbkg.GetAxis(2).SetRange(cos_min,cos_max)

                        hmass = Tdata.Projection(0).Clone(f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_cos_{cos_min_edge}_{cos_max_edge}")
                        hmass.GetXaxis().SetTitle(f"M(K#pi#pi) - M(K#pi) GeV/#it{{c}}^{{2}}")
                        hmass.GetYaxis().SetTitle("Entries")
                        hmass.Write("",ROOT.TObject.kOverwrite)

                        hbkg = Tbkg.Projection(0).Clone(f"hrotbkg_{frame}_pt_{pt_min_edge}_{pt_max_edge}_cos_{cos_min_edge}_{cos_max_edge}")
                        hbkg.GetXaxis().SetTitle(f"M(K#pi#pi) - M(K#pi) GeV/#it{{c}}^{{2}}")
                        hbkg.GetYaxis().SetTitle("Entries")
                        hbkg.Write("",ROOT.TObject.kOverwrite)

                        plot = self.Compare(hbkg,hmass)
                        plot.Write("",ROOT.TObject.kOverwrite)

        outfile.Close()
    
    def write_mc(self,mc):

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        outfile_name = os.path.join(self.out_dir, "Analysis-root", self.config.Analysis["Ana_name"]+".root")

        outfile = ROOT.TFile(outfile_name, "UPDATE")
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            self.logger.info(f"Writing efficiency of {frame} framework...")
            Reco_prompt, Reco_nonprompt, Gen_prompt, Gen_nonprompt = self.mc_ops.get_mc_sparse(mc,frame)

            Reco_total = Reco_prompt.Clone("Reco_total")
            Reco_total.Add(Reco_nonprompt)

            Gen_total = Gen_prompt.Clone("Gen_total")
            Gen_total.Add(Gen_nonprompt)

            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            y_min = Gen_prompt.GetAxis(2).FindBin(-0.8+1e-6)
            y_max = Gen_prompt.GetAxis(2).FindBin(0.8-1e-6)
            Gen_prompt.GetAxis(2).SetRange(y_min,y_max)
            Gen_nonprompt.GetAxis(2).SetRange(y_min,y_max)
            Gen_total.GetAxis(2).SetRange(y_min,y_max)

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue
                self.logger.info(f"    Working in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                pt_bin_dir.cd()

                pt_min = Reco_prompt.GetAxis(1).FindBin(pt_min_edge+1e-6)
                pt_max = Reco_prompt.GetAxis(1).FindBin(pt_max_edge-1e-6)
                Reco_prompt.GetAxis(1).SetRange(pt_min,pt_max)
                Reco_nonprompt.GetAxis(1).SetRange(pt_min,pt_max)
                Gen_prompt.GetAxis(0).SetRange(pt_min,pt_max)
                Gen_nonprompt.GetAxis(0).SetRange(pt_min,pt_max)
                Reco_total.GetAxis(1).SetRange(pt_min,pt_max)
                Gen_total.GetAxis(0).SetRange(pt_min,pt_max)

                D0_min = Reco_prompt.GetAxis(4).FindBin(pt_bin_set["D0Mass"][0]+1e-6)
                D0_max = Reco_prompt.GetAxis(4).FindBin(pt_bin_set["D0Mass"][1]-1e-6)
                Reco_prompt.GetAxis(4).SetRange(D0_min,D0_max)
                Reco_nonprompt.GetAxis(4).SetRange(D0_min,D0_max)
                Reco_total.GetAxis(4).SetRange(D0_min,D0_max)

                bkg_max = Reco_prompt.GetAxis(6).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
                Reco_prompt.GetAxis(6).SetRange(0, bkg_max)
                Reco_nonprompt.GetAxis(6).SetRange(0, bkg_max)
                Reco_total.GetAxis(6).SetRange(0, bkg_max)
                
                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]

                fd_bin = np.array(fd_edges)
                
                hReco_prompt_fd = ROOT.TH1F("hReco_prompt_fd",f"hReco_prompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Entries",len(fd_bin)-1,fd_bin)
                hGen_prompt_fd = ROOT.TH1F("hGen_prompt_fd",f"hGen_prompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Entries",len(fd_bin)-1,fd_bin)
                hEff_prompt_fd = ROOT.TH1F("hEff_prompt_fd",f"hEff_prompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Efficiency",len(fd_bin)-1,fd_bin)

                hReco_nonprompt_fd = ROOT.TH1F("hReco_nonprompt_fd",f"hReco_nonprompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Entries",len(fd_bin)-1,fd_bin)
                hGen_nonprompt_fd = ROOT.TH1F("hGen_nonprompt_fd",f"hGen_nonprompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Entries",len(fd_bin)-1,fd_bin)
                hEff_nonprompt_fd = ROOT.TH1F("hEff_nonprompt_fd",f"hEff_nonprompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Efficiency",len(fd_bin)-1,fd_bin)
                
                hReco_total_fd = ROOT.TH1F("hReco_total_fd",f"hReco_total_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Entries",len(fd_bin)-1,fd_bin)
                hGen_total_fd = ROOT.TH1F("hGen_total_fd",f"hGen_total_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Entries",len(fd_bin)-1,fd_bin)
                hEff_total_fd = ROOT.TH1F("hEff_total_fd",f"hEff_total_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};FD;Efficiency",len(fd_bin)-1,fd_bin)

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    cos_max = Reco_prompt.GetAxis(5).FindBin(1-1e-6)
                    Reco_prompt.GetAxis(5).SetRange(0,cos_max)
                    Reco_nonprompt.GetAxis(5).SetRange(0,cos_max)
                    Gen_prompt.GetAxis(3).SetRange(0,cos_max)
                    Gen_nonprompt.GetAxis(3).SetRange(0,cos_max)
                    Reco_total.GetAxis(5).SetRange(0,cos_max)
                    Gen_total.GetAxis(3).SetRange(0,cos_max)

                    fd_min = Reco_prompt.GetAxis(7).FindBin(fd_min_edge+1e-6)
                    fd_max = Reco_prompt.GetAxis(7).FindBin(fd_max_edge-1e-6)
                    Reco_prompt.GetAxis(7).SetRange(fd_min,fd_max)
                    Reco_nonprompt.GetAxis(7).SetRange(fd_min,fd_max)
                    Reco_total.GetAxis(7).SetRange(fd_min,fd_max)

                    if len(fd_edges) > 2:
                        bin_num = ifd
                    else:
                        bin_num = ifd+1

                    hReco_prompt_fd = self.Entries_fill(hReco_prompt_fd,Reco_prompt,bin_num)
                    hReco_nonprompt_fd = self.Entries_fill(hReco_nonprompt_fd,Reco_nonprompt,bin_num)
                    hReco_total_fd =  self.Entries_fill(hReco_total_fd,Reco_total,bin_num)

                    hGen_prompt_fd = self.Entries_fill(hGen_prompt_fd,Gen_prompt,bin_num)
                    hGen_nonprompt_fd = self.Entries_fill(hGen_nonprompt_fd,Gen_nonprompt,bin_num)
                    hGen_total_fd = self.Entries_fill(hGen_total_fd,Gen_total,bin_num)

                    cos_bin = np.array(cos_edges)
                
                    hReco_prompt_cos = ROOT.TH1F("hReco_prompt_cos",f"hReco_prompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin)
                    hGen_prompt_cos = ROOT.TH1F("hGen_prompt_cos",f"hGen_prompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin)
                    hEff_prompt_cos = ROOT.TH1F("hEff_prompt_cos",f"hEff_prompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Efficiency",len(cos_bin)-1,cos_bin)

                    hReco_nonprompt_cos = ROOT.TH1F("hReco_nonprompt_cos",f"hReco_nonprompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin)
                    hGen_nonprompt_cos = ROOT.TH1F("hGen_nonprompt_cos",f"hGen_nonprompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin)
                    hEff_nonprompt_cos = ROOT.TH1F("hEff_nonprompt_cos",f"hEff_nonprompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Efficiency",len(cos_bin)-1,cos_bin)
             
                    hReco_total_cos = ROOT.TH1F("hReco_total_cos",f"hReco_total_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin)
                    hGen_total_cos = ROOT.TH1F("hGen_total_cos",f"hGen_total_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin)
                    hEff_total_cos = ROOT.TH1F("hEff_total_cos",f"hEff_total_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge};Cos#vartheta*;Efficiency",len(cos_bin)-1,cos_bin)


                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        cos_dir = fd_dir.mkdir(f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", "", ROOT.kTRUE)
                        cos_dir.cd()

                        cos_min = Reco_prompt.GetAxis(5).FindBin(cos_min_edge+1e-6)
                        cos_max = Reco_prompt.GetAxis(5).FindBin(cos_max_edge-1e-6)
                        Reco_prompt.GetAxis(5).SetRange(cos_min,cos_max)
                        Reco_nonprompt.GetAxis(5).SetRange(cos_min,cos_max)
                        Gen_prompt.GetAxis(3).SetRange(cos_min,cos_max)
                        Gen_nonprompt.GetAxis(3).SetRange(cos_min,cos_max)
                        Reco_total.GetAxis(5).SetRange(cos_min,cos_max)
                        Gen_total.GetAxis(3).SetRange(cos_min,cos_max)

                        hReco_prompt_cos = self.Entries_fill(hReco_prompt_cos,Reco_prompt,icos+1)
                        hReco_nonprompt_cos = self.Entries_fill(hReco_nonprompt_cos,Reco_nonprompt,icos+1)
                        hReco_total_cos = self.Entries_fill(hReco_total_cos,Reco_total,icos+1)

                        hGen_prompt_cos = self.Entries_fill(hGen_prompt_cos,Gen_prompt,icos+1)
                        hGen_nonprompt_cos = self.Entries_fill(hGen_nonprompt_cos,Gen_nonprompt,icos+1)
                        hGen_total_cos = self.Entries_fill(hGen_total_cos,Gen_total,icos+1)

                    fd_dir.cd()

                    hEff_prompt_cos.Divide(hReco_prompt_cos,hGen_prompt_cos,1,1,"B")
                    hEff_prompt_cos.Write("",ROOT.TObject.kOverwrite)

                    hEff_nonprompt_cos.Divide(hReco_nonprompt_cos,hGen_nonprompt_cos,1,1,"B")
                    hEff_nonprompt_cos.Write("",ROOT.TObject.kOverwrite)

                    hEff_total_cos.Divide(hReco_total_cos,hGen_total_cos,1,1,"B")
                    hEff_total_cos.Write("",ROOT.TObject.kOverwrite)

                    # Sligihtly different if calculate efficiency from cos's entries
                    # Reco_error = ctypes.c_double()
                    # hReco_Enties_fd.SetBinContent(ifd,Reco_mc.Projection(0).IntegralAndError(1,Reco_mc.Projection(0).GetNbinsX(),Reco_error))
                    # hReco_Enties_fd.SetBinError(ifd,Reco_error.value)

                    # Gen_error = ctypes.c_double()
                    # hGen_Enties_fd.SetBinContent(ifd,hGen_Enties_cos.IntegralAndError(1,hGen_Enties_cos.GetNbinsX(),Gen_error))
                    # hGen_Enties_fd.SetBinError(ifd,Gen_error.value)

                pt_bin_dir.cd()
              
                hEff_prompt_fd.Divide(hReco_prompt_fd,hGen_prompt_fd,1,1,"B")
                hEff_prompt_fd.Write("",ROOT.TObject.kOverwrite)

                hEff_nonprompt_fd.Divide(hReco_nonprompt_fd,hGen_nonprompt_fd,1,1,"B")
                hEff_nonprompt_fd.Write("",ROOT.TObject.kOverwrite)

                hEff_total_fd.Divide(hReco_total_fd,hGen_total_fd,1,1,"B")
                hEff_total_fd.Write("",ROOT.TObject.kOverwrite)

        outfile.Close()

    def read_frac(self):
        
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        outfile_name = os.path.join(self.out_dir, "Analysis-root", self.config.Analysis["Ana_name"]+".root")

        outfile = ROOT.TFile(outfile_name, "UPDATE")
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            self.logger.info(f"Reading fraction of {frame} framework...")
            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue
                self.logger.info(f"    Working in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                pt_bin_dir.cd()

                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]

                heff_prompt_fd = pt_bin_dir.Get("hEff_prompt_fd")
                heff_nonprompt_fd = pt_bin_dir.Get(f"hEff_nonprompt_fd")

                hfrac_np_fd = heff_prompt_fd.Clone("hfrac_np_fd")

                frac_file_dir = os.path.join(self.out_dir,f"Cut-variation/{self.config.Analysis['Ana_name']}/pt_{pt_min_edge:0d}_{pt_max_edge:0d}/fraction","CutVar_"+self.config.Analysis['Ana_name']+".root")
                frac_file = ROOT.TFile(frac_file_dir,"READ")
                corr_yield_prompt = frac_file.Get("hCorrYieldsPrompt")
                corr_yield_nonprompt = frac_file.Get("hCorrYieldsNonPrompt")
                Cov_prompt_nonprompt = frac_file.Get("hCovPromptNonPrompt")

                hEff_total = pt_bin_dir.Get("fd_0.00_1.00/hEff_total_cos")
                
                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    if len(fd_edges) > 2:
                        bin_num = ifd
                    else:
                        bin_num = ifd+1

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    raw_np_frc, raw_np_frc_error = self.get_nonprompt_frac(heff_prompt_fd.GetBinContent(bin_num),
                                                                heff_nonprompt_fd.GetBinContent(bin_num),
                                                                corr_yield_prompt.GetBinContent(1),
                                                                corr_yield_nonprompt.GetBinContent(1),
                                                                corr_yield_prompt.GetBinError(1),
                                                                corr_yield_nonprompt.GetBinError(1),
                                                                Cov_prompt_nonprompt.GetBinContent(1))
 
                    hfrac_np_fd.SetBinContent(bin_num, raw_np_frc)
                    hfrac_np_fd.SetBinError(bin_num, raw_np_frc_error)

                    hraw_yield = fd_dir.Get("hraw_yield")
                    hEff_cos = fd_dir.Get("hEff_total_cos")
                    # hp_eff_cos = fd_dir.Get("hEff_prompt_cos")
                    # hnp_eff_cos = fd_dir.Get("hEff_nonprompt_cos")
                    hfrac_cos = hraw_yield.Clone("hfrac_cos")

                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        hfrac_cos.SetBinContent(icos+1, raw_np_frc)
                        hfrac_cos.SetBinError(icos+1,raw_np_frc_error)
  
                    hcorr_yield_cos = hraw_yield.Clone("hcorr_yield_cos")
                    hcorr_yield_cos.Divide(hEff_cos)
                    hcorr_yield_cos.Write("hcorr_yield_cos_old",ROOT.TObject.kOverwrite)

                    # hEff_frac , hcorr_yield_frac = self.get_total_eff(hraw_yield,hfrac_cos,hp_eff_cos,hnp_eff_cos)
                    hcorr_yield_total = hraw_yield.Clone("hcorr_yield_total")
                    hcorr_yield_total.Divide(hEff_total)
                    hcorr_yield_total.GetYaxis().SetTitle("Corrected Yield")
                    hcorr_yield_total.SetTitle("Corrected Yield")

                    hEff_total.Write("hEff_total",ROOT.TObject.kOverwrite)
                    hcorr_yield_total.Write("hcorr_yield",ROOT.TObject.kOverwrite)

                pt_bin_dir.cd()
                hfrac_np_fd.Write("hfrac",ROOT.TObject.kOverwrite)

        outfile.Close()
        
    def Reduceing(self, file):

        reduced_file = file.replace(".root", "_reduced.root")
        
        infile_name = os.path.join(self.config.Directories["InputDir"], "Data/Primary", file)
        outfile_name = os.path.join(self.config.Directories["InputDir"], "Data/Reduced", reduced_file)

        axis_to_keep = self.config.Data_keep_axis
        
        infile = ROOT.TFile(infile_name)
        outfile = ROOT.TFile(outfile_name, "RECREATE")
        outfile.cd()
        for frame in tqdm(["Helicity", "Production", "Random", "Beam"]):
            sparse = infile.Get(f"task-polarisation-charm-hadrons/h{frame}")

            pt_edges = self.config.BinSet["pt_bin_edges"]

            for pt_min, pt_max in zip(pt_edges[:-1], pt_edges[1:]):
                pt_bin_min = sparse.GetAxis(1).FindBin(pt_min+1e-6)
                pt_bin_max = sparse.GetAxis(1).FindBin(pt_max-1e-6)
                sparse.GetAxis(1).SetRange(pt_bin_min, pt_bin_max)
                sparse.GetAxis(11).SetRange(1, 1)
                sparse_signal_reduced = sparse.Projection(len(axis_to_keep), axis_to_keep)
                sparse_signal_reduced.SetName(f"h{frame}_pt_{pt_min:.0f}_{pt_max:.0f}_data")
                sparse_signal_reduced.Write()
                sparse.GetAxis(11).SetRange(2, 2)
                sparse_bkg_reduced = sparse.Projection(len(axis_to_keep), axis_to_keep)
                sparse_bkg_reduced.SetName(f"h{frame}_pt_{pt_min:.0f}_{pt_max:.0f}_rotbkg")
                sparse_bkg_reduced.Write()

        outfile.Close()
        return outfile_name
 
    def get_sparse(self, file_list, name):
        
        sparse = []
        for file_dir in file_list:
            file = ROOT.TFile(file_dir, "READ")
            sparse.append(file.Get(name))
            file.Close()
      
        for i in range(1, len(sparse)):
            try:
                sparse[0].Add(sparse[i])
                self.logger.info(f"Merging {i}_{sparse[i].GetName()}")
                sparse[i].Delete()
            except:
                self.logger.error(f"Error when mergeing {name}")
                sys.exit(1)
            
        return sparse[0]
        
    def Compare(self,hbkg,hmass):

        hbkg.Sumw2()
        hmass.Sumw2()
        hbkg.Scale(1./hbkg.Integral())
        hmass.Scale(1./hmass.Integral())

        hbkg_for_norm = hbkg.Clone("hbkg_for_norm")
        numbers = hmass.Integral()

        hbkg_for_norm.Divide(hmass)

        func = ROOT.TF1("constant_func", constant_func, 0.16, 0.178, 1)
        
        hbkg_for_norm.Fit("constant_func", "R")

        hbkg.Scale(1./func.GetParameter(0))

        hbkg.Scale(numbers)
        hmass.Scale(numbers)    

        canvas = ROOT.TCanvas("data_v.s_bkg","data_v.s_bkg",800,600)
        canvas.cd()
        hmass.SetLineColor(ROOT.kBlack)
        hmass.SetMarkerColor(ROOT.kBlack)
        hmass.SetMarkerStyle(20)
        hmass.SetMarkerSize(0.5)
        hmass.Draw("PE")
        hbkg.SetLineColor(ROOT.kRed)
        hbkg.SetLineWidth(2)
        hbkg.Draw("HIST SAME")
        legend = ROOT.TLegend(0.7,0.7,0.9,0.9)

        legend.AddEntry(hmass,"Data","lep")
        legend.AddEntry(hbkg,"Rotated bkg","l")
        legend.Draw("SAME")

        canvas.Update()
        return canvas
    
    def Entries_fill(self, hist,sparse,ibin):

        error = ctypes.c_double()
        hist.SetBinContent(ibin,sparse.Projection(0).IntegralAndError(1,sparse.Projection(0).GetNbinsX(),error))
        hist.SetBinError(ibin,error.value)

        return hist

    def get_nonprompt_frac(self, effacc_p,effacc_np,corry_p,corry_np,cov_pp,cov_npnp,cov_pnp):
        
        rawy_p = effacc_p * corry_p
        rawy_np = effacc_np * corry_np
        f_np = rawy_np / (rawy_p + rawy_np)
    
        # derivatives of prompt fraction wrt corr yields
        d_p = (effacc_p * (rawy_p + rawy_np) - effacc_p**2 * corry_p) / (rawy_p + rawy_np) ** 2
        d_np = -effacc_np * rawy_p / (rawy_p + rawy_np) ** 2
    
        f_p_unc = np.sqrt(
            d_p**2 * cov_pp**2
            + d_np**2 * cov_npnp**2
            + 2 * d_p *d_np* cov_pnp
        )
        return f_np, f_p_unc

    def get_total_eff(self, hraw_yield,hfrac_np,hp_eff,hnp_eff):
        
        hfrac_p = hfrac_np.Clone("hfrac_p")
        for i in range(hfrac_p.GetNbinsX()):
            hfrac_p.SetBinContent(i+1, 1-hfrac_np.GetBinContent(i+1))
            hfrac_p.SetBinError(i+1, hfrac_np.GetBinError(i+1))

        hraw_prompt = hraw_yield.Clone("hraw_prompt")
        hraw_nonprompt = hraw_yield.Clone("hraw_nonprompt")
        hEff_total = hraw_yield.Clone("hEff_total")

        hraw_nonprompt.Multiply(hraw_nonprompt,hfrac_np,1,1)
        hraw_prompt.Multiply(hraw_prompt,hfrac_p,1,1)
       
        hcorr_prompt = hraw_prompt.Clone("hcorr_prompt")
        hcorr_nonprompt = hraw_nonprompt.Clone("hcorr_nonprompt")
        hcorr_prompt.Divide(hp_eff)
        hcorr_nonprompt.Divide(hnp_eff)

        hcorr_yield= hcorr_nonprompt.Clone("hcorr_yield")
        hcorr_yield.Add(hcorr_yield,hcorr_prompt,1,1)

        hEff_total.Divide(hEff_total,hcorr_yield,1,1,"B")

        # print("---------------------------------------")
        # print("non-prompt frac = " , hfrac_np.GetBinContent(1),"error = ",hfrac_np.GetBinError(1))
        # print("prompt fraction = ", hfrac_p.GetBinContent(1),"error = ",hfrac_p.GetBinError(1))
        # print("raw_yield = ",hraw_yield.GetBinContent(1),"error = ",hraw_yield.GetBinError(1)/hraw_yield.GetBinContent(1))
        # print("raw_prompt = ", hraw_prompt.GetBinContent(1),"error = ",hraw_prompt.GetBinError(1)/hraw_prompt.GetBinContent(1))
        # print("raw_nonprompt = ",hraw_nonprompt.GetBinContent(1),"error = ",hraw_nonprompt.GetBinError(1)/hraw_nonprompt.GetBinContent(1))
        # print("eff_prompt = ", hp_eff.GetBinContent(1))
        # print("eff_nonprompt = ", hnp_eff.GetBinContent(1))
        # print("corr_prompt = " ,hcorr_prompt.GetBinContent(1),"error = ",hcorr_prompt.GetBinError(1)/hcorr_prompt.GetBinContent(1))
        # print("corr_nonprompt = ",hcorr_nonprompt.GetBinContent(1),"error = ",hcorr_nonprompt.GetBinError(1)/hcorr_nonprompt.GetBinContent(1))
        # print("corr_yield = ", hcorr_yield.GetBinContent(1),"error = ",hcorr_yield.GetBinError(1)/hcorr_yield.GetBinContent(1))
        # print("total_eff = ", hEff_total.GetBinContent(1))
        # print("total yield = ",  hcorr_yield.Integral())

        return hEff_total, hcorr_yield
    
    def write_mc_pars(self,mc):

        pt_edges = self.config.BinSet["pt_bin_edges"] 
        os.makedirs(os.path.join(os.path.join(os.getcwd(),"Output","Mc_pars"), "Analysis-root"), exist_ok=True) 
        outfile_name = os.path.join(os.path.join(os.getcwd(),"Output","Mc_pars"), "Analysis-root", "Mc_pars.root")

        outfile = ROOT.TFile(outfile_name, "UPDATE")
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            self.logger.info(f"Writing MC inv-mass of {frame} framework...")
            Reco_prompt, Reco_nonprompt, Gen_prompt, Gen_nonprompt = self.mc_ops.get_mc_sparse(mc,frame)

            Reco_total = Reco_prompt.Clone("Reco_total")
            Reco_total.Add(Reco_nonprompt)

            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue
                self.logger.info(f"    Working in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                pt_bin_dir.cd()

                D0_min = Reco_prompt.GetAxis(4).FindBin(pt_bin_set["D0Mass"][0]+1e-6)
                D0_max = Reco_prompt.GetAxis(4).FindBin(pt_bin_set["D0Mass"][1]-1e-6)
                Reco_prompt.GetAxis(4).SetRange(D0_min,D0_max)
                Reco_nonprompt.GetAxis(4).SetRange(D0_min,D0_max)
                Reco_total.GetAxis(4).SetRange(D0_min,D0_max)

                pt_min = Reco_prompt.GetAxis(1).FindBin(pt_min_edge+1e-6)
                pt_max = Reco_prompt.GetAxis(1).FindBin(pt_max_edge-1e-6)
                Reco_prompt.GetAxis(1).SetRange(pt_min,pt_max)
                Reco_nonprompt.GetAxis(1).SetRange(pt_min,pt_max)
                Reco_total.GetAxis(1).SetRange(pt_min,pt_max)

                bkg_max = Reco_prompt.GetAxis(6).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
                Reco_prompt.GetAxis(6).SetRange(0, bkg_max)
                Reco_nonprompt.GetAxis(6).SetRange(0, bkg_max)
                Reco_total.GetAxis(6).SetRange(0, bkg_max)
                
                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    cos_max = Reco_prompt.GetAxis(5).FindBin(1-1e-6)
                    Reco_prompt.GetAxis(5).SetRange(0,cos_max)
                    Reco_nonprompt.GetAxis(5).SetRange(0,cos_max)
                    Reco_total.GetAxis(5).SetRange(0,cos_max)

                    fd_min = Reco_prompt.GetAxis(7).FindBin(fd_min_edge+1e-6)
                    fd_max = Reco_prompt.GetAxis(7).FindBin(fd_max_edge-1e-6)
                    Reco_prompt.GetAxis(7).SetRange(fd_min,fd_max)
                    Reco_nonprompt.GetAxis(7).SetRange(fd_min,fd_max)
                    Reco_total.GetAxis(7).SetRange(fd_min,fd_max)

                    hmass_prompt = Reco_prompt.Projection(0)
                    hmass_prompt.Write(f"hprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_prompt",ROOT.TObject.kOverwrite)
                    hmass_nonprompt = Reco_nonprompt.Projection(0)
                    hmass_nonprompt.Write(f"hnonprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_nonprompt",ROOT.TObject.kOverwrite)
                    hmass_total = Reco_total.Projection(0)
                    hmass_total.Write(f"htotal_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_total",ROOT.TObject.kOverwrite)
                
                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        cos_dir = fd_dir.mkdir(f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", "", ROOT.kTRUE)
                        cos_dir.cd()

                        cos_min = Reco_prompt.GetAxis(5).FindBin(cos_min_edge+1e-6)
                        cos_max = Reco_prompt.GetAxis(5).FindBin(cos_max_edge-1e-6)
                        Reco_prompt.GetAxis(5).SetRange(cos_min,cos_max)
                        Reco_nonprompt.GetAxis(5).SetRange(cos_min,cos_max)
                        Reco_total.GetAxis(5).SetRange(cos_min,cos_max)

                        hmass_prompt = Reco_prompt.Projection(0)
                        hmass_prompt.Write(f"hprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}_prompt",ROOT.TObject.kOverwrite)
                        hmass_nonprompt = Reco_nonprompt.Projection(0)
                        hmass_nonprompt.Write(f"hnonprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}_nonprompt",ROOT.TObject.kOverwrite)
                        hmass_total = Reco_total.Projection(0)
                        hmass_total.Write(f"htotal_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}_total",ROOT.TObject.kOverwrite)

        outfile.Close()