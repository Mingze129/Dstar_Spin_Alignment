import os
import sys
import numpy as np
import ctypes
import ROOT
ROOT.gROOT.SetBatch(True)

from .DataOps import DataOps
from .McOps import McOps
from .FitOps import FitOps

sys.path.append("../")
import frc_test_config as config


class FracOps(object):

    def __init__(self, config,logger):
        self.config = config
        self.out_dir = os.path.join(self.config.Directories["OutputDir"],self.config.Analysis["Task_Name"],"Cut-variation")
        self.data_ops = DataOps(self.config,logger)
        self.data, self.mc = self.data_ops.load_data()
        self.mc_ops = McOps(config)
        self.fit_ops = FitOps(config,logger)

        self.logger = logger
        self.logger.info("Fracion Opreator Init...")

    def get_input(self):
        
        pt_edges = self.config.BinSet["pt_bin_edges"]

        TReco_prompt, TReco_nonprompt, TGen_prompt, TGen_nonprompt = self.mc_ops.get_mc_sparse(self.mc, "Helicity")
 
        for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

            pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
            if not pt_bin_set["doing"]:
                continue
            print(f"Do cut-variation in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

            os.makedirs(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/raw_yield"), exist_ok=True)
            os.makedirs(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/efficiency"), exist_ok=True)
            os.makedirs(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/fraction"), exist_ok=True)

            Tdata = self.data_ops.get_sparse(self.data, f"hHelicity_pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}_data")
            Tbkg = self.data_ops.get_sparse(self.data, f"hHelicity_pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}_rotbkg")
            
            # nbin = Tdata.GetAxis(2).FindBin(0.99)
            # Tdata.GetAxis(2).SetRange(nbin,nbin)
            # Tbkg.GetAxis(2).SetRange(nbin,nbin)
            # nbin = TReco_prompt.GetAxis(5).FindBin(0.99)
            # TReco_prompt.GetAxis(5).SetRange(nbin,nbin)
            # TReco_nonprompt.GetAxis(5).SetRange(nbin,nbin)
            # TGen_prompt.GetAxis(3).SetRange(nbin,nbin)
            # TGen_nonprompt.GetAxis(3).SetRange(nbin,nbin)

            pt_min = TReco_prompt.GetAxis(1).FindBin(pt_min_edge+1e-6)
            pt_max = TReco_prompt.GetAxis(1).FindBin(pt_max_edge-1e-6)
            TReco_prompt.GetAxis(1).SetRange(pt_min, pt_max)
            TReco_nonprompt.GetAxis(1).SetRange(pt_min, pt_max)
            TGen_prompt.GetAxis(0).SetRange(pt_min, pt_max)
            TGen_nonprompt.GetAxis(0).SetRange(pt_min, pt_max)

            bkg_max = TReco_prompt.GetAxis(6).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
            TReco_prompt.GetAxis(6).SetRange(0, bkg_max)
            TReco_nonprompt.GetAxis(6).SetRange(0, bkg_max)

            bkg_max = Tdata.GetAxis(3).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
            Tdata.GetAxis(3).SetRange(0, bkg_max)
            Tbkg.GetAxis(3).SetRange(0, bkg_max)

            fd_edges = pt_bin_set["var_fd_range"]

            pars_dict = {}

            for icut, fd_min_edge in enumerate(fd_edges[:-1]):

                print(f"Cut-variation, icut = {icut}, fd_min_edge = {fd_min_edge:.2f}")
                try:
                    fd_min = Tdata.GetAxis(4).FindBin(fd_min_edge+1e-6)
                    fd_max = Tdata.GetAxis(4).FindBin(1.0-1e-6)
                    Tdata.GetAxis(4).SetRange(fd_min, fd_max)

                    raw_yield_file = ROOT.TFile(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/raw_yield/{icut}_raw_yield_{fd_min_edge:.3f}_{fd_edges[icut+1]:.3f}.root"),"RECREATE")
                    hmass = Tdata.Projection(0).Clone("hmass")
                    raw_yield_file.cd()
                    hmass.Write("",ROOT.TObject.kOverwrite)
                    
                    file = os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/raw_yield/{icut}_raw_yield_{fd_min_edge:.2f}_{fd_edges[icut+1]:.2f}.root")

                    if icut == 0 :
                        Tbkg.GetAxis(4).SetRange(fd_min, fd_max)
                        hbkg = Tbkg.Projection(0).Clone("hbkg")
                        raw_yield_file.cd()
                        hbkg.Write("",ROOT.TObject.kOverwrite)
                        raw_yield_file.Close()

                        bkg_dir = "hbkg"
                        task_name = f"{icut}_hbkg_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge:.2f}_{fd_edges[icut+1]:.2f}"
                        fit_set = {"signal_func": ["nosignal"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "fix_pars":[False],
                                    "out_dir": os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}")
                                }
                        
                        raw_yield, raw_yield_error, par_dict_bkg = self.fit_ops.fit_inv_mass(file, bkg_dir, task_name, fit_set)

                        data_dir = "hmass"
                        task_name = f"{icut}_hmass_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge:.2f}_{fd_edges[icut+1]:.2f}"

                        fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "fix_pars":[True, par_dict_bkg,"power", "c1", "c2", "c3"],
                                    "out_dir": os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}")
                                }
                        
                        raw_yield, raw_yield_error, pars_dict_data = self.fit_ops.fit_inv_mass(file, data_dir, task_name, fit_set)
                        pars_dict.update(pars_dict_data)

                        raw_yield_file = ROOT.TFile(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/raw_yield/{icut}_raw_yield_{fd_min_edge:.3f}_{fd_edges[icut+1]:.3f}.root"),"UPDATE")
                        raw_yield_file.cd()
                        hraw_yield = ROOT.TH1F("hraw_yield","",1,pt_min_edge,pt_max_edge)
                        hraw_yield.SetBinContent(1,raw_yield)
                        hraw_yield.SetBinError(1,raw_yield_error)
                        hraw_yield.Write("",ROOT.TObject.kOverwrite)
                        raw_yield_file.Close()

                    else:
                        raw_yield_file.Close()
                        data_dir = "hmass"
                        task_name = f"{icut}_hmass_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge:.2f}_{fd_edges[icut+1]:.2f}"

                        fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "fix_pars":[True, pars_dict,"nl","nr","alphal","alphar"],
                                    "out_dir": os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}")
                                }
                        
                        raw_yield, raw_yield_error, pars_dict_no = self.fit_ops.fit_inv_mass(file, data_dir, task_name, fit_set)
                        
                        raw_yield_file = ROOT.TFile(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/raw_yield/{icut}_raw_yield_{fd_min_edge:.3f}_{fd_edges[icut+1]:.3f}.root"),"UPDATE")
                        raw_yield_file.cd()
                        hraw_yield = ROOT.TH1F("hraw_yield","",1,pt_min_edge,pt_max_edge)
                        hraw_yield.SetBinContent(1,raw_yield)
                        hraw_yield.SetBinError(1,raw_yield_error)
                        hraw_yield.Write("",ROOT.TObject.kOverwrite)
                        raw_yield_file.Close()

                except Exception as e:
                    self.logger.error(f"Raw yield extraction error in {pt_min_edge:.0f}-{pt_max_edge:.0f} bin {icut} cut {fd_min_edge:.2f}-{fd_edges[icut+1]:.2f}:")
                    self.logger.error("-"*50)
                    self.logger.error(e)
                    self.logger.error("-"*50)

                try:
                    fd_Score_min = TReco_prompt.GetAxis(7).FindBin(fd_min_edge+1e-6)
                    fd_Score_max = TReco_prompt.GetAxis(7).FindBin(1.0-1e-6)
                    TReco_prompt.GetAxis(7).SetRange(fd_Score_min, fd_Score_max)
                    TReco_nonprompt.GetAxis(7).SetRange(fd_Score_min, fd_Score_max)

                    eff_file = ROOT.TFile(os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/efficiency/{icut}_eff_{fd_min_edge:.3f}_{fd_edges[icut+1]:.3f}.root"),"RECREATE")
                    
                    heff_prompt = ROOT.TH1F("heff_prompt","",1,pt_min_edge,pt_max_edge)
                    heff_nonprompt = ROOT.TH1F("heff_nonprompt","",1,pt_min_edge,pt_max_edge)
                    heff_prompt_gen = ROOT.TH1F("heff_prompt_gen","",1,pt_min_edge,pt_max_edge)
                    heff_nonprompt_gen = ROOT.TH1F("heff_nonprompt_gen","",1,pt_min_edge,pt_max_edge)

                    error = ctypes.c_double(0)
                    heff_prompt.SetBinContent(1,TReco_prompt.Projection(0).IntegralAndError(1,TReco_prompt.Projection(0).GetNbinsX(),error))
                    heff_prompt.SetBinError(1,error.value)
                    heff_nonprompt.SetBinContent(1,TReco_nonprompt.Projection(0).IntegralAndError(1,TReco_nonprompt.Projection(0).GetNbinsX(),error))
                    heff_nonprompt.SetBinError(1,error.value)
                    heff_prompt_gen.SetBinContent(1,TGen_prompt.Projection(0).IntegralAndError(1,TGen_prompt.Projection(0).GetNbinsX(),error))
                    heff_prompt_gen.SetBinError(1,error.value)
                    heff_nonprompt_gen.SetBinContent(1,TGen_nonprompt.Projection(0).IntegralAndError(1,TGen_nonprompt.Projection(0).GetNbinsX(),error))
                    heff_nonprompt_gen.SetBinError(1,error.value)

                    heff_prompt.Divide(heff_prompt,heff_prompt_gen,1,1,"B")
                    heff_nonprompt.Divide(heff_nonprompt,heff_nonprompt_gen,1,1,"B")

                    eff_file.cd()
                    heff_prompt.Write()
                    heff_nonprompt.Write()
                    eff_file.Close()

                except Exception as e:
                    self.logger.error(f"Efficience extraction error in {pt_min_edge:.0f}-{pt_max_edge:.0f} bin {icut} cut {fd_min_edge:.2f}-{fd_edges[icut+1]:.2f}:")
                    self.logger.error("-"*50)
                    self.logger.error(e)
                    self.logger.error("-"*50)
       
    def get_fraction(self):
        
        pass

