import os
import sys
import numpy as np
import re
import ctypes
import ROOT
ROOT.gROOT.SetBatch(True)

from .DataOps import DataOps
from .McOps import McOps
from .FitOps import FitOps


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

    def get_yield_input(self):
        
        pt_edges = self.config.BinSet["pt_bin_edges"]

        for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

            pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
            if not pt_bin_set["doing"]:
                continue
            self.logger.info(f"Get raw yield for cut-variation in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

            raw_dir = os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/{self.config.Analysis['Ana_name']}/raw_yield")
            frac_dir = os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/{self.config.Analysis['Ana_name']}/fraction")
            os.makedirs(raw_dir, exist_ok=True)
            os.makedirs(frac_dir, exist_ok=True)

            Tdata = self.data_ops.get_sparse(self.data, f"hHelicity_pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}_data")
            Tbkg = self.data_ops.get_sparse(self.data, f"hHelicity_pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}_rotbkg")

            Tdata.GetAxis(5).SetRange(self.config.Analysis["Min_eta_track"],100)
            Tbkg.GetAxis(5).SetRange(self.config.Analysis["Min_eta_track"],100)
            Tdata.GetAxis(6).SetRange(self.config.Analysis["Min_cls_ITS"],100)
            Tbkg.GetAxis(6).SetRange(self.config.Analysis["Min_cls_ITS"],100)
            Tdata.GetAxis(7).SetRange(self.config.Analysis["Min_cls_TPC"],100)
            Tbkg.GetAxis(7).SetRange(self.config.Analysis["Min_cls_TPC"],100)

            bkg_max = Tdata.GetAxis(3).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
            Tdata.GetAxis(3).SetRange(0, bkg_max)
            Tbkg.GetAxis(3).SetRange(0, bkg_max)

            fd_edges = pt_bin_set["var_fd_range"]

            pars_dict = {}

            for icut, fd_min_edge in enumerate(fd_edges[:-1]):

                cos_max = Tdata.GetAxis(2).FindBin(1-1e-6)
                Tdata.GetAxis(2).SetRange(0,cos_max)
                Tbkg.GetAxis(2).SetRange(0,cos_max)

                fd_min = Tdata.GetAxis(4).FindBin(fd_min_edge+1e-6)
                fd_max = Tdata.GetAxis(4).FindBin(1.0-1e-6)
                Tdata.GetAxis(4).SetRange(fd_min, fd_max)
                Tbkg.GetAxis(4).SetRange(fd_min, fd_max)

                self.logger.info(f"Cut-variation, icut = {icut}, fd_min_edge = {fd_min_edge:.3f}")
                try:
                    
                    raw_yield_file = ROOT.TFile(os.path.join(raw_dir, f"{icut}_raw_yield_fd-cut_{fd_min_edge:.3f}.root"),"RECREATE")
                    hmass = Tdata.Projection(0).Clone("hmass")
                    raw_yield_file.cd()
                    hmass.Write("",ROOT.TObject.kOverwrite)
                    
                    file = os.path.join(raw_dir, f"{icut}_raw_yield_fd-cut_{fd_min_edge:.3f}.root")

                    if icut == 0 :
                        
                        hbkg = Tbkg.Projection(0).Clone("hbkg")
                        raw_yield_file.cd()
                        hbkg.Write("",ROOT.TObject.kOverwrite)
                        raw_yield_file.Close()

                        bkg_dir = "hbkg"
                        task_name = f"{icut}_hbkg_pt_{pt_min_edge}_{pt_max_edge}_fd-cut_{fd_min_edge:.3f}"
                        fit_set = {"signal_func": ["nosignal"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "init_pars":[False],
                                    "fix_pars":[False],
                                    "out_dir": os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}",self.config.Analysis["Ana_name"]) 
                                }
                        
                        raw_yield, raw_yield_error, par_dict_bkg = self.fit_ops.fit_inv_mass(file, bkg_dir, task_name, fit_set)

                        data_dir = "hmass"
                        task_name = f"{icut}_hmass_pt_{pt_min_edge}_{pt_max_edge}_fd-cut_{fd_min_edge:.3f}"

                        fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "init_pars":[False,par_dict_bkg],
                                    "fix_pars":[True, "power", "c1", "c2", "c3"],
                                    "out_dir": os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}",self.config.Analysis["Ana_name"])
                                }
                        
                        raw_yield, raw_yield_error, pars_dict_data = self.fit_ops.fit_inv_mass(file, data_dir, task_name, fit_set)
                        pars_dict.update(pars_dict_data)

                        raw_yield_file = ROOT.TFile(os.path.join(raw_dir, f"{icut}_raw_yield_fd-cut_{fd_min_edge:.3f}.root"),"UPDATE")
                        raw_yield_file.cd()
                        hraw_yield = ROOT.TH1F("hraw_yield","",1,pt_min_edge,pt_max_edge)
                        hraw_yield.SetBinContent(1,raw_yield)
                        hraw_yield.SetBinError(1,raw_yield_error)
                        hraw_yield.Write("",ROOT.TObject.kOverwrite)
                        raw_yield_file.Close()

                    else:
                        raw_yield_file.Close()
                        data_dir = "hmass"
                        task_name = f"{icut}_hmass_pt_{pt_min_edge}_{pt_max_edge}_fd-cut_{fd_min_edge:.3f}"

                        fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "init_pars":[False,pars_dict],
                                    "fix_pars":[True, "nl","nr","alphal","alphar","sigma"],
                                    "out_dir": os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}",self.config.Analysis["Ana_name"])
                                }
                        
                        raw_yield, raw_yield_error, pars_dict_no = self.fit_ops.fit_inv_mass(file, data_dir, task_name, fit_set)
                        
                        raw_yield_file = ROOT.TFile(os.path.join(raw_dir, f"{icut}_raw_yield_fd-cut_{fd_min_edge:.3f}.root"),"UPDATE")
                        raw_yield_file.cd()
                        hraw_yield = ROOT.TH1F("hraw_yield","",1,pt_min_edge,pt_max_edge)
                        hraw_yield.SetBinContent(1,raw_yield)
                        hraw_yield.SetBinError(1,raw_yield_error)
                        hraw_yield.Write("",ROOT.TObject.kOverwrite)
                        raw_yield_file.Close()

                except Exception as e:
                    self.logger.error(f"Raw yield extraction error in {pt_min_edge:.0f}-{pt_max_edge:.0f} bin {icut} cut {fd_min_edge:.3f}-{fd_edges[icut+1]:.3f}:")
                    self.logger.error("-"*50)
                    self.logger.error(e)
                    self.logger.error("-"*50)

    def get_eff_input(self):
        
        pt_edges = self.config.BinSet["pt_bin_edges"]

        TReco_prompt, TReco_nonprompt, TGen_prompt, TGen_nonprompt = self.mc_ops.get_mc_sparse(self.mc, "Helicity")
 
        for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

            pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
            if not pt_bin_set["doing"]:
                continue
            self.logger.info(f"Get efficiency for cut-variation in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

            eff_dir = os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/{self.config.Analysis['Ana_name']}/efficiency")
            frac_dir = os.path.join(self.out_dir,f"pt_{pt_min_edge:0d}_{pt_max_edge:0d}/{self.config.Analysis['Ana_name']}/fraction")
            os.makedirs(eff_dir, exist_ok=True)
            os.makedirs(frac_dir, exist_ok=True)

            pt_min = TReco_prompt.GetAxis(1).FindBin(pt_min_edge+1e-6)
            pt_max = TReco_prompt.GetAxis(1).FindBin(pt_max_edge-1e-6)
            TReco_prompt.GetAxis(1).SetRange(pt_min, pt_max)
            TReco_nonprompt.GetAxis(1).SetRange(pt_min, pt_max)
            TGen_prompt.GetAxis(0).SetRange(pt_min, pt_max)
            TGen_nonprompt.GetAxis(0).SetRange(pt_min, pt_max)

            y_min = TGen_prompt.GetAxis(2).FindBin(-0.8+1e-6)
            y_max = TGen_prompt.GetAxis(2).FindBin(0.8-1e-6)
            TGen_prompt.GetAxis(2).SetRange(y_min,y_max)
            TGen_nonprompt.GetAxis(2).SetRange(y_min,y_max)

            bkg_max = TReco_prompt.GetAxis(6).FindBin(pt_bin_set["Bkg_cut"]-1e-6)
            TReco_prompt.GetAxis(6).SetRange(0, bkg_max)
            TReco_nonprompt.GetAxis(6).SetRange(0, bkg_max)

            fd_edges = pt_bin_set["var_fd_range"]

            for icut, fd_min_edge in enumerate(fd_edges[:-1]):

                self.logger.info(f"Cut-variation, icut = {icut}, fd_min_edge = {fd_min_edge:.3f}")
               
                try:
                    fd_Score_min = TReco_prompt.GetAxis(7).FindBin(fd_min_edge+1e-6)
                    fd_Score_max = TReco_prompt.GetAxis(7).FindBin(1.0-1e-6)
                    TReco_prompt.GetAxis(7).SetRange(fd_Score_min, fd_Score_max)
                    TReco_nonprompt.GetAxis(7).SetRange(fd_Score_min, fd_Score_max)

                    eff_file = ROOT.TFile(os.path.join(eff_dir, f"{icut}_efficiency_fd-cut_{fd_min_edge:.3f}.root"),"RECREATE")
                    
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
                    self.logger.error(f"Efficience extraction error in {pt_min_edge:.0f}-{pt_max_edge:.0f} bin {icut} cut fd-cut_{fd_min_edge:.3f}:")
                    self.logger.error("-"*50)
                    self.logger.error(e)
                    self.logger.error("-"*50)
       
    
    def get_fraction(self):
        
        self.logger.info("Get fraction from raw yield and efficiency...")

        pt_min_set = re.compile(r"pt_min\s*=\s*\d+\n")    
        pt_max_set = re.compile(r"pt_max\s*=\s*\d+\n")
        start_bin_set = re.compile(r"start_bin\s*=\s*\d+\n")
        end_bin_set = re.compile(r"end_bin\s*=\s*\d+\n")
        remove_bin_set = re.compile(r"remove_bin\s*=\s*\[.*\]\n")
        config_file_set = re.compile(r"import\s*.*\s*as\s*config\n")
        

        pt_edges = self.config.BinSet["pt_bin_edges"]
        for ipt , (pt_min_edge, pt_max_edge) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
             
            pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
            if not pt_bin_set["doing"]:
                continue
            self.logger.info(f"Running in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

            cutvar_config_file = os.path.join(self.config.working_dir,"cut-variation/config_cutvar.py")
            cutvar_config = open(cutvar_config_file,"r+")
            config_content = cutvar_config.read()

            try:
                pt_min = pt_min_set.search(config_content).group()[:-1]
                pt_max = pt_max_set.search(config_content).group()[:-1]
                start_bin = start_bin_set.search(config_content).group()[:-1]
                end_bin = end_bin_set.search(config_content).group()[:-1]
                remove_bin = remove_bin_set.search(config_content).group()[:-1]
                config_file = config_file_set.search(config_content).group()[:-1]
            except:
                self.logger.error(f"Error in {pt_min_edge:.0f}-{pt_max_edge:.0f} bin:")
                self.logger.error("-"*50)
                self.logger.error("Can not find the pt_min, pt_max, start_bin, end_bin, remove_bin in config_cutvar.py")
                self.logger.error("-"*50)
                continue

            file_name = self.config.__file__.split("/")[-1].split(".")[0]
            new_content = config_content.replace(pt_min,f"pt_min = {pt_min_edge:.0f}")
            new_content = new_content.replace(pt_max,f"pt_max = {pt_max_edge:.0f}")
            new_content = new_content.replace(start_bin,f"start_bin = {pt_bin_set['frac_min_bin']}")
            new_content = new_content.replace(end_bin,f"end_bin = {pt_bin_set['frac_max_bin']}")
            new_content = new_content.replace(remove_bin,f"remove_bin = {pt_bin_set['frac_remove_bin']}")
            new_content = new_content.replace(config_file,f"import {file_name} as config")

            cutvar_config.seek(0)
            cutvar_config.write(new_content)
            cutvar_config.truncate()
            cutvar_config.close()

            try:
                os.chdir(os.path.join(self.config.working_dir,"cut-variation"))
                os.system(f"python3 compute_fraction_cutvar.py")
                os.chdir(self.config.working_dir)
            except Exception as e:
                self.logger.error(f"Fraction extraction error in {pt_min_edge:.0f}-{pt_max_edge:.0f} bin:")
                self.logger.error("-"*50)
                self.logger.error(e)
                self.logger.error("-"*50)
            