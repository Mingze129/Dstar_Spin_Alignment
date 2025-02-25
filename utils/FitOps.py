import os 
import sys
import numpy as np
from tqdm import tqdm
import ctypes
import zfit
import ROOT
ROOT.gROOT.SetBatch(True)

from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from particle import Particle


class FitOps(object):

    def __init__(self,config, logger):
        self.logger = logger
        self.config = config
        self.out_dir = os.path.join(self.config.Directories["OutputDir"], self.config.Analysis["Task_Name"])
        if config.Analysis["Cre_fit_root"]:
            fit_root = ROOT.TFile(os.path.join(os.path.join(self.out_dir, "Mass-Fit"), 'AnalysisFit.root'), "RECREATE")
            fit_root.Close()
        self.ax_title = r"$M(\mathrm{K\pi\pi}) - M(\mathrm{K\pi})$ GeV$/c^2$"
        self.func_pars = {
            "Init_value": {"sigma": 0.0005,
                           "frac": 0.1,
                           "nr": 20,
                           "nl":20,
                           "alphal": 1.5,
                           "alphar": 1.5
                           },
            "expopowext": ["power","c1","c2","c3"],
            "expopow": ["lam"],
            "doublecb": ["frac","mu","sigma","alphal","alphar","nl","nr"],
            "nosignal": []
        }    
        
    def get_raw_yield(self):

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "AnalysisSpinAlignment.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            self.logger.info(f"Start fitting raw-yield for {frame} frame...")

            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                self.logger.info(f"     Fitting raw-yield for pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}...")
                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue

                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = fd_edges[:-1]+[0.0]
                    fd_max_edges = fd_edges[1:]+[1.0]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]
                
                if not pt_bin_set["doing"]:
                    continue
                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    bkg_fd_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hrotbkg_{frame}_pt_{pt_min_edge}_{pt_max_edge}")
                    task_name = f"hrotbkg_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"

                    fit_set = {"signal_func": ["nosignal"],
                                "bkg_func": pt_bin_set["Bkg_func"],
                                "mass_range": pt_bin_set["Mass_range"],
                                "rebin": pt_bin_set["Rebin"],
                                "bin_counting": [True, 0.1396, 0.165],
                                "fix_pars":[False],
                                "out_dir": self.out_dir
                            }
                    raw_yield, raw_yield, par_dict_bkg = self.fit_inv_mass(outfile_name, bkg_fd_dir, task_name, fit_set)

                    data_fd_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}")
                    task_name = f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"

                    fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                "bkg_func": pt_bin_set["Bkg_func"],
                                "mass_range": pt_bin_set["Mass_range"],
                                "rebin": pt_bin_set["Rebin"],
                                "bin_counting": [True, 0.1396, 0.165],
                                "fix_pars":[True, par_dict_bkg,"power", "c1", "c2", "c3"],
                                "out_dir": self.out_dir
                            }
                    
                    raw_yield, raw_yield_error, par_dict_data = self.fit_inv_mass(outfile_name, data_fd_dir, task_name, fit_set)

                    cos_bin = np.array(cos_edges)
                    hraw_yield = ROOT.TH1F(f"hraw_yield", "hraw_yield_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge};Cos#vartheta*;raw_yield", len(cos_bin)-1, cos_bin)
                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        cos_dir = fd_dir.mkdir(f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", "", ROOT.kTRUE)
                        cos_dir.cd()

                        fit_task = f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}"
                        data_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}" , f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_cos_{cos_min_edge}_{cos_max_edge}")

                        fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": [True, 0.1396, 0.165],
                                    "fix_pars":[True, par_dict_data,"nl","nr","alphal","alphar"],
                                    "out_dir": self.out_dir
                                }
                        raw_yield , raw_yield_error, par_dict_cos = self.fit_inv_mass(outfile_name, data_dir,  fit_task, fit_set)
                        hraw_yield.SetBinContent(icos+1, raw_yield)
                        hraw_yield.SetBinError(icos+1, raw_yield_error)
                    fd_dir.cd()
                    hraw_yield.Write("", ROOT.TObject.kOverwrite)

        outfile.Close()
                 
    def fit_inv_mass(self, outfile_name, data_dir,  task_name , fit_set):

        os.makedirs(os.path.join(fit_set["out_dir"], "Mass-Fit/Figure"), exist_ok=True)
        figure_dir = os.path.join(fit_set["out_dir"], "Mass-Fit/Figure")
        fit_dir = os.path.join(fit_set["out_dir"], "Mass-Fit")

        data_hdl = DataHandler(data=outfile_name, 
                                histoname=data_dir,
                                limits=fit_set["mass_range"],
                                rebin=fit_set["rebin"])       

        fitter = F2MassFitter(data_hdl, name_signal_pdf=fit_set["signal_func"], signal_at_threshold=[False], 
                                name_background_pdf=fit_set["bkg_func"],
                                name=f"{task_name}_fit", 
                                chi2_loss=False,verbosity=7, tol=1.e-1)

        mass_init = (Particle.from_pdgid(413).mass - Particle.from_pdgid(421).mass)*1.e-3
        fitter.set_particle_mass(0, mass=mass_init, limits=[mass_init*0.95, mass_init*1.05])

        for par in self.func_pars[fit_set["bkg_func"][0]]:
            if par in self.func_pars["Init_value"]:
                fitter.set_background_initpar(0,par,self.func_pars["Init_value"][par])

        for par in self.func_pars[fit_set["signal_func"][0]]:
             if par in self.func_pars["Init_value"]:
                fitter.set_signal_initpar(0,par,self.func_pars["Init_value"][par])

        if fit_set["fix_pars"][0]:
            for par in fit_set["fix_pars"][2:]:
                if par in self.func_pars[fit_set["bkg_func"][0]]:
                    fitter.set_background_initpar(0, par, fit_set["fix_pars"][1][par],fix = True)
                elif par in self.func_pars[fit_set["signal_func"][0]]:
                    fitter.set_signal_initpar(0, par, fit_set["fix_pars"][1][par],fix = True)

        fit_result = fitter.mass_zfit()

        if fit_result.converged:
            try:
                fig , axs= fitter.plot_mass_fit(style="ATLAS", show_extra_info=True,
                                            figsize=(8, 8), extra_info_loc=["upper left", "right"],
                                            axis_title=self.ax_title)
                fig.savefig(os.path.join(figure_dir, f"{task_name}.pdf"))
                if not os.path.exists(os.path.join(fit_dir, 'AnalysisFit.root')):
                    fit_root = ROOT.TFile(os.path.join(fit_dir, 'AnalysisFit.root'), "RECREATE")
                    fit_root.Close()
                fitter.dump_to_root(filename = f"{os.path.join(fit_dir, 'AnalysisFit.root')}", option="update",
                                        suffix=f"_{task_name}_bkg")
            except:
                fig , axs= fitter.plot_mass_fit(style="ATLAS", show_extra_info=False,
                                            figsize=(8, 8), extra_info_loc=["upper left", "right"],
                                            axis_title=self.ax_title)
                fig.savefig(os.path.join(figure_dir, f"{task_name}.pdf"))
                if not os.path.exists(os.path.join(fit_dir, 'AnalysisFit.root')):
                    fit_root = ROOT.TFile(os.path.join(fit_dir, 'AnalysisFit.root'), "RECREATE")
                    fit_root.Close()
                fitter.dump_to_root(filename = f"{os.path.join(fit_dir, 'AnalysisFit.root')}", option="update",
                                        suffix=f"_{task_name}_bkg")
        else:
            self.logger.error(f"Fit for {task_name} is not converged!")

        par_dict = {}
        for par in self.func_pars[fit_set["bkg_func"][0]]:
            par_dict[par] = fitter.get_background_parameter(0, par)[0]
            par_dict[par+"_error"] = fitter.get_background_parameter(0, par)[1]
        
        for par in self.func_pars[fit_set["signal_func"][0]]:
            par_dict[par] = fitter.get_signal_parameter(0, par)[0]
            par_dict[par+"_error"] = fitter.get_signal_parameter(0, par)[1]

        if fit_set["bin_counting"][0]:
            raw_yield, raw_yield_err = fitter.get_raw_yield_bincounting(0, min=fit_set["bin_counting"][1], max=fit_set["bin_counting"][2])
        else:
            raw_yield, raw_yield_err = fitter.get_raw_yield(0)

        return raw_yield, raw_yield_err, par_dict
    