import os 
import sys
import gc
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
        if config.Cre_fit_root:
            self.fit_file = os.path.join(os.path.join(self.out_dir, "Mass-Fit"), self.config.Analysis["Ana_name"]+'_Fit.root')
            fit_root = ROOT.TFile(self.fit_file, "RECREATE")
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
            "nosignal": [],
            "nobkg": []
        }    
        
    def get_raw_yield(self):

        outfile_name = os.path.join(self.out_dir, "Analysis-root", self.config.Analysis["Ana_name"]+".root")
        infile_name = os.path.join(self.out_dir, "Analysis-root", "Data_And_Efficiency.root")
        outfile_name = os.path.join(self.out_dir, "Analysis-root", "RawYield_Extraction.root")

        infile = ROOT.TFile(infile_name, "READ")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        in_ana_dir = infile.Get(self.config.Analysis["Ana_name"])
        ana_dir = outfile.mkdir(self.config.Analysis["Ana_name"],"",ROOT.kTRUE)
        ana_dir.cd()

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]
        ratio_pars = ["sigma","alphal","alphar"]

        # mc_pars_file_dir = os.path.join(os.path.join(os.getcwd(),"Output","Mc_pars"), "Analysis-root", "Mc_pars.root")
        # mc_pars_file = ROOT.TFile(mc_pars_file_dir, "READ")
        for frame in frame_list:

            self.logger.info(f"Start fitting raw-yield for {frame} frame...")

            in_frame_dir = in_ana_dir.Get(frame)
            type_dir = ana_dir.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):
                gc.collect()
                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue
                self.logger.info(f"     Fitting raw-yield for pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}...")

                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]
                
                if not pt_bin_set["doing"]:
                    continue
                in_pt_dir = in_frame_dir.Get(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}")
                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)

                fd_bin_pars = {}
                cos_dict_list = []
                # fd_simga_pars = mc_pars_file.Get(f"{frame}/pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}/sigma_total_fd")
                # fd_alphal_pars = mc_pars_file.Get(f"{frame}/pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}/alphal_total_fd")
                # fd_alphar_pars = mc_pars_file.Get(f"{frame}/pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}/alphar_total_fd")

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    if pt_bin_set["with_bkg"]:
                        bkg_fd_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hrotbkg_{frame}_pt_{pt_min_edge}_{pt_max_edge}")
                        task_name = f"hrotbkg_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"

                        fit_set = {"signal_func": ["nosignal"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "chi2_loss": pt_bin_set["chi2_loss"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": pt_bin_set["bin_counting"],
                                    "init_pars": [False],
                                    "fix_pars":[False],
                                    "Custom_pars":[False],
                                    "threshold": pt_bin_set["threshold"],
                                    "corr_bkg": [False],
                                    "out_dir": self.out_dir
                                }
                        raw_yield, raw_yield, par_dict_bkg = self.fit_inv_mass(outfile_name, bkg_fd_dir, task_name, fit_set)
                        if pt_bin_set["Bkg_func"] == ["expopowext"]:
                            pars_set = [True,"power", "c1", "c2", "c3"]
                        elif pt_bin_set["Bkg_func"] == ["expopow"]:
                            pars_set = [True,"lam"]
                        init_set = [False,par_dict_bkg]
                    elif ifd == 0:
                        par_dict_bkg = {}
                        init_set = [False,par_dict_bkg]
                        pars_set = [False]

                    else:
                        par_dict_bkg = fd_bin_pars
                        init_set = [False,par_dict_bkg]
                        pars_set = pt_bin_set["fix_pars"]

                    data_fd_dir = os.path.join(self.config.Analysis["Ana_name"], frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}")
                    task_name = f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"
                    if pt_bin_set["corr_bkg"][0]:
                        factor_hist = in_pt_dir.Get("hist_Norm_cost")
                        fd_factor = factor_hist.GetBinContent(1)
                        corr_set = pt_bin_set["corr_bkg"] + [infile_name] + [os.path.join(self.config.Analysis["Ana_name"], frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hist_corrbkg_integral_cos")] + [fd_factor]
                    else:
                        corr_set = [False]

                    fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                "bkg_func": pt_bin_set["Bkg_func"],
                                "chi2_loss": pt_bin_set["chi2_loss"],
                                "mass_range": pt_bin_set["Mass_range"],
                                "rebin": pt_bin_set["Rebin"],
                                "bin_counting": pt_bin_set["bin_counting"],
                                "init_pars": init_set,
                                "fix_pars":pars_set,
                                "Custom_pars":[False],
                                "threshold": pt_bin_set["threshold"],
                                "corr_bkg": corr_set,
                                "out_dir": self.out_dir
                            }
                    
                    raw_yield, raw_yield_error, par_dict_data = self.fit_inv_mass(infile_name, data_fd_dir, task_name, fit_set)
                    
                    if ifd == 0:
                        fd_bin_pars.update(par_dict_data)

                    cos_bin = np.array(cos_edges)
                    hraw_yield = ROOT.TH1F(f"hraw_yield", "hraw_yield_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge};Cos#vartheta*;raw_yield", len(cos_bin)-1, cos_bin)
                    
                    # cos_sigma_pars = mc_pars_file.Get(f"{frame}/pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}/fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}/sigma_total_cos")
                    # cos_alphal_pars = mc_pars_file.Get(f"{frame}/pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}/fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}/alphal_total_cos")
                    # cos_alphar_pars = mc_pars_file.Get(f"{frame}/pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}/fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}/alphar_total_cos")

                    # cos_sigma_pars.Scale(1/fd_simga_pars.GetBinContent(ifd+1))
                    # cos_alphal_pars.Scale(1/fd_alphal_pars.GetBinContent(ifd+1))
                    # cos_alphar_pars.Scale(1/fd_alphar_pars.GetBinContent(ifd+1))

                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        fit_task = f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}"
                        data_dir = os.path.join(self.config.Analysis["Ana_name"], frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}" , f"hmass_{frame}_pt_{pt_min_edge}_{pt_max_edge}_cos_{cos_min_edge}_{cos_max_edge}")

                        if ifd == 0:
                            init_pars_set = [False, par_dict_data]
                            fix_pars_set = pt_bin_set["fix_pars"]
                        else:
                            init_pars_set = [False, cos_dict_list[icos]]
                            fix_pars_set = [True,"nl","nr","alphal","alphar","sigma"]

                        # for par in par_dict_data:
                        #     cos_dict_pars[par] = par_dict_data[par]
                        
                        # cos_dict_pars["sigma"] = cos_sigma_pars.GetBinContent(icos+1)*par_dict_data["sigma"]
                        # cos_dict_pars["alphal"] = cos_alphal_pars.GetBinContent(icos+1)*par_dict_data["alphal"]
                        # cos_dict_pars["alphar"] = cos_alphar_pars.GetBinContent(icos+1)*par_dict_data["alphar"]
                        if pt_bin_set["corr_bkg"][0]:
                            factor_func = in_pt_dir.Get("decay_plateau_func")
                            cos_bin_center = (cos_min_edge+cos_max_edge)/2
                            fd_factor = factor_func.Eval(cos_bin_center)
                            corr_set = pt_bin_set["corr_bkg"] + [infile_name] + [os.path.join(self.config.Analysis["Ana_name"], frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", f"hist_corrbkg_template")] + [fd_factor]
                        else:
                            corr_set = [False]


                        fit_set = {"signal_func": pt_bin_set["Signal_func"],
                                    "bkg_func": pt_bin_set["Bkg_func"],
                                    "chi2_loss": pt_bin_set["chi2_loss"],
                                    "mass_range": pt_bin_set["Mass_range"],
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": pt_bin_set["bin_counting"],
                                    "init_pars": init_pars_set,
                                    "fix_pars": fix_pars_set,
                                    "Custom_pars":[False],
                                    "threshold": pt_bin_set["threshold"],
                                    "corr_bkg": corr_set,
                                    "out_dir": self.out_dir
                                }
                        raw_yield , raw_yield_error, par_dict_cos = self.fit_inv_mass(infile_name, data_dir,  fit_task, fit_set)
                        if ifd == 0:
                            cos_dict_list.append(par_dict_cos)

                        hraw_yield.SetBinContent(icos+1, raw_yield)
                        hraw_yield.SetBinError(icos+1, raw_yield_error)
                    fd_dir.cd()
                    hraw_yield.Write("", ROOT.TObject.kOverwrite)

        outfile.Close()
                 
    def fit_inv_mass(self, infile_name, data_dir,  task_name , fit_set):

        figure_dir = os.path.join(fit_set["out_dir"], "Mass-Fit", self.config.Analysis["Ana_name"]+"-Figure")
        fit_dir = os.path.join(fit_set["out_dir"], "Mass-Fit")
        os.makedirs(figure_dir, exist_ok=True)
        os.makedirs(fit_dir, exist_ok=True)

        if fit_set["corr_bkg"][0]:
            corr_bkg_hdl = DataHandler(data=fit_set["corr_bkg"][1], 
                                        histoname=fit_set["corr_bkg"][2],
                                        limits=fit_set["mass_range"],
                                        rebin=fit_set["rebin"])

        data_hdl = DataHandler(data=infile_name, 
                                histoname=data_dir,
                                limits=fit_set["mass_range"],
                                rebin=fit_set["rebin"])       

        if fit_set["corr_bkg"][0]:
            fitter = F2MassFitter(data_hdl, name_signal_pdf=fit_set["signal_func"], signal_at_threshold=fit_set["threshold"],
                                name_background_pdf=["hist"]+fit_set["bkg_func"],
                                name=f"{task_name}_fit", 
                                chi2_loss=fit_set["chi2_loss"],verbosity=7, tol=1.e-1,
                                corr_bkg_hdl=corr_bkg_hdl)
            fitter.set_background_template(0,corr_bkg_hdl)
            fitter.fix_bkg_frac_to_signal_pdf(0,0,fit_set["corr_bkg"][3])
        else:
            fitter = F2MassFitter(data_hdl, name_signal_pdf=fit_set["signal_func"], signal_at_threshold=fit_set["threshold"], 
                                name_background_pdf=fit_set["bkg_func"],
                                name=f"{task_name}_fit", 
                                chi2_loss=fit_set["chi2_loss"],verbosity=7, tol=1.e-1)

        mass_init = (Particle.from_pdgid(413).mass - Particle.from_pdgid(421).mass)*1.e-3
        fitter.set_particle_mass(0, mass=mass_init, limits=[mass_init*0.95, mass_init*1.05])

        if fit_set["init_pars"][0]:
            for par in self.func_pars[fit_set["bkg_func"][0]]:
                if par in self.func_pars["Init_value"]:
                    par_value = fit_set["init_pars"][1][par]
                    par_error = np.abs(fit_set["init_pars"][1][par+"_error"])
                    if fit_set["corr_bkg"][0]:
                        fitter.set_background_initpar(1,par,par_value,limits=[par_value-par_error-abs(par_value)*0.2,par_value+par_error+abs(par_value)*0.2])
                    else:
                        fitter.set_background_initpar(0,par,par_value,limits=[par_value-par_error-abs(par_value)*0.2,par_value+par_error+abs(par_value)*0.2])

            for par in self.func_pars[fit_set["signal_func"][0]] and par in fit_set["init_pars"][1]:
                if par in self.func_pars["Init_value"]:
                    par_value = fit_set["init_pars"][1][par]
                    par_error = np.abs(fit_set["init_pars"][1][par+"_error"])
                    fitter.set_signal_initpar(0,par,par_value,limits=[par_value-par_error-abs(par_value)*0.2,par_value+par_error+abs(par_value)*0.2])
        else:
            
            for par in self.func_pars[fit_set["signal_func"][0]]:
                if par in self.func_pars["Init_value"]:
                    fitter.set_signal_initpar(0,par,self.func_pars["Init_value"][par])

            for par in self.func_pars[fit_set["bkg_func"][0]]:
                if par in self.func_pars["Init_value"]:
                    if fit_set["corr_bkg"][0]:
                        fitter.set_background_initpar(1,par,self.func_pars["Init_value"][par])
                    else:
                        fitter.set_background_initpar(0,par,self.func_pars["Init_value"][par])

        if fit_set["fix_pars"][0]:
            for par in fit_set["fix_pars"][1:]:
                if par in self.func_pars[fit_set["bkg_func"][0]] and par in fit_set["init_pars"][1]:
                    if fit_set["corr_bkg"][0]:
                        fitter.set_background_initpar(1, par, fit_set["init_pars"][1][par],fix = True)
                    else:
                        fitter.set_background_initpar(0, par, fit_set["init_pars"][1][par],fix = True)
                elif par in self.func_pars[fit_set["signal_func"][0]]:
                    if par in fit_set["init_pars"][1]:
                        fitter.set_signal_initpar(0, par, fit_set["init_pars"][1][par],fix = True)

        if fit_set["Custom_pars"][0]:
            for par in fit_set["Custom_pars"][1]:
                if par in self.func_pars[fit_set["bkg_func"][0]]:
                    if fit_set["corr_bkg"][0]:
                        fitter.set_background_initpar(1, par, fit_set["Custom_pars"][1][par], limits=[fit_set["Custom_pars"][1][par]-fit_set["Custom_pars"][1][par+"_err"],fit_set["Custom_pars"][1][par]+fit_set["Custom_pars"][1][par+"_err"]])
                    else:
                        fitter.set_background_initpar(0, par, fit_set["Custom_pars"][1][par], limits=[fit_set["Custom_pars"][1][par]-fit_set["Custom_pars"][1][par+"_err"],fit_set["Custom_pars"][1][par]+fit_set["Custom_pars"][1][par+"_err"]])
                elif par in self.func_pars[fit_set["signal_func"][0]]:
                    if par in fit_set["init_pars"][1]:
                        fitter.set_signal_initpar(0, par, fit_set["Custom_pars"][1][par], limits=[fit_set["Custom_pars"][1][par]-fit_set["Custom_pars"][1][par+"_err"],fit_set["Custom_pars"][1][par]+fit_set["Custom_pars"][1][par+"_err"]])
        
        try:
            fit_result = fitter.mass_zfit()
        except Exception as e:
            self.logger.error(f"Failed to fit {task_name}!")
            self.logger.error(f"Error: {e}")
            return 0, 0, {}

        if fit_result.converged:
            try:
                fig , axs= fitter.plot_mass_fit(style="ATLAS", show_extra_info=True,
                                            figsize=(8, 8), extra_info_loc=["upper left", "right"],
                                            axis_title=self.ax_title)
                fig.savefig(os.path.join(figure_dir, f"{task_name}.pdf"))
                if not os.path.exists(os.path.join(fit_dir, self.config.Analysis["Ana_name"]+'_Fit.root')):
                    fit_root = ROOT.TFile(os.path.join(fit_dir, self.config.Analysis["Ana_name"]+'_Fit.root'), "RECREATE")
                    fit_root.Close()
                fitter.dump_to_root(filename = f"{os.path.join(fit_dir, self.fit_file)}", option="update",
                                        suffix=f"_{task_name}")
            except:
                fig , axs= fitter.plot_mass_fit(style="ATLAS", show_extra_info=False,
                                            figsize=(8, 8), extra_info_loc=["upper left", "right"],
                                            axis_title=self.ax_title)
                fig.savefig(os.path.join(figure_dir, f"{task_name}.pdf"))
                if not os.path.exists(os.path.join(fit_dir, self.config.Analysis["Ana_name"]+'_Fit.root')):
                    fit_root = ROOT.TFile(os.path.join(fit_dir, self.config.Analysis["Ana_name"]+'_Fit.root'), "RECREATE")
                    fit_root.Close()
                try:
                    fitter.dump_to_root(filename = f"{os.path.join(fit_dir, self.fit_file)}", option="update",
                                        suffix=f"_{task_name}")
                except:
                    self.logger.error(f"Failed to save fit results for {task_name}!")
        else:
            try:
                fig , axs= fitter.plot_mass_fit(style="ATLAS", show_extra_info=False,
                                                figsize=(8, 8), extra_info_loc=["upper left", "right"],
                                                axis_title=self.ax_title)
                fig.savefig(os.path.join(figure_dir, f"{task_name}.pdf"))
                if not os.path.exists(os.path.join(fit_dir, self.config.Analysis["Ana_name"]+'_Fit.root')):
                    fit_root = ROOT.TFile(os.path.join(fit_dir, self.config.Analysis["Ana_name"]+'_Fit.root'), "RECREATE")
                    fit_root.Close()
                self.logger.error(f"Fit for {task_name} is not converged!")
            except:
                self.logger.error(f"Fit for {task_name} is not converged! and can't save figure!")

        par_dict = {}
        for par in self.func_pars[fit_set["bkg_func"][0]]:
            if fit_set["corr_bkg"][0] :
                par_dict[par] = fitter.get_background_parameter(1, par)[0]
                par_dict[par+"_error"] = fitter.get_background_parameter(1, par)[1]
            else:
                par_dict[par] = fitter.get_background_parameter(0, par)[0]
                par_dict[par+"_error"] = fitter.get_background_parameter(0, par)[1]
        
        for par in self.func_pars[fit_set["signal_func"][0]]:
            if fit_set["bkg_func"][0] == "nobkg" and par == "frac":
                continue
            par_dict[par] = fitter.get_signal_parameter(0, par)[0]
            par_dict[par+"_error"] = fitter.get_signal_parameter(0, par)[1]

        if fit_set["bin_counting"][0]:
            raw_yield, raw_yield_err = fitter.get_raw_yield_bincounting(0, min=fit_set["bin_counting"][1], max=fit_set["bin_counting"][2])
        else:
            raw_yield, raw_yield_err = fitter.get_raw_yield(0)

        return raw_yield, raw_yield_err, par_dict
    
    def get_mc_pars(self):
                
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        outfile_name = os.path.join(os.path.join(os.getcwd(),"Output","Mc_pars"), "Analysis-root", "Mc_pars.root")

        outfile = ROOT.TFile(outfile_name, "UPDATE")
        frame_list = self.config.Analysis["Framework"]

        save_pars = ["sigma", "nl", "nr", "mu", "alphal", "alphar"]
        mass_range = [0.1396, 0.16]

        for frame in frame_list:

            self.logger.info(f"Fitting MC inv-mass at {frame} framework...")
           
            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):
                gc.collect()
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

                fd_bin = np.array(fd_edges)
                
                hpars_prompt_fd = []
                hpars_nonprompt_fd = []
                hpars_total_fd = []

                for par in save_pars:
                    hpars_prompt_fd.append(ROOT.TH1F(f"{par}_prompt_fd",f"{par}_prompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};fd;Entries",len(fd_bin)-1,fd_bin))
                    hpars_nonprompt_fd.append(ROOT.TH1F(f"{par}_nonprompt_fd",f"{par}_nonprompt_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};fd;Entries",len(fd_bin)-1,fd_bin))
                    hpars_total_fd.append(ROOT.TH1F(f"{par}_total_fd",f"{par}_total_fd_{frame}_pt_{pt_min_edge}_{pt_max_edge};fd;Entries",len(fd_bin)-1,fd_bin))

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    hmass_prompt_fd_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_prompt")
                    task_name = f"hprompt_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"

                    fit_set = {"signal_func": ["doublecb"],
                                "bkg_func": ["nobkg"],
                                "chi2_loss": pt_bin_set["chi2_loss"],
                                "mass_range": mass_range,
                                "rebin": pt_bin_set["Rebin"],
                                "bin_counting": pt_bin_set["bin_counting"],
                                "init_pars": [False],
                                "fix_pars":[False],
                                "threshold": [False],
                                "out_dir": os.path.join(os.getcwd(),"Output","Mc_pars")
                            }
                    try:
                        raw_yield, raw_yield, par_dict_prompt_fd = self.fit_inv_mass(outfile_name, hmass_prompt_fd_dir, task_name, fit_set)
                    except:
                        self.logger.error(f"Failed to fit {task_name}!")
                        continue

                    hmass_nonprompt_fd_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"hnonprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_nonprompt")
                    task_name = f"hnonprompt_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"

                    fit_set = {"signal_func": ["doublecb"],
                                "bkg_func": ["nobkg"],
                                "chi2_loss": pt_bin_set["chi2_loss"],
                                "mass_range": mass_range,
                                "rebin": pt_bin_set["Rebin"],
                                "bin_counting": pt_bin_set["bin_counting"],
                                "init_pars": [False],
                                "fix_pars":[False],
                                "threshold": [False],
                                "out_dir": os.path.join(os.getcwd(),"Output","Mc_pars")
                            }
                    try:
                        raw_yield, raw_yield, par_dict_nonprompt_fd = self.fit_inv_mass(outfile_name, hmass_nonprompt_fd_dir, task_name, fit_set)
                    except:
                        self.logger.error(f"Failed to fit {task_name}!")
                        continue

                    hmass_total_fd_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"htotal_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_total")
                    task_name = f"htotal_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}"

                    fit_set = {"signal_func": ["doublecb"],
                                "bkg_func": ["nobkg"],
                                "chi2_loss": pt_bin_set["chi2_loss"],
                                "mass_range": mass_range,
                                "rebin": pt_bin_set["Rebin"],
                                "bin_counting": pt_bin_set["bin_counting"],
                                "init_pars": [False],
                                "fix_pars":[False],
                                "threshold": [False],
                                "out_dir": os.path.join(os.getcwd(),"Output","Mc_pars")
                            }
                    try:
                        raw_yield, raw_yield, par_dict_total_fd = self.fit_inv_mass(outfile_name, hmass_total_fd_dir, task_name, fit_set)
                    except:
                        self.logger.error(f"Failed to fit {task_name}!")
                        continue
                    
                    for par in save_pars:
                        hpars_prompt_fd[save_pars.index(par)].SetBinContent(ifd+1, par_dict_prompt_fd[par])
                        hpars_prompt_fd[save_pars.index(par)].SetBinError(ifd+1, par_dict_prompt_fd[par+"_error"])

                        hpars_nonprompt_fd[save_pars.index(par)].SetBinContent(ifd+1, par_dict_nonprompt_fd[par])
                        hpars_nonprompt_fd[save_pars.index(par)].SetBinError(ifd+1, par_dict_nonprompt_fd[par+"_error"])

                        hpars_total_fd[save_pars.index(par)].SetBinContent(ifd+1, par_dict_total_fd[par])
                        hpars_total_fd[save_pars.index(par)].SetBinError(ifd+1, par_dict_total_fd[par+"_error"])

                    cos_bin = np.array(cos_edges)
                
                    hpars_prompt_cos = []
                    hpars_nonprompt_cos = []
                    hpars_total_cos = []

                    for par in save_pars:
                        hpars_prompt_cos.append(ROOT.TH1F(f"{par}_prompt_cos",f"{par}_prompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin))
                        hpars_nonprompt_cos.append(ROOT.TH1F(f"{par}_nonprompt_cos",f"{par}_nonprompt_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin))
                        hpars_total_cos.append(ROOT.TH1F(f"{par}_total_cos",f"{par}_total_cos_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge};Cos#vartheta*;Entries",len(cos_bin)-1,cos_bin))

                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        cos_dir = fd_dir.mkdir(f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", "", ROOT.kTRUE)
                        cos_dir.cd()

                        hmass_total_cos_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", f"htotal_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}_total")
                        task_name = f"htotal_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}"

                        fit_set = {"signal_func": ["doublecb"],
                                    "bkg_func": ["nobkg"],
                                    "chi2_loss": pt_bin_set["chi2_loss"],
                                    "mass_range": mass_range,
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": pt_bin_set["bin_counting"],
                                    "init_pars": [False],
                                    "fix_pars":[False],
                                    "threshold": [False],
                                    "out_dir": os.path.join(os.getcwd(),"Output","Mc_pars")
                                }
                        try:
                            raw_yield, raw_yield, par_dict_total_cos = self.fit_inv_mass(outfile_name, hmass_total_cos_dir, task_name, fit_set)
                            for par in save_pars:
                                hpars_total_cos[save_pars.index(par)].SetBinContent(icos+1, par_dict_total_cos[par])
                                hpars_total_cos[save_pars.index(par)].SetBinError(icos+1, par_dict_total_cos[par+"_error"])
                        
                        except:
                            self.logger.error(f"Failed to fit {task_name}!")
                            continue
                    

                        hmass_prompt_cos_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", f"hprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}_prompt")
                        task_name = f"hprompt_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}"

                        fit_set = {"signal_func": ["doublecb"],
                                    "bkg_func": ["nobkg"],
                                    "chi2_loss": pt_bin_set["chi2_loss"],
                                    "mass_range": mass_range,
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": pt_bin_set["bin_counting"],
                                    "init_pars": [False],
                                    "fix_pars":[False],
                                    "threshold": [False],
                                    "out_dir": os.path.join(os.getcwd(),"Output","Mc_pars")
                                }
                        try:
                            raw_yield, raw_yield, par_dict_prompt_cos = self.fit_inv_mass(outfile_name, hmass_prompt_cos_dir, task_name, fit_set)
                            for par in save_pars:
                                hpars_prompt_cos[save_pars.index(par)].SetBinContent(icos+1, par_dict_prompt_cos[par])
                                hpars_prompt_cos[save_pars.index(par)].SetBinError(icos+1, par_dict_prompt_cos[par+"_error"])
                        except:
                            self.logger.error(f"Failed to fit {task_name}!")
                            continue
                        
                        hmass_nonprompt_cos_dir = os.path.join(frame, f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}", f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}", f"cos_{cos_min_edge:.1f}_{cos_max_edge:.1f}", f"hnonprompt_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}_nonprompt")
                        task_name = f"hnonprompt_{frame}_pt_{pt_min_edge}_{pt_max_edge}_fd_{fd_min_edge}_{fd_max_edge}_cos_{cos_min_edge}_{cos_max_edge}"

                        fit_set = {"signal_func": ["doublecb"],
                                    "bkg_func": ["nobkg"],
                                    "chi2_loss": pt_bin_set["chi2_loss"],
                                    "mass_range": mass_range,
                                    "rebin": pt_bin_set["Rebin"],
                                    "bin_counting": pt_bin_set["bin_counting"],
                                    "init_pars": [False],
                                    "fix_pars":[False],
                                    "threshold": [False],
                                    "out_dir": os.path.join(os.getcwd(),"Output","Mc_pars")
                                }
                        try:
                            raw_yield, raw_yield, par_dict_nonprompt_cos = self.fit_inv_mass(outfile_name, hmass_nonprompt_cos_dir, task_name, fit_set)
                            for par in save_pars:
                                hpars_nonprompt_cos[save_pars.index(par)].SetBinContent(icos+1, par_dict_nonprompt_cos[par])
                                hpars_nonprompt_cos[save_pars.index(par)].SetBinError(icos+1, par_dict_nonprompt_cos[par+"_error"])
                        except:
                            self.logger.error(f"Failed to fit {task_name}!")
                            continue

                    fd_dir.cd()

                    for par in save_pars:
                        hpars_prompt_cos[save_pars.index(par)].Write("",ROOT.TObject.kOverwrite)
                        hpars_nonprompt_cos[save_pars.index(par)].Write("",ROOT.TObject.kOverwrite)
                        hpars_total_cos[save_pars.index(par)].Write("",ROOT.TObject.kOverwrite)

                pt_bin_dir.cd()
              
                for par in save_pars:
                    hpars_prompt_fd[save_pars.index(par)].Write("",ROOT.TObject.kOverwrite)
                    hpars_nonprompt_fd[save_pars.index(par)].Write("",ROOT.TObject.kOverwrite)
                    hpars_total_fd[save_pars.index(par)].Write("",ROOT.TObject.kOverwrite)

        outfile.Close()