import os
import sys
import numpy as np

import ROOT
ROOT.gROOT.SetBatch(True)

class SpinOps(object):

    def __init__ (self,config, logger):
        self.config = config
        self.out_dir = os.path.join(self.config.Directories["OutputDir"], self.config.Analysis["Task_Name"])
        self.logger = logger

    def get_rho(self):

        file_name = os.path.join(self.out_dir, "Analysis-root", "Frac_And_Rho.root")
        outfile = ROOT.TFile(file_name, "UPDATE")
        ana_dir = outfile.mkdir(self.config.Analysis["Ana_name"],"",ROOT.kTRUE)
        ana_dir.cd()

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            type_dir = ana_dir.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]

                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]
                
                if not pt_bin_set["doing"]:
                    continue
                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)

                fd_bins = np.array(fd_edges)
                hraw_rho = ROOT.TH1F("hraw_rho", "hraw_rho;fd_score;#rho", len(fd_bins)-1, fd_bins)
                hcorr_rho = ROOT.TH1F("hcorr_rho", "hcorr_rho;fd_score;#rho", len(fd_bins)-1, fd_bins)

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    if len(fd_edges) > 2:
                        bin_num = ifd
                    else:
                        bin_num = ifd+1

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    raw_yield_hist = fd_dir.Get("hraw_yield")
                    corr_yield_hist = fd_dir.Get("hcorr_yield")
                    if not raw_yield_hist or not corr_yield_hist:
                        print(f"Raw yield or efficiency histogram not found in {fd_dir.GetName()}")
                        continue

                    raw_fit_plot, raw_fit_rho = self.get_rho_plot(raw_yield_hist)
                    raw_fit_plot.Write("raw_yield_fit", ROOT.TObject.kOverwrite)
                    hraw_rho.SetBinContent(bin_num, raw_fit_rho[0])
                    hraw_rho.SetBinError(bin_num, raw_fit_rho[1])

                    corr_fit_plot, corr_fit_pars = self.get_rho_plot(corr_yield_hist)
                    corr_fit_plot.Write("corr_yield_fit", ROOT.TObject.kOverwrite)
                    hcorr_rho.SetBinContent(bin_num, corr_fit_pars[0])
                    hcorr_rho.SetBinError(bin_num, corr_fit_pars[1])

                pt_bin_dir.cd()
                hraw_rho.Write("", ROOT.TObject.kOverwrite)
                hcorr_rho.Write("", ROOT.TObject.kOverwrite)

                if frame == "Helicity":
                    frame_axis = 1
                elif frame == "Production":
                    frame_axis = 2
                elif frame == "EP":
                    frame_axis = 2 # Temporary use production frame for EP frame
                else:
                    continue
                simu_dir = pt_bin_dir.mkdir("Simulation","",ROOT.kTRUE)

                simu_file_list = [os.path.join(self.config.Directories["InputDir"], "MC", file_name) for file_name in self.config.Files["simulation"]]
                TPrompt = self.get_sparse(simu_file_list, "TPrompt")
                TNonPrompt = self.get_sparse(simu_file_list, "TNonPrompt")

                pt_min = TPrompt.GetAxis(0).FindBin(pt_min_edge+1e-6)
                pt_max = TPrompt.GetAxis(0).FindBin(pt_max_edge-1e-6)
                TPrompt.GetAxis(0).SetRange(pt_min, pt_max)
                TNonPrompt.GetAxis(0).SetRange(pt_min, pt_max)

                prompt_hist = TPrompt.Projection(frame_axis).Clone("prompt_hist")
                prompt_hist.Write("Prompt_yield", ROOT.TObject.kOverwrite)
                nonprompt_hist = TNonPrompt.Projection(frame_axis).Clone("nonprompt_hist")
                nonprompt_hist.Write("NonPrompt_yield", ROOT.TObject.kOverwrite)

                pro_plot, pro_rho = self.get_rho_plot(prompt_hist)
                simu_dir.cd()
                pro_plot.Write("Prompt_yield_fit", ROOT.TObject.kOverwrite)
                nonpro_plot, nonpro_rho = self.get_rho_plot(nonprompt_hist)
                simu_dir.cd()
                nonpro_plot.Write("NonPrompt_yield_fit", ROOT.TObject.kOverwrite)


                # 将pro_rho和nonpro_rho作为在0和1处的两个点写入文件
                rho_list = np.array([pro_rho[0],nonpro_rho[0]],dtype=np.float64)
                rho_error_list = np.array([pro_rho[1],nonpro_rho[1]],dtype=np.float64)
                simu_rho = ROOT.TGraphErrors(2,np.array([0,1],dtype=np.float64), rho_list, np.array([0,0],dtype=np.float64), rho_error_list)
                
                simu_dir.cd()
                simu_rho.Write("simu_rho", ROOT.TObject.kOverwrite) 

        outfile.Close()

    def write_simu_rho(self):

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "Frac_And_Rho.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        ana_dir = outfile.mkdir(self.config.Analysis["Ana_name"],"",ROOT.kTRUE)
        ana_dir.cd()

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            type_dir = ana_dir.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                
                if not pt_bin_set["doing"]:
                    continue
                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)

                if frame == "Helicity":
                    frame_axis = 1
                elif frame == "Production":
                    frame_axis = 2
                elif frame == "EP":
                    frame_axis = 2 #It should be random axis but we don't have the info in the simulation
                else:
                    continue
                simu_dir = pt_bin_dir.mkdir("Simulation","",ROOT.kTRUE)

                simu_file_list = [os.path.join(self.config.Directories["InputDir"], "MC", file_name) for file_name in self.config.Files["simulation"]]
                TPrompt = self.get_sparse(simu_file_list, "TPrompt")
                TNonPrompt = self.get_sparse(simu_file_list, "TNonPrompt")

                pt_min = TPrompt.GetAxis(0).FindBin(pt_min_edge+1e-6)
                pt_max = TPrompt.GetAxis(0).FindBin(pt_max_edge-1e-6)
                TPrompt.GetAxis(0).SetRange(pt_min, pt_max)
                TNonPrompt.GetAxis(0).SetRange(pt_min, pt_max)

                prompt_hist = TPrompt.Projection(frame_axis).Clone("prompt_hist")
                prompt_hist.Write("Prompt_yield", ROOT.TObject.kOverwrite)
                nonprompt_hist = TNonPrompt.Projection(frame_axis).Clone("nonprompt_hist")
                nonprompt_hist.Write("NonPrompt_yield", ROOT.TObject.kOverwrite)

                pro_plot, pro_rho = self.get_rho_plot(prompt_hist)
                simu_dir.cd()
                pro_plot.Write("Prompt_yield_fit", ROOT.TObject.kOverwrite)
                nonpro_plot, nonpro_rho = self.get_rho_plot(nonprompt_hist)
                simu_dir.cd()
                nonpro_plot.Write("NonPrompt_yield_fit", ROOT.TObject.kOverwrite)

                rho_list = np.array([pro_rho[0],nonpro_rho[0]],dtype=np.float64)
                rho_error_list = np.array([pro_rho[1],nonpro_rho[1]],dtype=np.float64)
                simu_rho = ROOT.TGraphErrors(2,np.array([0,1],dtype=np.float64), rho_list, np.array([0,0],dtype=np.float64), rho_error_list)
                
                simu_dir.cd()
                simu_rho.Write("simu_rho", ROOT.TObject.kOverwrite) 
  
        outfile.Close()

    def get_rho_plot(self, hraw_yield):
                    
        func = ROOT.TF1("frho00","[0]*((1-[1])+(3*[1]-1)*x*x)",-1,1)
        func.SetParameter(0,1)
        func.SetParameter(1,0.3)
        func.SetParLimits(1,0,0.8)

        hraw_yield.Fit(func,"R,S")

        chi2 = func.GetChisquare()
        ndf = func.GetNDF()
        chi2ndf = chi2/ndf
        canvas = ROOT.TCanvas("rho_fit","rho_fit",800,600)
        canvas.cd()

        hraw_yield.SetLineColor(ROOT.kBlack)
        hraw_yield.SetMarkerColor(ROOT.kBlack)
        hraw_yield.SetMarkerStyle(20)
        hraw_yield.SetMarkerSize(0.5)
        hraw_yield.SetStats(0)
        hraw_yield.GetYaxis().SetRangeUser(0,(hraw_yield.GetMaximum()+hraw_yield.GetBinError(1))*1.1)
        hraw_yield.Draw("PE")
        func.SetLineColor(ROOT.kRed)
        func.SetLineWidth(2)
        func.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetTextSize(0.03)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.6,0.86,f"#chi^{{2}}/ndf = {chi2ndf:.3f}")
        latex.DrawLatexNDC(0.6,0.8,f"N_{{0}} = {func.GetParameter(0):.3f} #pm {func.GetParError(0):.3f}")
        latex.DrawLatexNDC(0.6,0.75,f"#rho_{{00}} = {func.GetParameter(1):.3f} #pm {func.GetParError(1):.3f}")
        latex.Draw("same")

        rho = [func.GetParameter(1),func.GetParError(1)]
        return canvas, rho

    def extro_rho(self):

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "Frac_And_Rho.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        ana_dir = outfile.mkdir(self.config.Analysis["Ana_name"],"",ROOT.kTRUE)
        ana_dir.cd()

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            type_dir = ana_dir.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                
                if not pt_bin_set["doing"]:
                    continue
                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                pt_bin_dir.cd()

                fd_edges = pt_bin_set["fd_edges"]

                hraw_frc = pt_bin_dir.Get("hfrac")
                hcorr_rho = pt_bin_dir.Get("hcorr_rho")
                
                frc_list = []
                frc_error = []
                rho_list = []
                rho_error = []

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_edges[:-1], fd_edges[1:])):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    frc_list.append(hraw_frc.GetBinContent(ifd+1))
                    frc_error.append(hraw_frc.GetBinError(ifd+1))
                    rho_list.append(hcorr_rho.GetBinContent(ifd+1)) 
                    rho_error.append(hcorr_rho.GetBinError(ifd+1))
                
                if len(fd_edges) < 3:
                    frc_list.append(1)
                    frc_error.append(0)
                    hHardrho = pt_bin_dir.Get("Simulation/simu_rho")
                    rho_list.append(hHardrho.GetY()[1])
                    rho_error.append(hHardrho.GetErrorY(1))
        
                frc_list = np.array(frc_list,dtype=np.float64)
                frc_error = np.array(frc_error,dtype=np.float64)
                rho_list = np.array(rho_list,dtype=np.float64)
                rho_error = np.array(rho_error,dtype=np.float64)
                
                pt_bin_dir.cd()
                canvas = ROOT.TCanvas(f"prompt_fraction_{pt_min_edge}_{pt_max_edge}",f"prompt_fraction_{pt_min_edge}_{pt_max_edge}",800,600)
                canvas.cd().DrawFrame(
                    0., 0.2, 1., 0.6,
                    ';#it{f}_{non-prompt};#it{#rho}_{00}'
                )

                TGraph = ROOT.TGraphErrors(len(frc_list),frc_list,rho_list,frc_error,rho_error)
                TGraph.SetName(f"rho_prompt_{pt_min_edge}_{pt_max_edge}")
                TGraph.SetTitle(f"Non-prompt fraction vs rho_{pt_min_edge}_{pt_max_edge}")
                TGraph.GetXaxis().SetTitle("f_{non-prompt}")
                TGraph.GetYaxis().SetTitle("#rho_{00}")
                TGraph.SetMarkerStyle(20)
                TGraph.SetMarkerSize(0.5)
                TGraph.Write(f"rhovsfrc_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

                linearFit = ROOT.TF1("linearFit", "[0] + [1]*x", -0.1, 1)
                TGraph.Fit(linearFit)

                ConInt3 = ROOT.TH1D("ConInt3","ConInt3",1000,-0.1,1)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(ConInt3, 0.99)

                ConInt3.SetLineColor(ROOT.kAzure)
                ConInt3.SetFillColorAlpha(ROOT.kAzure,0.3)
                ConInt3.SetFillStyle(1001)
                ConInt3.SetMarkerStyle(0)
                ConInt3.DrawCopy("SAME")
                ConInt3.Write(f"ConInt3_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

                ConInt2 = ROOT.TH1D("ConInt2","ConInt2",1000,-0.1,1)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(ConInt2, 0.95)

                ConInt2.SetLineColor(ROOT.kAzure+4)
                ConInt2.SetFillColorAlpha(ROOT.kAzure+4,0.3)
                ConInt2.SetFillStyle(1001)
                ConInt2.SetMarkerStyle(0)
                ConInt2.DrawCopy("SAME")
                ConInt2.Write(f"ConInt2_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

                ConInt = ROOT.TH1D("ConInt","ConInt",1000,-0.1,1)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(ConInt, 0.68)
            
                TGraph.GetXaxis().SetRangeUser(0,1)
                ConInt.SetLineColor(ROOT.kAzure+8)
                ConInt.SetFillColorAlpha(ROOT.kAzure+8,0.3)
                ConInt.SetFillStyle(1001)
                ConInt.SetMarkerStyle(0)
                ConInt.DrawCopy("SAME")
                ConInt.Write(f"ConInt_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

                TGraph.Draw("EPSAME")
                TGraph.Write(f"fit_rhovsfrac",ROOT.TObject.kOverwrite)
            
                linearFit.Draw("LSAME")
                linearFit.Write(f"linearFit_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)
                
                pt_bin_dir.cd()
                canvas.Write("",ROOT.TObject.kOverwrite)

                hrho = ROOT.TH1F("hrho","#rho_{00};f_nonprompt;#rho_{00}",2,-0.5,1.5)
                hrho.SetBinContent(1,linearFit.GetParameter(0))
                hrho.SetBinError(1,linearFit.GetParError(0))
                hrho.SetBinContent(2,linearFit.GetParameter(0)+linearFit.GetParameter(1))
                hrho.SetBinError(2,linearFit.GetParError(0)+linearFit.GetParError(1))
                hrho.Write(f"rho_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

        outfile.Close()
        
    def plot_rho(self):

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "Frac_And_Rho.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        ana_dir = outfile.Get(self.config.Analysis["Ana_name"])
        ana_dir.cd()

        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            self.logger.info(f"Plotting rho vs pt in {frame} frame")
            type_dir = ana_dir.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            hprompt_rho = ROOT.TH1F("hprompt_rho","hprompt_rho;fd_score;#rho", len(pt_edges)-1, np.array(pt_edges,dtype=np.float64))
            hnonprompt_rho = ROOT.TH1F("hnonprompt_rho","hnonprompt_rho;fd_score;#rho", len(pt_edges)-1, np.array(pt_edges,dtype=np.float64))

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                
                if not pt_bin_set["doing"]:
                    continue
                pt_bin_dir = type_dir.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)

                pt_rho = pt_bin_dir.Get(f"rho_{pt_min_edge:.0f}_{pt_max_edge:.0f}")
                hprompt_rho.SetBinContent(hprompt_rho.FindBin(pt_min_edge+1e-6), pt_rho.GetBinContent(1))
                hprompt_rho.SetBinError(hprompt_rho.FindBin(pt_min_edge+1e-6), pt_rho.GetBinError(1))
                hnonprompt_rho.SetBinContent(hnonprompt_rho.FindBin(pt_min_edge+1e-6), pt_rho.GetBinContent(2))
                hnonprompt_rho.SetBinError(hnonprompt_rho.FindBin(pt_min_edge+1e-6), pt_rho.GetBinError(2))

            canves = ROOT.TCanvas(f"rho_{frame}",f"rho_{frame}",800,600)

            hnonprompt_rho.SetMarkerStyle(20)
            hnonprompt_rho.SetMarkerSize(0.5)
            hnonprompt_rho.SetMarkerColor(ROOT.kRed+1)
            hnonprompt_rho.SetLineWidth(2)
            hnonprompt_rho.SetLineColor(ROOT.kRed+1)
            hnonprompt_rho.GetYaxis().SetRangeUser(0.1,0.7)
            hnonprompt_rho.Draw("P")

            hprompt_rho.SetMarkerStyle(20)
            hprompt_rho.SetMarkerSize(0.5)
            hprompt_rho.SetMarkerColor(ROOT.kBlue+1)
            hprompt_rho.SetLineWidth(2)
            hprompt_rho.SetLineColor(ROOT.kBlue+1)
            hprompt_rho.Draw("PSAME")

            line = ROOT.TLine(hprompt_rho.GetXaxis().GetXmin(),1/3,hprompt_rho.GetXaxis().GetXmax(),1/3)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("SAME")

            leg = ROOT.TLegend(0.6, 0.2, 0.85, 0.3)
            leg.AddEntry(hprompt_rho, "prompt D*^{+}", "p")
            leg.AddEntry(hnonprompt_rho, "non-prompt D*^{+}", "p")
            leg.SetBorderSize(0)
            leg.Draw("same")

            hprompt_rho.Write(f"hprompt_rho_{frame}",ROOT.TObject.kOverwrite)
            hnonprompt_rho.Write(f"hnonprompt_rho_{frame}",ROOT.TObject.kOverwrite)
            canves.Write(f"rho{frame}",ROOT.TObject.kOverwrite)
 
    def get_sparse(self, file_list, name):
        
        sparse = []
        for file_dir in file_list:
            file = ROOT.TFile(file_dir, "READ")
            sparse.append(file.Get(name))

        for i in range(1, len(sparse)):
            try:
                sparse[0].Add(sparse[i])
                print(sparse[i].GetName())
                sparse[i].GetAxis(0).SetRange(0,0)
                sparse[i].Delete()
            except:
                print(f"Error when mergeing {name}")
                sys.exit(1)
            
        return sparse[0]
        
    def read_frac(self):
        
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        infile_eff_name = os.path.join(self.out_dir, "Analysis-root", "Data_And_Efficiency.root")
        infile_yield_name = os.path.join(self.out_dir, "Analysis-root", "RawYield_Extraction.root")
        outfile_name = os.path.join(self.out_dir, "Analysis-root", "Frac_And_Rho.root")
        # self.config.Analysis["Ana_name"]+".root")

        infile_eff = ROOT.TFile(infile_eff_name, "READ")
        ana_dir_eff = infile_eff.Get(self.config.Analysis["Ana_name"])
        infile_yield = ROOT.TFile(infile_yield_name, "READ")
        ana_dir_yield = infile_yield.Get(self.config.Analysis["Ana_name"])
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        ana_dir_out = outfile.mkdir(self.config.Analysis["Ana_name"],"",ROOT.kTRUE)
        ana_dir_out.cd()

        frame_list = self.config.Analysis["Framework"]
        for frame in frame_list:

            self.logger.info(f"Reading fraction of {frame} framework...")
            type_dir_eff = ana_dir_eff.mkdir(frame,"",ROOT.kTRUE)
            type_dir_yield = ana_dir_yield.Get(frame)
            type_dir_out = ana_dir_out.mkdir(frame,"",ROOT.kTRUE)
            type_dir_out.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]
                if not pt_bin_set["doing"]:
                    continue
                self.logger.info(f"    Working in pt bin {pt_min_edge:.0f}-{pt_max_edge:.0f}...")

                pt_bin_dir_eff = type_dir_eff.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                pt_bin_dir_yield = type_dir_yield.Get(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}")
                pt_bin_dir_out = type_dir_out.mkdir(f"pt_{pt_min_edge:.0f}_{pt_max_edge:.0f}","",ROOT.kTRUE)
                pt_bin_dir_out.cd()

                cos_edges = pt_bin_set["cos_bin_edges"]
                fd_edges = pt_bin_set["fd_edges"]
                if len(fd_edges) > 2:
                    fd_min_edges = [0.0]+fd_edges[:-1]
                    fd_max_edges = [1.0]+fd_edges[1:]
                else:
                    fd_min_edges = fd_edges[:-1]
                    fd_max_edges = fd_edges[1:]

                heff_prompt_fd = pt_bin_dir_eff.Get("hEff_prompt_fd")
                heff_nonprompt_fd = pt_bin_dir_eff.Get(f"hEff_nonprompt_fd")

                hfrac_np_fd = heff_prompt_fd.Clone("hfrac_np_fd")

                frac_file_dir = os.path.join(self.out_dir,f"Cut-variation/{self.config.Analysis['Ana_name']}/pt_{pt_min_edge:0d}_{pt_max_edge:0d}/fraction","CutVar_"+self.config.Analysis['Ana_name']+".root")
                frac_file = ROOT.TFile(frac_file_dir,"READ")
                corr_yield_prompt = frac_file.Get("hCorrYieldsPrompt")
                corr_yield_nonprompt = frac_file.Get("hCorrYieldsNonPrompt")
                Cov_prompt_nonprompt = frac_file.Get("hCovPromptNonPrompt")

                hEff_total = pt_bin_dir_eff.Get("fd_0.00_1.00/hEff_total_cos")
                
                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    if len(fd_edges) > 2:
                        bin_num = ifd
                    else:
                        bin_num = ifd+1

                    fd_dir_eff = pt_bin_dir_eff.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir_yield = pt_bin_dir_yield.Get(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}")
                    fd_dir_out = pt_bin_dir_out.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir_out.cd()

                    raw_np_frc, raw_np_frc_error = self.get_nonprompt_frac(heff_prompt_fd.GetBinContent(bin_num),
                                                                heff_nonprompt_fd.GetBinContent(bin_num),
                                                                corr_yield_prompt.GetBinContent(1),
                                                                corr_yield_nonprompt.GetBinContent(1),
                                                                corr_yield_prompt.GetBinError(1),
                                                                corr_yield_nonprompt.GetBinError(1),
                                                                Cov_prompt_nonprompt.GetBinContent(1))
 
                    hfrac_np_fd.SetBinContent(bin_num, raw_np_frc)
                    hfrac_np_fd.SetBinError(bin_num, raw_np_frc_error)

                    hraw_yield = fd_dir_yield.Get("hraw_yield")
                    hEff_cos = fd_dir_eff.Get("hEff_total_cos")
                    # hp_eff_cos = fd_dir.Get("hEff_prompt_cos")
                    # hnp_eff_cos = fd_dir.Get("hEff_nonprompt_cos")
                    hfrac_cos = hraw_yield.Clone("hfrac_cos")

                    for icos,(cos_min_edge,cos_max_edge) in enumerate(zip(cos_edges[:-1], cos_edges[1:])):

                        hfrac_cos.SetBinContent(icos+1, raw_np_frc)
                        hfrac_cos.SetBinError(icos+1,raw_np_frc_error)
  
                    # Corrected yield with efficiency in different fd bin
                    # hcorr_yield_cos = hraw_yield.Clone("hcorr_yield_cos")
                    # hcorr_yield_cos.Divide(hEff_cos)
                    # hcorr_yield_cos.Write("hcorr_yield_fd",ROOT.TObject.kOverwrite)

                    # hEff_frac , hcorr_yield_frac = self.get_total_eff(hraw_yield,hfrac_cos,hp_eff_cos,hnp_eff_cos)
                    hcorr_yield_total = hraw_yield.Clone("hcorr_yield_total")
                    hcorr_yield_total.Divide(hEff_total)
                    hcorr_yield_total.GetYaxis().SetTitle("Corrected Yield")
                    hcorr_yield_total.SetTitle("Corrected Yield")

                    hEff_cos.Write("hEff_fdbin",ROOT.TObject.kOverwrite)
                    hEff_total.Write("hEff_total",ROOT.TObject.kOverwrite)
                    hraw_yield.Write("hraw_yield",ROOT.TObject.kOverwrite)
                    hcorr_yield_total.Write("hcorr_yield",ROOT.TObject.kOverwrite)

                pt_bin_dir_out.cd()
                hfrac_np_fd.Write("hfrac",ROOT.TObject.kOverwrite)

        infile_eff.Close()
        infile_yield.Close()
        outfile.Close()

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