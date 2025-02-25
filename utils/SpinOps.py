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

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "AnalysisSpinAlignment.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
            type_dir.cd()

            for pt_min_edge, pt_max_edge in zip(pt_edges[:-1], pt_edges[1:]):

                pt_bin_set = self.config.BinSet["pt_bin_set"][f"{pt_min_edge:.0f}-{pt_max_edge:.0f}"]

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

                fd_bins = np.array(fd_edges)
                hraw_rho = ROOT.TH1F("hraw_rho", "hraw_rho;fd_score;#rho", len(fd_bins)-1, fd_bins)
                hcorr_rho = ROOT.TH1F("hcorr_rho", "hcorr_rho;fd_score;#rho", len(fd_bins)-1, fd_bins)

                for ifd, (fd_min_edge, fd_max_edge) in enumerate(zip(fd_min_edges, fd_max_edges)):

                    fd_dir = pt_bin_dir.mkdir(f"fd_{fd_min_edge:.2f}_{fd_max_edge:.2f}","",ROOT.kTRUE)
                    fd_dir.cd()

                    raw_yield_hist = fd_dir.Get("hraw_yield")
                    corr_yield_hist = fd_dir.Get("hcorr_yield")
                    if not raw_yield_hist or not corr_yield_hist:
                        print(f"Raw yield or efficiency histogram not found in {fd_dir.GetName()}")
                        continue

                    raw_fit_plot, raw_fit_rho = self.write_rho_plots(raw_yield_hist)
                    raw_fit_plot.Write("raw_yield_fit", ROOT.TObject.kOverwrite)
                    hraw_rho.SetBinContent(ifd+1, raw_fit_rho[0])
                    hraw_rho.SetBinError(ifd+1, raw_fit_rho[1])

                    corr_fit_plot, corr_fit_pars = self.write_rho_plots(corr_yield_hist)
                    corr_fit_plot.Write("corr_yield_fit", ROOT.TObject.kOverwrite)
                    hcorr_rho.SetBinContent(ifd+1, corr_fit_pars[0])
                    hcorr_rho.SetBinError(ifd+1, corr_fit_pars[1])

                pt_bin_dir.cd()
                hraw_rho.Write("", ROOT.TObject.kOverwrite)
                hcorr_rho.Write("", ROOT.TObject.kOverwrite)

                if frame == "Helicity":
                    frame_axis = 1
                elif frame == "Production":
                    frame_axis = 2
                else:
                    continue
                simu_dir = pt_bin_dir.mkdir("Simulation","",ROOT.kTRUE)

                simu_file_list = [os.path.join(self.config.Directories["InputDir"], "MC", file_name) for file_name in self.config.Target["simulation"]]
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

                pro_plot, pro_rho = self.write_rho_plots(prompt_hist)
                simu_dir.cd()
                pro_plot.Write("Prompt_yield_fit", ROOT.TObject.kOverwrite)
                nonpro_plot, nonpro_rho = self.write_rho_plots(nonprompt_hist)
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

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "AnalysisSpinAlignment.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
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
                else:
                    continue
                simu_dir = pt_bin_dir.mkdir("Simulation","",ROOT.kTRUE)

                simu_file_list = [os.path.join(self.config.Directories["InputDir"], "MC", file_name) for file_name in self.config.Target["simulation"]]
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

                pro_plot, pro_rho = self.write_rho_plots(prompt_hist)
                simu_dir.cd()
                pro_plot.Write("Prompt_yield_fit", ROOT.TObject.kOverwrite)
                nonpro_plot, nonpro_rho = self.write_rho_plots(nonprompt_hist)
                simu_dir.cd()
                nonpro_plot.Write("NonPrompt_yield_fit", ROOT.TObject.kOverwrite)

                rho_list = np.array([pro_rho[0],nonpro_rho[0]],dtype=np.float64)
                rho_error_list = np.array([pro_rho[1],nonpro_rho[1]],dtype=np.float64)
                simu_rho = ROOT.TGraphErrors(2,np.array([0,1],dtype=np.float64), rho_list, np.array([0,0],dtype=np.float64), rho_error_list)
                
                simu_dir.cd()
                simu_rho.Write("simu_rho", ROOT.TObject.kOverwrite) 
  
        outfile.Close()

    def write_rho_plots(self, hraw_yield):
                    
        func = ROOT.TF1("frho00","[0]*((1-[1])+(3*[1]-1)*x*x)",-1,1)
        func.SetParameter(0,1)
        func.SetParameter(1,0.3)
        func.SetParLimits(1,0,1)

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

        outfile_name = os.path.join(self.out_dir, "Analysis-root", "AnalysisSpinAlignment.root")
        outfile = ROOT.TFile(outfile_name, "UPDATE")
        pt_edges = self.config.BinSet["pt_bin_edges"]  
        frame_list = self.config.Analysis["Framework"]

        for frame in frame_list:

            type_dir = outfile.mkdir(frame,"",ROOT.kTRUE)
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

                ConInt = ROOT.TH1D("ConInt","ConInt",1000,-0.1,1)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(ConInt, 0.68)
            
                TGraph.GetXaxis().SetRangeUser(0,1)
                ConInt.SetLineColor(ROOT.kAzure+4)
                ConInt.SetFillColorAlpha(ROOT.kAzure+4,0.3)
                ConInt.SetFillStyle(1001)
                ConInt.SetMarkerStyle(0)
                ConInt.DrawCopy("SAME")
                ConInt.Write(f"ConInt_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

                ConInt2 = ROOT.TH1D("ConInt2","ConInt2",1000,-0.1,1)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(ConInt2, 0.95)

                ConInt2.SetLineColor(ROOT.kAzure+4)
                ConInt2.SetFillColorAlpha(ROOT.kAzure+4,0.3)
                ConInt2.SetFillStyle(1001)
                ConInt2.SetMarkerStyle(0)
                ConInt2.DrawCopy("SAME")
                ConInt2.Write(f"ConInt2_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

                ConInt3 = ROOT.TH1D("ConInt3","ConInt3",1000,-0.1,1)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(ConInt3, 0.99)

                ConInt3.SetLineColor(ROOT.kAzure+4)
                ConInt3.SetFillColorAlpha(ROOT.kAzure+4,0.3)
                ConInt3.SetFillStyle(1001)
                ConInt3.SetMarkerStyle(0)
                ConInt3.DrawCopy("SAME")
                ConInt3.Write(f"ConInt3_{pt_min_edge}_{pt_max_edge}",ROOT.TObject.kOverwrite)

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
        