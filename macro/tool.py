import os
import ROOT
import numpy as np

ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gROOT.SetBatch(True)

def plot_ratio(hist1, hist2, set1, set2, pt, title,outfile,scale = 0):

    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    outfile.cd()

    canv = ROOT.TCanvas(f"{title}_{pt}",f"{title}_{pt}",1000,600)
    canv.Divide(2,1)
    canv.cd(1)

    hist1.SetLineColor(ROOT.kRed)
    hist2.SetLineColor(ROOT.kBlue)

    hist1.SetMarkerColor(ROOT.kRed)
    hist2.SetMarkerColor(ROOT.kBlue)

    hist1.SetMarkerStyle(20)
    hist2.SetMarkerStyle(20)

    hist1.SetMarkerSize(0.5)
    hist2.SetMarkerSize(0.5)

    if scale == 1:
        print(f"Scaling {title} to 1")
        hist1.Scale(1/hist1.Integral())
        hist2.Scale(1/hist2.Integral())

    hist1.SetTitle(f"{title} vs Cut Number for Pt {pt}")

    if title == "prompt_pt_ratio" or title == "nonprompt_pt_ratio":
        hist1.GetXaxis().SetTitle("Pt")
    else: 
        hist1.GetXaxis().SetTitle("Cut Number")

    hist1.GetYaxis().SetTitle(f"{title}")
    hist1.SetStats(0)

    hist1.GetYaxis().SetRangeUser(0.8*np.min([hist1.GetMinimum(),hist2.GetMinimum()]),1.1*np.max([hist1.GetMaximum(),hist2.GetMaximum()]))

    hist1.Draw("E1")
    hist2.Draw("E1 same")

    leg_1 = ROOT.TLegend(0.4,0.7,0.8,0.9)
    leg_1.AddEntry(hist1,f"{set1}","l")
    leg_1.AddEntry(hist2,f"{set2}","l")
    leg_1.Draw("same")

    canv.cd(2)
    canv.Update()

    ratio_hist = hist1.Clone(f"ratio_hist")
    ratio_hist.Divide(ratio_hist,hist2,1,1,"B")

    ratio_hist.SetLineColor(ROOT.kBlack)
    ratio_hist.SetMarkerColor(ROOT.kBlack)
    ratio_hist.SetMarkerStyle(20)
    ratio_hist.SetMarkerSize(0.5)

    ratio_hist.SetTitle(f"{set1}/{set2} Ratio vs Cut Number for Pt {pt}")
    if title == "prompt_pt_ratio" or title == "nonprompt_pt_ratio":
        ratio_hist.GetXaxis().SetTitle("Pt")
    else:
        ratio_hist.GetXaxis().SetTitle("Cut Number")
    ratio_hist.GetYaxis().SetTitle(f"{set1}/{set2}")

    ratio_hist.GetYaxis().SetRangeUser(0.8*ratio_hist.GetMinimum(),1.1*ratio_hist.GetMaximum())
    ratio_hist.Draw("E1")

    leg_2 = ROOT.TLegend(0.4,0.7,0.8,0.9)
    leg_2.AddEntry(ratio_hist,f"{set1}/{set2}","l")
    leg_2.Draw("same")

    canv.Update()
    canv.Write()

