
from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, gStyle
import sys
import math

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)



s = [0.0092, 0.025, 0.023, 0.03]
b = [0.12, 0.0048, 0.0043, 0.028]
syst_1_s = [ 0, 0, 0, 0]
syst_1_b = [ 0.002, 0, 0, 0.0005]
data = [0, 0, 0, 0]
def prepare_input():
    nbins = len(s)
    file0 = TFile('input.root', 'recreate')
    h_s = TH1F('h_sig', '', nbins, 0, nbins)
    h_b = TH1F('h_bkg', '', nbins, 0, nbins)
    h_data = TH1F('h_data', '', nbins, 0, nbins)
    h_syst_1_s = TH1F('h_syst_1_sig', '', nbins, 0, nbins)
    h_syst_1_b = TH1F('h_syst_1_bkg', '', nbins, 0, nbins)
    for i in range(nbins):
        h_s.SetBinContent(i+1, s[i])
        h_b.SetBinContent(i+1, b[i])
        h_data.SetBinContent(i+1, data[i])
        h_syst_1_s.SetBinContent(i+1, syst_1_s[i])
        h_syst_1_b.SetBinContent(i+1, syst_1_b[i])
    h_s.Write()
    h_b.Write()
    h_data.Write()
    h_syst_1_s.Write()
    h_syst_1_b.Write()
    file0.Close()
    return

prepare_input()


