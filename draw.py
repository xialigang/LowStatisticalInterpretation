

from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, TGraph, gPad
import sys
import math
import os
from array import array


gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

dir_results = ''

def get_upperlimit(alpha=0.05, list_mus=[], list_CLs=[]):
    n = len(list_mus)
    if n==0:
        print 'no data?'
        return -1.
    mu = 0.
    for i in range(n-1):
        if list_CLs[i]>=alpha and list_CLs[i+1]<alpha:
            mu0 = list_mus[i]
            mu1 = list_mus[i+1]
            mu = mu0 + (mu1-mu0)*(alpha - list_CLs[i])/(list_CLs[i+1]-list_CLs[i])
            print 'upper limit at ',1-alpha,'is', mu
            break
    return mu
def plot_CLs(dir_results='',nametag='', list_mus=[], alpha=0.05):
    list_CLs = []
    n = len(list_mus)
    a_mu = array('f', [])
    a_CLs = array('f', [])
    for mu in list_mus:
        lines = []
        a_mu.append(float(mu))
        with open (dir_results+'/results_upperlimit_mu'+str(mu)+'.txt','r') as f:
            lines = f.readlines()
        for line in lines:
            if 'CLs = ' not in line:
                continue
            aline = line.strip().split(' ')
            CLs = float(aline[-1])
            a_CLs.append(CLs)
            list_CLs.append(CLs)
            break
    mu_ul = get_upperlimit(alpha, list_mus, list_CLs)
    g = TGraph(n, a_mu, a_CLs)
    Cs_g = TCanvas('Cs_g', '', 10, 10, 800, 600)
    g.Draw('APL')
    g.SetLineWidth(2)
    g.GetXaxis().SetTitle("#mu")
    g.GetYaxis().SetTitle("CL_{s}")
    gPad.Update()
    xmin = gPad.GetUxmin()
    xmax = gPad.GetUxmax()
    line_alpha = TLine(xmin, 0.05, xmax, 0.05)
    line_alpha.SetLineStyle(2)
    line_alpha.SetLineWidth(2)
    line_alpha.SetLineColor(4)
    line_alpha.Draw()
    leg = TLegend(0.5, 0.6, 0.95, 0.9)
    leg.AddEntry(line_alpha, 'U.L.(%i%%) = %.1f' % (100-alpha*100, mu_ul) ,'L')
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.Draw()
    plotname = 'Cs_CLs'
    if nametag != '':
        plotname += '_' + nametag
    Cs_g.SaveAs(plotname + '.png')
    Cs_g.SaveAs(plotname + '.pdf')
    return


list_mus = [10, 30, 50, 75, 100]
plot_CLs('pic_obs_lephad1', 'obs_lephad1', list_mus)

list_mus = [1, 10, 30, 50, 75]
plot_CLs('pic_obs_0', 'obs_0', list_mus)



