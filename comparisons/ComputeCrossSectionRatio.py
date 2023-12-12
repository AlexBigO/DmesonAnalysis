"""
Script to compute cross section ratio between two particle species
"""

import sys
import numpy as np
from ROOT import TGraphErrors # needed for normalisation syst unc.
from ROOT import TCanvas, TFile, TLegend, TLine, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kOpenSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins

inputdir = '/home/abigot/AnalysisNonPromptDplus/'
prompt_ratio_filename = 'Ratio_Dplus_D0_prompt_pPb.root'
inputfilenames = ['Run2pPb5Tev/4_Analysis/5_CrossSection/wptweights_cent/cross_section.root',
                  'CrossSection_D0_pPb5TeV_FD.root']  # CrossSectionDplus_pPb_nonprompt.root']
histonames = ['hCorrYield', 'hCrossSection']
graphnames = ['gCorrYieldSystTot', 'gCrossSectionSystTot']
colors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
linecolors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
markers = [kFullDiamond, kFullCircle]
markers_ratio = [kOpenDiamond, kOpenCircle]
legendnames = ['Non-prompt D^{+}', 'Non-prompt D^{0}']
legend_subcaption = ['', '']
outputsuffix = 'Dplus_D0'

prompt_ratio_file = TFile('%s/%s' % (inputdir, prompt_ratio_filename))
hRatioPrompt = prompt_ratio_file.Get('hRatio')
gRatioSystTotPrompt = prompt_ratio_file.Get('gRatioSystTot') # no BR


hCorrYield, gCorrYield, hCorrYieldRatio = ([] for _ in range(3))

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6,
               titlesize=0.045, labelsize=0.04)

leg = TLegend(0.42, 0.5, 0.82, 0.75)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

legRatio = TLegend(0.7, 0.2, 0.9, 0.4)
legRatio.SetFillStyle(0)
legRatio.SetBorderSize(0)
legRatio.SetTextSize(0.04)

for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hCorrYield.append(inputfile.Get(histonames[iFile]))

    gCorrYield.append(inputfile.Get(graphnames[iFile]))
    hCorrYield[iFile].SetDirectory(0)
    SetObjectStyle(hCorrYield[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile])  #, markersize=1.5)
    SetObjectStyle(gCorrYield[iFile], linecolor=colors[iFile], fillstyle=0)

# legend
for ileg, (leg_name, subcaption) in enumerate(zip(legendnames, legend_subcaption)):
    if ileg == 0:
        leg.AddEntry(hCorrYield[ileg], legendnames[ileg], 'p')
        # dummy = TObject(1)
        leg.AddEntry(hCorrYield[ileg], subcaption, '')
    else:
        leg.AddEntry(hCorrYield[ileg], legendnames[ileg], 'p')

kNum, kDenom = 0, 1

# Handle systematic uncertainties

# Uncorrelated: raw yield extraction, selection efficiency, pt shape, BR
# Fully correlated: tracking efficiency, luminosity

# as we have a ratio num/denom
# the relative uncertainty square is given by a^2 + b^2 - 2 rho*a*b
# note the "-" due to the ratio

graphNamesNum = ["gCorrYieldSystYieldExtr", "gCorrYieldSystSelEff", "gCorrYieldSystTrEff",
                 "gCorrYieldSystPIDEff", "gCorrYieldSystPtShape", "gCorrYieldSystFD",
                 "gCorrYieldSystNorm"]

graphNamesDenom = ["gCrossSectionSystYieldExtr", "gCrossSectionSystSelEff", "gCrossSectionSystTrEff",
                   "gCrossSectionSystPIDEff", "gCrossSectionSystPtShape", "gCrossSectionSystFD",
                   "gCorrYieldSystNorm"]

graphNamesRatio = ["gRatioSystYieldExtr", "gRatioSystSelEff", "gRatioSystTrEff",
                   "gRatioSystPIDEff", "gRatioSystPtShape", "gRatioSystFD",
                   "gCorrYieldSystNorm"]

gCorrYieldNum, gCorrYieldDenom, gRatio = [], [], []
for graphName in graphNamesNum:
    if 'Norm' not in graphName:
        inputfile = TFile('%s/%s' % (inputdir, inputfilenames[kNum]))
        gCorrYieldNum.append(inputfile.Get(graphName))
    else:
        gCorrYieldNum.append(TGraphErrors(1))
        gCorrYieldNum[-1].SetName(graphName)
        gCorrYieldNum[-1].SetTitle(
            ';#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} #times BR  (#mub GeV^{-1} #it{c})')

    # define gRatio graph
    gRatio.append(TGraphErrors(1))
    gRatio[-1].SetName(graphName)
    gRatio[-1].SetTitle(';#it{p}_{T} (GeV/#it{c}); D^{+}/D^{0}')

for graphName in graphNamesDenom:
    if 'Norm' not in graphName:
        inputfile = TFile('%s/%s' % (inputdir, inputfilenames[kDenom]))
        gCorrYieldDenom.append(inputfile.Get(graphName))
    else:
        gCorrYieldDenom.append(TGraphErrors(1))
        gCorrYieldDenom[-1].SetName(graphName)
        gCorrYieldDenom[-1].SetTitle(';#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} #times BR  (#mub GeV^{-1} #it{c})')


nBins = hCorrYield[0].GetNbinsX()

# configure normalisation TGraphErrors
systNormNum = 0.037
systNormDenom = 0.037
for iPt in range(nBins):
    iPt_denom = iPt+1
    x = gCorrYieldNum[0].GetPointX(iPt)
    y_num = gCorrYieldNum[0].GetPointY(iPt)
    y_denom = gCorrYieldDenom[0].GetPointY(iPt_denom)
    xerr = gCorrYieldNum[0].GetErrorX(iPt)
    rel_yerr_num, rel_yerr_denom = systNormNum, systNormDenom
    # numerator
    gCorrYieldNum[-1].SetPoint(iPt, x, y_num)
    gCorrYieldNum[-1].SetPointError(iPt, xerr, rel_yerr_num*y_num)
    # denominator
    gCorrYieldDenom[-1].SetPoint(iPt_denom, x, y_denom)
    gCorrYieldDenom[-1].SetPointError(iPt_denom, xerr, rel_yerr_denom*y_denom)

# configure gRatio
for iPt in range(nBins):
    iPt_denom = iPt+1
    for iGraph, _ in enumerate(gCorrYieldNum):
        x, xerr = gCorrYieldNum[iGraph].GetPointX(iPt), gCorrYieldNum[iGraph].GetErrorX(iPt)
        yerr_num, yerr_denom = gCorrYieldNum[iGraph].GetErrorY(iPt), gCorrYieldDenom[iGraph].GetErrorY(iPt_denom)
        # compute relative uncertainty
        yerr_num /= hCorrYield[kNum].GetBinContent(iPt+1) # do not scale to 1/BR yet as syst unc in input files are not
        yerr_denom /= hCorrYield[kDenom].GetBinContent(iPt_denom+1)

        if 'TrEff' not in graphNamesNum[iGraph]:
            yerr_ratio = np.sqrt(yerr_num*yerr_num + yerr_denom*yerr_denom)
        elif 'Norm' not in graphNamesNum[iGraph]:
            yerr_ratio = np.sqrt(yerr_num*yerr_num + yerr_denom*yerr_denom)
        else:
            yerr_ratio = np.abs(yerr_num - yerr_denom)
        # if iPt == 0:
        #     print('y err num and denom')
        #     print(yerr_num)
        #     print(yerr_denom)
        #     print(yerr_ratio)
        #     print(' ')
        y_dummy = -1
        gRatio[iGraph].SetPoint(iPt, x, y_dummy)
        gRatio[iGraph].SetPointError(iPt, xerr, yerr_ratio)


### DEBUG
# for graph in gRatio:
#     print('Errors in first pT bin of ratio')
#     print(graph.GetErrorY(0))

# CAREFUL: need to divide by BR
BR_dplus = 8.98 / 100
BR_dzero = 3.95 / 100

hCorrYield[0].Scale(1/BR_dplus)
hCorrYield[1].Scale(1/BR_dzero)

gCorrYield[0].Scale(1/BR_dplus)
gCorrYield[1].Scale(1/BR_dzero)


hRatio = ComputeRatioDiffBins(hCorrYield[0], hCorrYield[1])
hRatio.SetName('hRatio')
hRatio.SetTitle(';#it{p}_{T} (GeV/#it{c}); D^{+}/D^{0})')
SetObjectStyle(hRatio, linecolor=kRed+2, markercolor=kRed+2, markerstyle=kOpenCircle)

# set positions in TGraphErrors
for iPt in range(nBins):
    for graph in gRatio:
        x, xerr = graph.GetPointX(iPt), graph.GetErrorX(iPt)
        y = hRatio.GetBinContent(iPt+1)
        relSyst = graph.GetErrorY(iPt)
        graph.SetPoint(iPt, x, y)
        # if iPt == 0:
        #     print('relSyst')
        #     print(y)
        #     print(relSyst*y)
        graph.SetPointError(iPt, xerr, relSyst*y)

# compute total systematic uncertainty
gRatioTot = TGraphErrors(1)
gRatioTot.SetName('gRatioSystTot')
gRatioTot.SetTitle(';#it{p}_{T} (GeV/#it{c}); D^{+}/D^{0}')

for iPt in range(nBins):
    parsed = False
    yerr_tot = 0
    for graph in gRatio:
        if not parsed:
            x, y = graph.GetPointX(iPt), graph.GetPointY(iPt)
            xerr_tot = graph.GetErrorX(iPt)
            parsed = True
        yerr_tot += graph.GetErrorY(iPt)**2
    yerr_tot = np.sqrt(yerr_tot)
    gRatioTot.SetPoint(iPt, x, y)
    gRatioTot.SetPointError(iPt, xerr_tot, yerr_tot)


ptmin = 2  # hCorrYield[0].GetBinLowEdge(1)+2
ptmax = hCorrYield[0].GetBinLowEdge(hCorrYield[0].GetNbinsX())+hCorrYield[0].GetBinWidth(hCorrYield[0].GetNbinsX())

# lineatone = TLine(ptmin, 1., ptmax, 1.)
# lineatone.SetLineWidth(1)
# lineatone.SetLineColor(kBlack)
# lineatone.SetLineStyle(9)

cCorrYield = TCanvas('cCorrYield', '', 1000, 500)
cCorrYield.Divide(2, 1)
cCorrYield.cd(1).DrawFrame(ptmin, 1.e-1, ptmax, 1.e+4,
                           ';#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} (#mub GeV^{-1} #it{c})')
cCorrYield.cd(1).SetLogy()
# lineatone.Draw('same')
for iFile in range(len(inputfilenames)):
    gCorrYield[iFile].Draw('2')
    hCorrYield[iFile].Draw('same')
leg.Draw()


# ALICE header
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       -0.96 < #it{y} < 0.04')
latex.SetTextSize(0.02)
# latex.DrawLatex(0.22, 0.25, '#pm 1.7% BR unc. not shown')
# latex.DrawLatex(0.22, 0.2, '#pm 2.1% lumi. unc. not shown')

hCorrYieldRatio.append(hRatio)

gSystTot = []
gSystTot.append(gRatioTot)

cCorrYield.cd(2).DrawFrame(ptmin, 0.02, ptmax, 1.2,
                           ';#it{p}_{T} (GeV/#it{c}); D^{+} / D^{0}')

for iFile in range(len(hCorrYieldRatio)): #range(len(inputfilenames)):
    if iFile == 1:
        continue
    SetObjectStyle(hCorrYieldRatio[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                    markerstyle=markers_ratio[iFile]) #, markersize=1.5)
    hCorrYieldRatio[iFile].SetDirectory(0)
    # gCorrYield[iFile].Draw('2')
    gSystTot[iFile].Draw('2')
    SetObjectStyle(gSystTot[iFile], linecolor=colors[iFile], fillstyle=0)
    hCorrYieldRatio[iFile].Draw('same')

legRatio.Draw()

outFile = TFile('/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/comparisons/Ratio_Dplus_D0_FD_pPb.root', 'recreate')

hRatio.Write()
for graph in gRatio:
    graph.Write()
gRatioTot.Write()

cCorrYield.Write()

outFile.Close()
