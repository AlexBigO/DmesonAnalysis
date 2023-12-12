"""
Script to integrate quantity over pt (value, statistical unc, systematic unc)

TODO: make it more general (usable for cross section and Raa)
"""

import numpy as np
import ctypes
from ROOT import gROOT, TFile, TGraphErrors, TObject, TH1F, TCanvas, TLatex, TLegend
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor


SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.65, 0.1, 0.88, 0.35)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)

# kNum, kDenom = 0, 1


# inputFilePPbName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/VisibleCrossSection/VisibleCrossSection_Dplus_FD_pPb_pt_2_16.root'
# inputFilePPName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/VisibleCrossSection/VisibleCrossSection_Dplus_FD_pp_pt_2_16.root'
# inputFileNames = [inputFilePPbName, inputFilePPName]
# APb = 208

# inputFile, hVisCrossSec, gSystTot = [], [], []
# for fileName in inputFileNames:
#     inputFile.append(TFile.Open(fileName))
#     hVisCrossSec.append(inputFile[-1].Get('hVisCrossSec'))
#     gSystTot.append(inputFile[-1].Get('gVisCrossSecTotSys'))













###########################################
# Handle values + statistical uncertainties
###########################################

### Cross section [this is not a mean value as we integrate d^2sigma / dy dpt]

inputFileCrossSectionName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/wptweights_cent/cross_section_divided_by_BR.root'
inputFileCrossSection = TFile.Open(inputFileCrossSectionName)


hCrossSection = inputFileCrossSection.Get('hCorrYield')
nBinsCrossSection = hCrossSection.GetNbinsX()

# unc = 0
# for iPt in range(nBinsCrossSection):
#     print(hCrossSection.GetBinError(iPt+1))
#     unc += (hCrossSection.GetBinError(iPt+1) *  hCrossSection.GetBinWidth(iPt+1) )**2
# print(' ')
# print(np.sqrt(unc) / nBinsCrossSection)

err = ctypes.c_double()
crossSection = hCrossSection.IntegralAndError(1, nBinsCrossSection, err, 'width')
crossSectionUnc = err.value

axisTitleCrossSection = '; y; d#sigma/d#it{y} (#mub)'

hPtIntCrossSection = TH1F('hPtIntCrossSection', axisTitleCrossSection, 1, -0.96, 0.04)
hPtIntCrossSection.SetBinContent(1, crossSection)
hPtIntCrossSection.SetBinError(1, crossSectionUnc)
hPtIntCrossSection.SetDirectory(0)
hPtIntCrossSection.Draw()

leg.AddEntry(hPtIntCrossSection, 'Non-prompt D^{+}', 'p')




## Raa [this is a mean value over a pT range]

inputFileRaaName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/RpPb_nonprompt_Dplus.root'
inputFileRaa = TFile.Open(inputFileRaaName)


hRaa = inputFileRaa.Get('hRpPb')
nBinsRaa = hRaa.GetNbinsX()

deltaPtRaa = hRaa.GetBinLowEdge(nBinsRaa) + hRaa.GetBinWidth(nBinsRaa) - hRaa.GetBinLowEdge(1)
meanRaa, meanRaaUnc = 0, 0
for iPt in range(1, nBinsRaa+1):
    meanRaa += hRaa.GetBinContent(iPt) * hRaa.GetBinWidth(iPt)
    meanRaaUnc += hRaa.GetBinError(iPt)**2 * hRaa.GetBinWidth(iPt)**2

meanRaa /= deltaPtRaa
meanRaaUnc = np.sqrt(meanRaaUnc) / deltaPtRaa

hPtIntRaa = TH1F('hPtIntRaa', '; y; #it{p}_{T} integrated R_{pPb}', 1, -0.96, 0.04)
SetObjectStyle(hPtIntCrossSection, linecolor=kRed+2, markercolor=kRed+2,
                   markerstyle=kOpenCircle)
hPtIntRaa.SetBinContent(1, meanRaa)
hPtIntRaa.SetBinError(1, meanRaaUnc)


###########################################
# Handle systematic uncertainties
#
# raw yield syst: uncorrelated between pT bins
# others: fully correlated between pT bins
###########################################

### Cross section

graphNamesCrossSection = ["gCorrYieldSystYieldExtr", "gCorrYieldSystSelEff", "gCorrYieldSystTrEff",
                "gCorrYieldSystPIDEff", "gCorrYieldSystPtShape", "gCorrYieldSystFD"]

graphNamesPtIntCrossSection = ["gPtIntCrossSectionSystYieldExtr", "gPtIntCrossSectionSystSelEff", "gPtIntCrossSectionSystTrEff",
                "gPtIntCrossSectionSystPIDEff", "gPtIntCrossSectionSystPtShape", "gPtIntCrossSectionSystFD"]

gCrossSection, gPtIntCrossSection = [], []

for graphName in graphNamesCrossSection:
    gCrossSection.append(inputFileCrossSection.Get(graphName))
for graphName in graphNamesPtIntCrossSection:
    gPtIntCrossSection.append(TGraphErrors(1))
    gPtIntCrossSection[-1].SetName(graphName)
    gPtIntCrossSection[-1].SetTitle(axisTitleCrossSection)

for iGraph, graph in enumerate(gCrossSection):
    # relSyst = 0
    yerr = 0
    y_value = crossSection # it is the value of the pt integrated CrossSection
    x, y = 0, y_value
    xerr = 0.4
    if 'YieldExtr' not in graphNamesCrossSection[iGraph]: # fully correlated: linear sum
        for iPt in range(nBinsCrossSection):
            # relSyst += graph.GetErrorY(iPt)/graph.GetPointY(iPt)*hCrossSection.GetBinWidth(iPt+1)
            yerr += graph.GetErrorY(iPt)*hCrossSection.GetBinWidth(iPt+1)
    else: # uncorrelated: quadratic sum
        for iPt in range(nBinsCrossSection):
            # relSyst += (graph.GetErrorY(iPt)/graph.GetPointY(iPt)*hCrossSection.GetBinWidth(iPt+1))* \
            #             (graph.GetErrorY(iPt)/graph.GetPointY(iPt)*hCrossSection.GetBinWidth(iPt+1))
            yerr += (graph.GetErrorY(iPt)*hCrossSection.GetBinWidth(iPt+1))**2
        yerr = np.sqrt(yerr)

    # relSyst /= 22

    gPtIntCrossSection[iGraph].SetPoint(0, x, y)
    gPtIntCrossSection[iGraph].SetPointError(0, xerr, yerr)

yerr = 0
for graph in gPtIntCrossSection:
    yerr += graph.GetErrorY(0)*graph.GetErrorY(0)
yerr = np.sqrt(yerr)
gPtIntCrossSectionSystTot = TGraphErrors(1)
SetObjectStyle(gPtIntCrossSectionSystTot, linecolor=kRed+2, markercolor=kRed+2, fillstyle=0)
gPtIntCrossSectionSystTot.SetName('gPtIntCrossSectionSystTot')
gPtIntCrossSectionSystTot.SetTitle(axisTitleCrossSection)
y = crossSection # to set to pt int Raa value
gPtIntCrossSectionSystTot.SetPoint(0, 0.04-0.5, y)
gPtIntCrossSectionSystTot.SetPointError(0, 0.4, yerr)

lumiUnc = 0.037*crossSection
brUnc = 0.017*crossSection

print(f'Cross section sytematics: {gPtIntCrossSectionSystTot.GetErrorY(0)}')

print(f'Cross section = {crossSection:.1f} +- {crossSectionUnc:.1f} (stat) \
        +- {gPtIntCrossSectionSystTot.GetErrorY(0):.1f} (syst) +- {lumiUnc:.1f} (lumi) \
        +- {brUnc:.1f} (BR)')

cPtIntCrossSection = TCanvas('cPtIntCrossSection', '', 800, 800)
cPtIntCrossSection.DrawFrame(-0.96, 0.8*crossSection, 0.04, 1.2*crossSection, axisTitleCrossSection)
gPtIntCrossSectionSystTot.Draw('2')
hPtIntCrossSection.Draw('same')

latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       2 < #it{p}_{T} < 24 GeV/#it{c}')
latex.SetTextSize(0.025)
latex.DrawLatex(0.22, 0.25, '#pm 1.7% BR unc. not shown')
latex.DrawLatex(0.22, 0.2, '#pm 3.7% lumi. unc. not shown')

leg.Draw()

outFileCrossSectionName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/CrossSection_nonprompt_Dplus_pt_integrated.root'
outFileCrossSection = TFile.Open(outFileCrossSectionName, 'recreate')
hPtIntCrossSection.Write()
for graph in gPtIntCrossSection:
    graph.Write()
gPtIntCrossSectionSystTot.Write()

cPtIntCrossSection.Write()
outFileCrossSection.Close()

outdir = "/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection"
cPtIntCrossSection.SaveAs(f'{outdir}/PtIntegratedCrossSection_Dplus_FD.pdf')


### Raa

graphNamesRaa = ["gRaaSystYieldExtr", "gRaaSystSelEff", "gRaaSystTrEff",
                "gRaaSystPIDEff", "gRaaSystPtShape", "gRaaSystFD",
                "gRaaSystNorm"]

graphNamesPtIntRaa = ["gPtIntRaaSystYieldExtr", "gPtIntRaaSystSelEff", "gPtIntRaaSystTrEff",
                "gPtIntRaaSystPIDEff", "gPtIntRaaSystPtShape", "gPtIntRaaSystFD",
                "gRaaSystNorm"]

gRaa, gPtIntRaa = [], []

for graphName in graphNamesRaa:
    gRaa.append(inputFileRaa.Get(graphName))
for graphName in graphNamesPtIntRaa:
    gPtIntRaa.append(TGraphErrors(1))
    gPtIntRaa[-1].SetName(graphName)
    gPtIntRaa[-1].SetTitle('; ; #it{p}_{T} integrated R_{pPb}')


for iGraph, graph in enumerate(gRaa):
    relSyst = 0
    y_value = meanRaa # it is the value of the pt integrated Raa
    x, y = 0, y_value
    xerr = 0.4
    if 'YieldExtr' not in graphNamesRaa[iGraph]: # fully correlated: linear sum
        for iPt in range(nBinsRaa):
            if 'Norm' not in graphNamesRaa[iGraph]:
                relSyst += graph.GetErrorY(iPt)/graph.GetPointY(iPt)*hRaa.GetBinWidth(iPt+1)
            else:
                relSyst += graph.GetErrorY(0)/graph.GetPointY(0)*hRaa.GetBinWidth(iPt+1)
    else: # uncorrelated: quadratic sum
        for iPt in range(nBinsRaa):
            relSyst += (graph.GetErrorY(iPt)/graph.GetPointY(iPt)*hRaa.GetBinWidth(iPt+1))* \
                        (graph.GetErrorY(iPt)/graph.GetPointY(iPt)*hRaa.GetBinWidth(iPt+1))
        relSyst = np.sqrt(relSyst)

    relSyst /= deltaPtRaa

    gPtIntRaa[iGraph].SetPoint(0, x, y)
    gPtIntRaa[iGraph].SetPointError(0, xerr, relSyst*y)

yerr = 0
for graph in gPtIntRaa:
    yerr += graph.GetErrorY(0)*graph.GetErrorY(0)
yerr = np.sqrt(yerr)
gPtIntRaaSystTot = TGraphErrors(1)
gPtIntRaaSystTot.SetName('gPtIntRaaSystTot')
gPtIntRaaSystTot.SetTitle('; y; #it{p}_{T} integrated R_{pPb}')
y = meanRaa # to set to pt int Raa value
gPtIntRaaSystTot.SetPoint(0, 0.04-0.5, y)
gPtIntRaaSystTot.SetPointError(0, 0.4, yerr)


cPtIntRaa = TCanvas('cPtIntRaa', '', 800, 800)
cPtIntRaa.DrawFrame(-0.96, 0.8, 0.04, 1.6, '; y; R_{pPb}(y)')

SetObjectStyle(hPtIntRaa, linecolor=kRed+2, markercolor=kRed+2,
                   markerstyle=kOpenCircle)

hPtIntRaa.SetDirectory(0)
gPtIntRaaSystTot.Draw('2')
SetObjectStyle(gPtIntRaaSystTot, linecolor=kRed+2, markercolor=kRed+2, fillstyle=0)
hPtIntRaa.Draw('same')

latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       2 < #it{p}_{T} < 16 GeV/#it{c}')

leg.Draw()


outFileName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/RpPb_nonprompt_Dplus_pt_integrated.root'
outFile = TFile.Open(outFileName, 'recreate')
hPtIntRaa.Write()
for graph in gPtIntRaa:
    graph.Write()
gPtIntRaaSystTot.Write()

cPtIntRaa.Write()
outFile.Close()


outdir = "/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb"
cPtIntRaa.SaveAs(f'{outdir}/PtIntegratedRpPb_Dplus_FD.pdf')

input('Press enter to exit')