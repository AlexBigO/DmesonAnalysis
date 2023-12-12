"""
Script to integrate quantity over pt (value, statistical unc, systematic unc)

TODO: make it more general (usable for cross section and Raa)
"""

import numpy as np
import ctypes
from ROOT import gROOT, TFile, TGraphErrors, TObject, TH1F, TCanvas, TLatex, TLegend,TH1D
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor

doWoSystLumi = True
doExtrapolated = False

if doWoSystLumi and doExtrapolated:
    print('ERROR: extrapolated has to be with lumi unc.')

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.65, 0.1, 0.88, 0.35)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)

kNum, kDenom = 0, 1

if doExtrapolated:
    inputFilePPbName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/Extrapolation/NonPromptDplus_PtIntegr2-24_pPb_FONLLextrap.root'
    inputFilePPName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/Extrapolation/NonPromptDplus_PtIntegr2-16_FONLLextrap_pp_published.root'
else:
    inputFilePPbName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/VisibleCrossSection/VisibleCrossSection_Dplus_FD_pPb_pt_2_16.root'
    inputFilePPName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/VisibleCrossSection/VisibleCrossSection_Dplus_FD_pp_pt_2_16.root'
inputFileNames = [inputFilePPbName, inputFilePPName]
APb = 208

inputFile, hVisCrossSec, gSystTot, gSystExtrap = [], [], [], []
for fileName in inputFileNames:
    inputFile.append(TFile.Open(fileName))
    if doExtrapolated:
        hVisCrossSec.append(inputFile[-1].Get('hPtIntCrossSecStat'))
    else:
        hVisCrossSec.append(inputFile[-1].Get('hVisCrossSec'))
    if not doWoSystLumi:
        if doExtrapolated:
            gSystTot.append(inputFile[-1].Get('gPtIntCrossSecDataSyst'))
            gSystExtrap.append(inputFile[-1].Get('gPtIntCrossSecExtrapSyst'))
        else:
            gSystTot.append(inputFile[-1].Get('gVisCrossSecTotSys'))
    else:
        gSystTot.append(inputFile[-1].Get('gVisCrossSecSysWoBRAndLumi'))

hVisibleRaa_old = hVisCrossSec[kNum].Clone('hVisibleRaa')
axisTitle = ';y; #it{R}_{pPb}'
hVisibleRaa_old.SetTitle(axisTitle)
hVisibleRaa_old.Divide(hVisCrossSec[kDenom])
hVisibleRaa_old.Scale(1./APb)

hVisibleRaa = TH1D('hVisibleRaa', axisTitle, 1, -0.96, 0.04)
hVisibleRaa.SetBinContent(1, hVisibleRaa_old.GetBinContent(1))
hVisibleRaa.SetBinError(1, hVisibleRaa_old.GetBinError(1))
hVisibleRaa.SetDirectory(0)
SetObjectStyle(hVisibleRaa, linecolor=kRed+2, markercolor=kRed+2,
                   markerstyle=kOpenCircle)


visibleRaa = hVisibleRaa.GetBinContent(1)

# Handle systematics
gVisibleRaaSystTot = TGraphErrors(1)
gVisibleRaaSystTot.SetName('gVisibleRaaSystTot')
gVisibleRaaSystTot.SetTitle(axisTitle)
SetObjectStyle(gVisibleRaaSystTot, linecolor=kRed+2, markercolor=kRed+2, fillstyle=0)
relSystTot = 0
for graph in gSystTot:
    relSystTot += (graph.GetErrorY(0) / graph.GetPointY(0))**2

relSystTot = np.sqrt(relSystTot)

print(hVisibleRaa.GetBinContent(1))

gVisibleRaaSystTot.SetPoint(0, hVisibleRaa.GetBinCenter(1), visibleRaa)
gVisibleRaaSystTot.SetPointError(0, 0.3, relSystTot * visibleRaa)




cVisibleRaa = TCanvas('cVisibleRaa', '', 800, 800)
cVisibleRaa.DrawFrame(-0.96, 0.7*visibleRaa, 0.04, 1.3*visibleRaa, axisTitle)

gVisibleRaaSystTot.Draw('2')
hVisibleRaa.Draw('same')


latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
if doExtrapolated:
    latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       0 < #it{p}_{T} < 24 GeV/#it{c}')
else:
    latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       2 < #it{p}_{T} < 16 GeV/#it{c}')
if doWoSystLumi:
    lumiUncPP = 0.021
    lumiUncPPb = 0.037
    lumiUnc = np.sqrt(lumiUncPPb*lumiUncPPb + lumiUncPP*lumiUncPP)*100
    latex.SetTextSize(0.03)
    latex.DrawLatex(0.22, 0.2, f'#pm {lumiUnc:.1f}% lumi. unc. not shown')

leg.AddEntry(hVisibleRaa, 'Non-prompt D^{+}', 'p')
leg.Draw()


if doWoSystLumi:
    outFileName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/VisibleRpPb_Dplus_FD_pt_2_16_woLumi.root'
else:
    if doExtrapolated:
        outFileName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/VisibleRpPb_Dplus_FD_pt_0_24.root'
    else:
        outFileName = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/VisibleRpPb_Dplus_FD_pt_2_16.root'
outFile = TFile.Open(outFileName, 'recreate')
hVisibleRaa.Write()
gVisibleRaaSystTot.Write()
gVisibleRaaSystTot.Write()

cVisibleRaa.Write()
outFile.Close()

outFileNamePDF = outFileName.replace('.root', '.pdf')
cVisibleRaa.SaveAs(outFileNamePDF)

input('Press enter to exit')
