import sys
import numpy as np
from ROOT import TCanvas, TFile, TLegend, TLine, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kOpenSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins

compareD0Dplus = False # True

inputdir = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/'

if compareD0Dplus:
    inputfilenames = ['RpPb_np_D0_5TeV.root', 'RpPb_nonprompt_Dplus.root']
    histonames = ['RpPbFDD0', 'hRpPb']
    legendnames = ['Non-prompt D^{0}', 'Non-prompt D^{+}']
    graphnames = ['SystwoNorm', 'gRaaSystTot']
    legend_subcaption = ['(Mingyu, Preliminary)', '']
    outputsuffix = 'D0_Dplus_pPb_FD'

else:
    inputfilenames = ['HEPData_RpPbDplus_prompt.root', 'RpPb_nonprompt_Dplus.root'] # ['RpPb_np_D0_5TeV.root', 'RpPb_nonprompt_Dplus.root']
    histonames = ['hRaa', 'hRpPb'] # ['RpPbFDD0', 'hRpPb']
    graphnames = ['gRaaSystTot', 'gRaaSystTot'] # ['SystwoNorm', 'gRaaSystTot']
    legendnames = ['Prompt D^{+}', 'Non-prompt D^{+}'] # ['Non-prompt D^{0}', 'Non-prompt D^{+}']
    legend_subcaption = legend_subcaption = ['JHEP 12 (2019) 092, 2019', ''] # ['(Mingyu, Preliminary)', '']
    outputsuffix = 'Prompt_Dplus_FD_Dplus' #'D0_Dplus_pPb_FD'



colors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
linecolors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
markers = [kFullDiamond, kFullCircle]



hCorrYield, gCorrYield, hCorrYieldRatio = ([] for _ in range(3))

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

if compareD0Dplus:
    leg = TLegend(0.18, 0.58, 0.42, 0.8)
else:
    leg = TLegend(0.48, 0.58, 0.82, 0.78) #(0.18, 0.58, 0.42, 0.8)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.035)

legRatio = TLegend(0.7, 0.25, 0.9, 0.4)
legRatio.SetFillStyle(0)
legRatio.SetBorderSize(0)
legRatio.SetTextSize(0.04)


for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hCorrYield.append(inputfile.Get(histonames[iFile]))
    # if iFile == 0:
    gCorrYield.append(inputfile.Get(graphnames[iFile]))
    hCorrYield[iFile].SetDirectory(0)
    SetObjectStyle(hCorrYield[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile]) #, markersize=1.5)
    # if iFile == 0:
    SetObjectStyle(gCorrYield[iFile], linecolor=colors[iFile], fillstyle=0)
    # leg.AddEntry(hCorrYield[iFile], legendnames[iFile], 'p')

# legend
for ileg, (leg_name, subcaption) in enumerate(zip(legendnames, legend_subcaption)):
    if ileg == 0:
        leg.AddEntry(hCorrYield[ileg], legendnames[ileg], 'p')
        # dummy = TObject(1)
        leg.AddEntry(hCorrYield[ileg], subcaption, '')
    else:
        leg.AddEntry(hCorrYield[ileg], legendnames[ileg], 'p')


ptmin = 0 # hCorrYield[0].GetBinLowEdge(1)+2
ptmax = hCorrYield[0].GetBinLowEdge(hCorrYield[0].GetNbinsX())+hCorrYield[0].GetBinWidth(hCorrYield[0].GetNbinsX())

# lineatone = TLine(ptmin, 1., ptmax, 1.)
# lineatone.SetLineWidth(1)
# lineatone.SetLineColor(kBlack)
# lineatone.SetLineStyle(9)

cCorrYield = TCanvas('cCorrYield', '', 800, 800)
cCorrYield.DrawFrame(ptmin, 0.1, ptmax, 3.1,
                           ';#it{p}_{T} (GeV/#it{c}); #it{R}_{pPb}')
# cCorrYield.SetLogy()
# lineatone.Draw('same')
for iFile in range(len(inputfilenames)):
    # if iFile == 0:
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
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s} = 5.02 TeV       -0.96 < #it{y}_{cms} < 0.04')
latex.SetTextSize(0.025)
if not compareD0Dplus:
    latex.DrawLatex(0.22, 0.2, '#pm (4.3) 4.1% lumi. unc. not shown for (non)prompt')
else:
    lumiUncPPb = 0.037
    latex.DrawLatex(0.22, 0.2, f'#pm {lumiUncPPb*100:.1f}% lumi. unc. not show')


cCorrYield.SaveAs('%s/RpPb_Comparison_%s.pdf' % (inputdir, outputsuffix))

input('Press enter to exit')
