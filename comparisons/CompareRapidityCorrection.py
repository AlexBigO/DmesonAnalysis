import sys
from ROOT import TCanvas, TFile, TLegend, TLine, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kOpenSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins

inputdir = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/rapidity_correction/'
inputfilenames = ['rapidity_correction_rebinned.root']
histonames = ['hYCorrection_central', 'hYCorrection_min', 'hYCorrection_max']
# graphnames = ['gCorrYieldSystTot', 'gCorrYieldSystTot']
colors = [kBlack, kRed+1, kAzure+4, kBlue, kRed]
linecolors = [kBlack, kRed+1, kAzure+4, kBlue, kRed]
markers = [kFullCircle, kFullDiamond, kFullCircle]
legendnames = ['central', 'min', 'max']
# legend_subcaption = ['JHEP 12 (2019) 092, 2019', '']
outputsuffix = 'Dplus_pPb'


hYCorrection, gCorrYield, hYCorrectionRatio = ([] for _ in range(3))

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.42, 0.5, 0.82, 0.75)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

legRatio = TLegend(0.7, 0.2, 0.9, 0.4)
legRatio.SetFillStyle(0)
legRatio.SetBorderSize(0)
legRatio.SetTextSize(0.04)

# # pp ratio
# pp_ratio_file = TFile('%s/%s' % (inputdir, pp_ratio_filename))
# hRatio_pp = pp_ratio_file.Get('hRatioDplus')
# SetObjectStyle(hRatio_pp, linecolor=kAzure+2, markercolor=kAzure+2,
#                    markerstyle=kOpenDiamond)

inputfile = TFile('%s/%s' % (inputdir, inputfilenames[0]))
for iFile in range(3):
    hYCorrection.append(inputfile.Get(histonames[iFile]))
    # gCorrYield.append(inputfile.Get(graphnames[iFile]))
    hYCorrection[iFile].SetDirectory(0)
    SetObjectStyle(hYCorrection[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile]) #, markersize=1.5)
    # SetObjectStyle(gCorrYield[iFile], linecolor=colors[iFile], fillstyle=0)
    # leg.AddEntry(hYCorrection[iFile], legendnames[iFile], 'p')
    hYCorrectionRatio.append(hYCorrection[iFile].Clone("hYCorrection%d" % iFile))
    hYCorrectionRatio[iFile].SetDirectory(0)
    hYCorrectionRatio[iFile].Divide(hYCorrection[iFile], hYCorrection[0])

# legend
for ileg, leg_name in enumerate(legendnames):
    leg.AddEntry(hYCorrection[ileg], legendnames[ileg], 'p')




hRatio = ComputeRatioDiffBins(hYCorrection[1], hYCorrection[0])
hRatio.SetName('hRatio')
hRatio.SetTitle('Non-prompt / prompt ratio')
SetObjectStyle(hRatio, linecolor=kRed+2, markercolor=kRed+2,
                   markerstyle=kOpenCircle)

ptmin = hYCorrection[0].GetBinLowEdge(1)
ptmax = hYCorrection[0].GetBinLowEdge(hYCorrection[0].GetNbinsX())+hYCorrection[0].GetBinWidth(hYCorrection[0].GetNbinsX())


cCorrYield = TCanvas('cCorrYield', '', 1000, 500)
cCorrYield.Divide(2, 1)
cCorrYield.cd(1).DrawFrame(ptmin, 0.980, ptmax, 1.,
                           ';#it{p}_{T} (GeV/#it{c}); rapidity correction factor')
# cCorrYield.cd(1).SetLogy()
# lineatone.Draw('same')
for iHist in range(len(histonames)):
    # gCorrYield[iFile].Draw('2')
    hYCorrection[iHist].Draw('same')
leg.Draw()


# ALICE header
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s} = 5.02 TeV       -0.96 < #it{y} < 0.04')



cCorrYield.cd(2).DrawFrame(ptmin, 0.99, ptmax, 1.01,
                           ';#it{p}_{T} (GeV/#it{c}); ratio')
for iHist in range(len(histonames)):
    if iHist == 0:
        continue
    # gCorrYield[iFile].Draw('2')
    hYCorrectionRatio[iHist].Draw('same')
# cCorrYield.cd(2).SetLogy()
# hRatio.Draw('same')

# legRatio.AddEntry(hRatio_pp, 'D^{+} in pp', 'p')
# legRatio.AddEntry(hRatio, 'D^{+} in pPb', 'p')
# legRatio.Draw()

lineatone = TLine(ptmin, 1., ptmax, 1.)
lineatone.SetLineWidth(1)
lineatone.SetLineColor(kBlack)
lineatone.SetLineStyle(9)

cCorrYield.SaveAs('%s/RapidityCorrection_%s.pdf' % (inputdir, outputsuffix))

input('Press enter to exit')
