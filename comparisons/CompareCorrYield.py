import sys
from ROOT import TCanvas, TFile, TLegend, TLine, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kOpenSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins

inputdir = '/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/'
pp_ratio_filename = 'D0DplusDs_NonPromptOverPrompt_ratios_pp5TeV.root'
inputfilenames = ['CrossSectionDplus_pPb_prompt.root', 'wptweights_cent/cross_section_divided_by_BR.root'] #CrossSectionDplus_pPb_nonprompt.root']
histonames = ['hCorrYield', 'hCorrYield']
graphnames = ['gCorrYieldSystTot', 'gCorrYieldSystTot']
colors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
linecolors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
markers = [kFullDiamond, kFullCircle]
markers_ratio = [kOpenDiamond, kOpenCircle]
legendnames = ['Prompt D^{+}', 'Non-prompt D^{+}']
legend_subcaption = ['JHEP 12 (2019) 092, 2019', '']
outputsuffix = 'Dplus_pPb'


hCorrYield, gCorrYield, hCorrYieldRatio = ([] for _ in range(3))

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.42, 0.5, 0.82, 0.75)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

legRatio = TLegend(0.7, 0.2, 0.9, 0.4)
legRatio.SetFillStyle(0)
legRatio.SetBorderSize(0)
legRatio.SetTextSize(0.04)

# pp ratio
pp_ratio_file = TFile('%s/%s' % (inputdir, pp_ratio_filename))
hRatio_pp = pp_ratio_file.Get('hRatioDplus')
SetObjectStyle(hRatio_pp, linecolor=kAzure+2, markercolor=kAzure+2,
                   markerstyle=kOpenDiamond)

for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hCorrYield.append(inputfile.Get(histonames[iFile]))
    gCorrYield.append(inputfile.Get(graphnames[iFile]))
    hCorrYield[iFile].SetDirectory(0)
    SetObjectStyle(hCorrYield[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile]) #, markersize=1.5)
    SetObjectStyle(gCorrYield[iFile], linecolor=colors[iFile], fillstyle=0)
    # leg.AddEntry(hCorrYield[iFile], legendnames[iFile], 'p')
    # hCorrYieldRatio.append(hCorrYield[iFile].Clone("hCorrYield%d" % iFile))
    # hCorrYieldRatio[iFile].SetDirectory(0)
    # hCorrYieldRatio[iFile].Divide(hCorrYield[iFile], hCorrYield[0])

# legend
for ileg, (leg_name, subcaption) in enumerate(zip(legendnames, legend_subcaption)):
    if ileg == 0:
        leg.AddEntry(hCorrYield[ileg], legendnames[ileg], 'p')
        # dummy = TObject(1)
        leg.AddEntry(hCorrYield[ileg], subcaption, '')
    else:
        leg.AddEntry(hCorrYield[ileg], legendnames[ileg], 'p')




hRatio = ComputeRatioDiffBins(hCorrYield[1], hCorrYield[0])
hRatio.SetName('hRatio')
hRatio.SetTitle('Non-prompt / prompt ratio')
SetObjectStyle(hRatio, linecolor=kRed+2, markercolor=kRed+2,
                   markerstyle=kOpenCircle)

ptmin = hCorrYield[0].GetBinLowEdge(1)+2
ptmax = hCorrYield[0].GetBinLowEdge(hCorrYield[0].GetNbinsX())+hCorrYield[0].GetBinWidth(hCorrYield[0].GetNbinsX())

# lineatone = TLine(ptmin, 1., ptmax, 1.)
# lineatone.SetLineWidth(1)
# lineatone.SetLineColor(kBlack)
# lineatone.SetLineStyle(9)



cCorrYieldNoRatio = TCanvas('cCorrYieldNoRatio', '', 800, 800)
cCorrYieldNoRatio.DrawFrame(ptmin, 1.e-1, ptmax, 1.e+5,
                           ';#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} (#mub GeV^{-1} #it{c})')

cCorrYieldNoRatio.cd(1).SetLogy()
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
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       -0.96 < #it{y}_{cms} < 0.04')
latex.SetTextSize(0.02)
latex.DrawLatex(0.22, 0.25, '#pm 1.7% BR unc. not shown')
latex.DrawLatex(0.22, 0.2, '#pm 3.7% lumi. unc. not shown')

cCorrYieldNoRatio.SaveAs('%s/CorrYieldComparisonNoRatio_%s.pdf' % (inputdir, outputsuffix))


########################################


cCorrYield = TCanvas('cCorrYield', '', 1000, 500)
cCorrYield.Divide(2, 1)
cCorrYield.cd(1).DrawFrame(ptmin, 1.e-1, ptmax, 1.e+5,
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
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       -0.96 < #it{y}_{cms} < 0.04')
latex.SetTextSize(0.02)
latex.DrawLatex(0.22, 0.25, '#pm 1.7% BR unc. not shown')
latex.DrawLatex(0.22, 0.2, '#pm 3.7% lumi. unc. not shown')


hCorrYieldRatio.append(hRatio_pp)
hCorrYieldRatio.append(hRatio)

cCorrYield.cd(2).DrawFrame(ptmin, 0.02, ptmax, 0.16,
                           ';#it{p}_{T} (GeV/#it{c}); d^{2}#sigma/d#it{p}_{T}d#it{y} (non-prompt) / d^{2}#sigma/d#it{p}_{T}d#it{y} (prompt')
legRatio.AddEntry(hRatio_pp, 'D^{+} in pp', 'p')
legRatio.AddEntry(hRatio, 'D^{+} in pPb', 'p')
for iFile in range(len(inputfilenames)):
    # if iFile == 0:
    #     continue
    SetObjectStyle(hCorrYieldRatio[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers_ratio[iFile]) #, markersize=1.5)
    hCorrYieldRatio[iFile].SetDirectory(0)
    gCorrYield[iFile].Draw('2')
    hCorrYieldRatio[iFile].Draw('same')
# cCorrYield.cd(2).SetLogy()
# hRatio.Draw('same')
# hRatio_pp.Draw('same')
legRatio.Draw()

cCorrYield.SaveAs('%s/CorrYieldComparison_%s.pdf' % (inputdir, outputsuffix))

input('Press enter to exit')
