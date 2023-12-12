import sys
from ROOT import TCanvas, TFile, TLegend, TLine, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kOpenSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins

inputdir = '/home/abigot/AnalysisNonPromptDplus/'
inputfilenames = ['Ratio_Dplus_D0_prompt_pPb.root',
                'Run2pPb5Tev/4_Analysis/5_CrossSection/comparisons/Ratio_Dplus_D0_FD_pPb.root'] #CrossSectionDplus_pPb_nonprompt.root']
histonames = ['hRatio', 'hRatio']
graphnames = ['gRatioSystTot', 'gRatioSystTot']
colors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
linecolors = [kRed+1, kAzure+4, kBlack, kBlue, kRed]
markers = [kFullDiamond, kFullCircle]
markers_ratio = [kOpenDiamond, kOpenCircle]
legendnames = ['Prompt D^{+}/D^{0}', 'Non-prompt D^{+}/D^{0}']
legend_subcaption = ['JHEP 12 (2019) 092, 2019', '']
outputsuffix = 'Dplus_D0'

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.5, 0.58, 0.82, 0.75)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.035)

hRatios, gRatios = [], []

for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hRatios.append(inputfile.Get(histonames[iFile]))

    gRatios.append(inputfile.Get(graphnames[iFile]))
    hRatios[iFile].SetDirectory(0)
    SetObjectStyle(hRatios[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile]) #, markersize=1.5)
    SetObjectStyle(gRatios[iFile], linecolor=colors[iFile], fillstyle=0)

# legend
for ileg, (leg_name, subcaption) in enumerate(zip(legendnames, legend_subcaption)):
    if ileg == 0:
        leg.AddEntry(hRatios[ileg], legendnames[ileg], 'p')
        # dummy = TObject(1)
        leg.AddEntry(hRatios[ileg], subcaption, '')
    else:
        leg.AddEntry(hRatios[ileg], legendnames[ileg], 'p')

ptmin = 0 #hRatios[1].GetBinLowEdge(1)
ptmax = hRatios[0].GetBinLowEdge(hRatios[0].GetNbinsX())+hRatios[0].GetBinWidth(hRatios[0].GetNbinsX())

cRatio = TCanvas('cRatio', '', 800, 800)
cRatio.DrawFrame(ptmin, 0.0, ptmax, 1.5, ';#it{p}_{T} (GeV/#it{c}); D^{+}/D^{0}')

for iFile in range(len(inputfilenames)):
    gRatios[iFile].Draw('2')
    hRatios[iFile].Draw('same')
leg.Draw()


# ALICE header
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
latex.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s_{NN}} = 5.02 TeV       -0.96 < #it{y}_{cms} < 0.04')
latex.SetTextSize(0.025)
latex.DrawLatex(0.22, 0.2, '#pm 1.9% BR unc. not shown')  
# latex.DrawLatex(0.22, 0.2, '#pm 2.1% lumi. unc. not shown')

outdir = "/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/comparisons"
cRatio.SaveAs('%s/Ratio_%s.pdf' % (outdir, outputsuffix))

input('Press enter to exit')
