'''
python script to compare measured cross sections with FONLL
run: python CompareCrossSecToFONLL.py FONLL.root outFileName.pdf [--Dplus] [--Ds] [--prompt CrossSecPrompt.root]
                                                                 [--FD CrossSecFD.root] [--logx] [--syst]
Either Dplus or Ds must be chosen
Either prompt or FD (or both) must be set
'''

import sys
import argparse
import numpy as np
from ROOT import TFile, TCanvas, TGraphAsymmErrors, TLegend, TLine, TLatex, TPad # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ScaleGraph, DivideGraphByHisto, ComputeRatioGraph #pylint: disable=wrong-import-position,import-error

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('FONLLFileName', metavar='text', default='FONLL.root', help='root file FONLL predictions')
parser.add_argument('outFileName', metavar='text', default='outFile.pdf', help='pdf output file')
parser.add_argument('--Dplus', action='store_true', default=False, help='enable comparison for D+')
parser.add_argument('--Ds', action='store_true', default=False, help='enable comparison for Ds')
parser.add_argument('--prompt', metavar='text', default=None,
                    help='enable comparison for prompt D and pass Cross section file name')
parser.add_argument('--FD', metavar='text', default=None,
                    help='enable comparison for FD D and pass Cross section file name')
parser.add_argument('--logx', action='store_true', default=False, help='enable log scale for x axis')
parser.add_argument('--syst', action='store_true', default=False, help='plot also graph of syst unc')
args = parser.parse_args()

if not args.Dplus and not args.Ds:
    print('ERROR: you should enable the comparison for either D+ or Ds! Exit')
    sys.exit()
elif args.Dplus and args.Ds:
    print('ERROR: you cannot enable the comparison for both D+ and Ds! Exit')
    sys.exit()

if args.Dplus:
    BR = 1.
    ApPb = 208
    mesonName = 'D^{+}'
elif args.Ds:
    BR = 0.0227
    mesonName = 'D_{s}^{+}'

SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.6, maxdigits=2, optstat=0)

if args.prompt:
    infilePrompt = TFile.Open(args.prompt)
    hCorrYieldPrompt = infilePrompt.Get('hCorrYield')
    if args.syst:
        gCorrYieldPrompt = infilePrompt.Get('gCorrYieldSystTot')
    if not hCorrYieldPrompt:
        hCorrYieldPrompt = infilePrompt.Get('histoSigmaCorr')
        hCorrYieldPrompt.Scale(1.e-6 / BR)
        hCorrYieldPrompt.SetStats(0)
        if args.syst:
            gCorrYieldPrompt = infilePrompt.Get('gSigmaCorr')
            gCorrYieldPrompt.RemovePoint(0)
            ScaleGraph(gCorrYieldPrompt, 1.e-6 / BR)
            for iPt in range(hCorrYieldPrompt.GetNbinsX()):
                gCorrYieldPrompt.SetPointEXhigh(iPt, hCorrYieldPrompt.GetBinWidth(iPt+1)/2.5)
                gCorrYieldPrompt.SetPointEXlow(iPt, hCorrYieldPrompt.GetBinWidth(iPt+1)/2.5)
    hCorrYieldPrompt.SetName('hCorrYieldPrompt')
    hCorrYieldPrompt.SetDirectory(0)
    SetObjectStyle(hCorrYieldPrompt, color=kBlack, markerstyle=kFullCircle)
    if args.syst:
        SetObjectStyle(gCorrYieldPrompt, color=kBlack, fillstyle=0)
        gCorrYieldPrompt.SetName('gCorrYieldPrompt')
    infilePrompt.Close()

    ptLimitsPrompt = np.array(hCorrYieldPrompt.GetXaxis().GetXbins(), 'd')
    ptMinPrompt = list(ptLimitsPrompt)[0]
    ptMaxPrompt = list(ptLimitsPrompt)[-1]
    ptMin = ptMinPrompt
    ptMax = ptMaxPrompt
    sigmaMin = hCorrYieldPrompt.GetMinimum()*0.2
    sigmaMax = hCorrYieldPrompt.GetMaximum()*5

    infileFONLL = TFile.Open(args.FONLLFileName)
    if args.Dplus:
        hFONLLPromptCentral = infileFONLL.Get('hDpluskpipipred_central')
        hFONLLPromptMin = infileFONLL.Get('hDpluskpipipred_min')
        hFONLLPromptMax = infileFONLL.Get('hDpluskpipipred_max')
    elif args.Ds:
        hFONLLPromptCentral = infileFONLL.Get('hDsPhipitoKkpipred_central')
        hFONLLPromptMin = infileFONLL.Get('hDsPhipitoKkpipred_min')
        hFONLLPromptMax = infileFONLL.Get('hDsPhipitoKkpipred_max')
    hFONLLPromptCentral = hFONLLPromptCentral.Rebin(
        hCorrYieldPrompt.GetNbinsX(), 'hFONLLPromptCentral', ptLimitsPrompt)
    hFONLLPromptMin = hFONLLPromptMin.Rebin(hCorrYieldPrompt.GetNbinsX(), 'hFONLLPromptMin', ptLimitsPrompt)
    hFONLLPromptMax = hFONLLPromptMax.Rebin(hCorrYieldPrompt.GetNbinsX(), 'hFONLLPromptMax', ptLimitsPrompt)
    hFONLLPromptCentral.Scale(1.e-6 / BR , 'width')
    hFONLLPromptMin.Scale(1.e-6 / BR , 'width')
    hFONLLPromptMax.Scale(1.e-6 / BR , 'width')
    hFONLLPromptCentral.SetDirectory(0)
    hFONLLPromptMin.SetDirectory(0)
    hFONLLPromptMax.SetDirectory(0)
    hFONLLPromptCentral.SetStats(0)
    hFONLLPromptMin.SetStats(0)
    hFONLLPromptMax.SetStats(0)
    gFONLLPrompt = TGraphAsymmErrors(0)
    for iPt in range(hFONLLPromptCentral.GetNbinsX()):
        gFONLLPrompt.SetPoint(iPt, hFONLLPromptCentral.GetBinCenter(iPt+1), hFONLLPromptCentral.GetBinContent(iPt+1))
        gFONLLPrompt.SetPointError(iPt, hFONLLPromptCentral.GetBinWidth(iPt+1)/2,
                                   hFONLLPromptCentral.GetBinWidth(iPt+1)/2,
                                   hFONLLPromptCentral.GetBinContent(iPt+1)-hFONLLPromptMin.GetBinContent(iPt+1),
                                   hFONLLPromptMax.GetBinContent(iPt+1)-hFONLLPromptCentral.GetBinContent(iPt+1))

    gFONLLPromptUnc = ComputeRatioGraph(gFONLLPrompt, gFONLLPrompt, False)
    SetObjectStyle(hFONLLPromptCentral, color=kRed+1, markerstyle=0, fillstyle=0)
    SetObjectStyle(gFONLLPrompt, color=kRed+1, fillalpha=0.2, fillstyle=1000)
    SetObjectStyle(gFONLLPromptUnc, color=kRed+1, fillalpha=0.2, fillstyle=1000)
    infileFONLL.Close()

    hRatioPromptOverFONLL = hCorrYieldPrompt.Clone('hRatioPromptOverFONLL')
    hRatioPromptOverFONLL.Divide(hRatioPromptOverFONLL, hFONLLPromptCentral)
    hRatioPromptOverFONLL.GetYaxis().SetTitle('Data / FONLL')

    if args.syst:
        gRatioPromptOverFONLL = DivideGraphByHisto(gCorrYieldPrompt, hFONLLPromptCentral, False)
        SetObjectStyle(gRatioPromptOverFONLL, color=kBlack, fillstyle=0)

    lineFONLLPrompt = TLine(ptMinPrompt, 1., ptMaxPrompt, 1.)
    lineFONLLPrompt.SetLineColor(kRed+1)
    lineFONLLPrompt.SetLineWidth(2)

if args.FD:
    infileFD = TFile.Open(args.FD)
    hCorrYieldFD = infileFD.Get('hCorrYield')
    if args.syst:
        gCorrYieldFD = infileFD.Get('gCorrYieldSystTot')
    if not hCorrYieldFD:
        print("Not hCorrYieldFD")
        hCorrYieldFD = infileFD.Get('histoSigmaCorr')
        hCorrYieldFD.SetStats(0)
        hCorrYieldFD.Scale(1.e-6 / BR)
        if args.syst:
            gCorrYieldFD = infileFD.Get('gSigmaCorr')
            ScaleGraph(gCorrYieldFD, 1.e-6 / BR)
            for iPt in range(hCorrYieldFD.GetNbinsX()):
                gCorrYieldFD.SetPointEXhigh(iPt, hCorrYieldFD.GetBinWidth(iPt+1)/2.5)
                gCorrYieldFD.SetPointEXlow(iPt, hCorrYieldFD.GetBinWidth(iPt+1)/2.5)
    hCorrYieldFD.SetName('hCorrYieldFD')
    hCorrYieldFD.SetDirectory(0)
    SetObjectStyle(hCorrYieldFD, color=kBlack, markerstyle=kFullSquare)
    if args.syst:
        SetObjectStyle(gCorrYieldFD, color=kBlack, fillstyle=0)
        gCorrYieldFD.SetName('gCorrYieldFD')
    infileFD.Close()

    ptLimitsFD = np.array(hCorrYieldFD.GetXaxis().GetXbins(), 'd')
    ptMinFD = list(ptLimitsFD)[0]
    ptMaxFD = list(ptLimitsFD)[-1]
    if not args.prompt:
        ptMin = ptMinFD
        ptMax = ptMaxFD
        sigmaMin = hCorrYieldFD.GetBinContent(hCorrYieldFD.GetNbinsX())*0.2
        sigmaMax = hCorrYieldFD.GetBinContent(1)*5
    else:
        if ptMinFD < ptMin:
            ptMin = ptMinFD
        if ptMaxFD > ptMax:
            ptMax = ptMaxFD
        if hCorrYieldFD.GetBinContent(hCorrYieldFD.GetNbinsX())*0.2 < sigmaMin:
            sigmaMin = hCorrYieldFD.GetBinContent(hCorrYieldFD.GetNbinsX())*0.2
        if hCorrYieldFD.GetBinContent(1)*5 > sigmaMax:
            sigmaMax = hCorrYieldFD.GetBinContent(1)*5

    infileFONLL = TFile.Open(args.FONLLFileName)
    if args.Dplus:
        hFONLLFDCentral = infileFONLL.Get('hDpluskpipifromBpred_central_corr')
        hFONLLFDMin = infileFONLL.Get('hDpluskpipifromBpred_min_corr')
        hFONLLFDMax = infileFONLL.Get('hDpluskpipifromBpred_max_corr')
    elif args.Ds:
        hFONLLFDCentral = infileFONLL.Get('hDsPhipitoKkpifromBpred_central_corr')
        hFONLLFDMin = infileFONLL.Get('hDsPhipitoKkpifromBpred_min_corr')
        hFONLLFDMax = infileFONLL.Get('hDsPhipitoKkpifromBpred_max_corr')
    hFONLLFDCentral = hFONLLFDCentral.Rebin(hCorrYieldFD.GetNbinsX(), 'hFONLLFDCentral', ptLimitsFD)
    hFONLLFDMin = hFONLLFDMin.Rebin(hCorrYieldFD.GetNbinsX(), 'hFONLLFDMin', ptLimitsFD)
    hFONLLFDMax = hFONLLFDMax.Rebin(hCorrYieldFD.GetNbinsX(), 'hFONLLFDMax', ptLimitsFD)
    hFONLLFDCentral.Scale(1.e-6 / BR * ApPb /20, 'width')
    hFONLLFDMin.Scale(1.e-6 / BR * ApPb/20, 'width')
    hFONLLFDMax.Scale(1.e-6 / BR * ApPb/20, 'width')
    hFONLLFDCentral.SetDirectory(0)
    hFONLLFDMin.SetDirectory(0)
    hFONLLFDMax.SetDirectory(0)
    hFONLLFDCentral.SetStats(0)
    hFONLLFDMin.SetStats(0)
    hFONLLFDMax.SetStats(0)
    gFONLLFD = TGraphAsymmErrors(1)
    
    print(hFONLLFDCentral.GetNbinsX())

    for iPt in range(hFONLLFDCentral.GetNbinsX()):
        gFONLLFD.SetPoint(iPt, hFONLLFDCentral.GetBinCenter(iPt+1), hFONLLFDCentral.GetBinContent(iPt+1))
        gFONLLFD.SetPointError(iPt, hFONLLFDCentral.GetBinWidth(iPt+1)/2,
                               hFONLLFDCentral.GetBinWidth(iPt+1)/2,
                               hFONLLFDCentral.GetBinContent(iPt+1)-hFONLLFDMin.GetBinContent(iPt+1),
                               hFONLLFDMax.GetBinContent(iPt+1)-hFONLLFDCentral.GetBinContent(iPt+1))

    gFONLLFDUnc = ComputeRatioGraph(gFONLLFD, gFONLLFD, False)
    SetObjectStyle(hFONLLFDCentral, color=kAzure+4, markerstyle=0, fillstyle=0)
    SetObjectStyle(gFONLLFD, color=kAzure+4, fillalpha=0.2, linecolor=kAzure+4, fillstyle=1001)
    SetObjectStyle(gFONLLFDUnc, color=kAzure+4, fillalpha=0.2, fillstyle=1000)
    infileFONLL.Close()

    hRatioFDOverFONLL = hCorrYieldFD.Clone('hRatioFDOverFONLL')
    hRatioFDOverFONLL.Divide(hRatioFDOverFONLL, hFONLLFDCentral)
    hRatioFDOverFONLL.GetYaxis().SetTitle('Data / FONLL + Pythia8 (e^{+}e^{-} FF)')

    if args.syst:
        gRatioFDOverFONLL = DivideGraphByHisto(gCorrYieldFD, hFONLLFDCentral, False)
        SetObjectStyle(gRatioFDOverFONLL, color=kBlack, fillstyle=0)

    lineFONLLFD = TLine(ptMinFD, 1., ptMaxFD, 1.)
    lineFONLLFD.SetLineColor(kAzure+4)
    lineFONLLFD.SetLineWidth(2)

# protection for log scale
if sigmaMin <= 0:
    sigmaMin = 1.e-3
if args.logx and ptMin <= 0:
    print('WARNING: disabling log scale for x axis because minimum pT <= 0!')
    args.logx = False



cCrossSec = TCanvas('cCrossSec', '', 700, 800)
hFrame = cCrossSec.DrawFrame(ptMin, sigmaMin, ptMax, sigmaMax,
                             ';#it{p}_{T} (GeV/#it{c});d#sigma #times BR / d#it{p}_{T} (#mub GeV^{-1} #it{c})')
hFrame.GetXaxis().SetMoreLogLabels()
cCrossSec.SetLogy()
if args.logx:
    cCrossSec.SetLogx()

if args.logx:
    legPrompt = TLegend(0.65, 0.75, 0.85, 0.9)
    legFD = TLegend(0.25, 0.2, 0.45, 0.35)
else:
    legPrompt = TLegend(0.6, 0.8, 0.8, 0.95)
    legFD = TLegend(0.47, 0.62, 0.88, 0.8)
legPrompt.SetTextSize(0.04)
legPrompt.SetBorderSize(0)
legPrompt.SetFillStyle(0)
legFD.SetTextSize(0.045)
legFD.SetBorderSize(0)
legFD.SetFillStyle(0)





# def createCanvasPads():
#     c = TCanvas("c", "canvas", 800, 1000)
#     # Upper histogram plot is pad1
#     pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
#     pad1.SetBottomMargin(0)  # joins upper and lower plot
#     pad1.SetGridx()
#     pad1.SetLogy()
#     pad1.Range(ptMin, sigmaMin, ptMax, sigmaMax)
#     pad1.Draw()
#     # Lower ratio plot is pad2
#     c.cd()  # returns to main canvas before defining pad2
#     pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
#     pad2.SetTopMargin(0)  # joins upper and lower plot
#     pad2.SetBottomMargin(0.1)
#     pad2.SetGridx()
#     pad2.Draw()

#     return c, pad1, pad2

# doRatioPlot = False
# if doRatioPlot:
#     c, pad1, pad2 = createCanvasPads()
#     pad1.cd()
#     gFONLLFD.Draw('2')
#     hFONLLFDCentral.Draw('same')
#     hCorrYieldFD.Draw('esame')


#     pad2.cd()

#     hRatioFDOverFONLL.Draw('same')
#     print('hello')
# else:

if args.prompt:
    gFONLLPrompt.Draw('2')
    hFONLLPromptCentral.DrawCopy('same')
    hCorrYieldPrompt.DrawCopy('esame')
    if args.syst:
        gCorrYieldPrompt.Draw('2')
    legPrompt.AddEntry('', f'Prompt {mesonName}', '')
    legPrompt.AddEntry(hCorrYieldPrompt, 'Data', 'p')
    legPrompt.AddEntry(gFONLLPrompt, 'FONLL', 'f')
    legPrompt.Draw()
if args.FD:
    gFONLLFD.Draw('2')
    hFONLLFDCentral.DrawCopy('same')
    hCorrYieldFD.DrawCopy('esame')
    if args.syst:
        gCorrYieldFD.Draw('2')
    legFD.SetTextSize(0.04)
    legFD.AddEntry('', f'Non-prompt {mesonName}', '')
    legFD.AddEntry(hCorrYieldFD, 'Data', 'p')
    legFD.AddEntry(gFONLLFD, 'FONLL scaled by A_{Pb}', 'f')
    legFD.Draw()


lat = TLatex()
lat.SetNDC()
lat.SetTextSize(0.04)
lat.SetTextColor(kBlack)
lat.SetTextFont(42)

lat.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
lat.DrawLatex(0.22, 0.87, 'p-Pb, #sqrt{s} = 5.02 TeV       -0.96 < #it{y}_{cms} < 0.04')
lat.SetTextSize(0.02)
lat.DrawLatex(0.22, 0.25, '#pm 1.7% BR unc. not shown')
lat.DrawLatex(0.22, 0.2, '#pm 3.7% lumi. unc. not shown')
lat.SetTextSize(0.04)

cCrossSec.Update()
cCrossSec.SaveAs(args.outFileName)

if args.FD and args.prompt:
    cRatioToFONLL = TCanvas('cRatioToFONLL', '', 500, 1000)
    cRatioToFONLL.Divide(1, 2)
    hFramePrompt = cRatioToFONLL.cd(1).DrawFrame(ptMinPrompt, 0., ptMaxPrompt, 5.,
                                                ';#it{p}_{T} (GeV/#it{c});Data / FONLL')
    hFramePrompt.GetYaxis().SetDecimals()
    if args.logx:
        cRatioToFONLL.cd(1).SetLogx()
    gFONLLPromptUnc.Draw('2')
    lineFONLLPrompt.Draw('same')
    hRatioPromptOverFONLL.DrawCopy('same')
    if args.syst:
        gRatioPromptOverFONLL.Draw('2')
    lat.DrawLatex(0.7, 0.85, f'Prompt {mesonName}')
    hFrameFD = cRatioToFONLL.cd(2).DrawFrame(ptMinFD, 0., ptMaxFD, 3.,
                                            ';#it{p}_{T} (GeV/#it{c});Data / FONLL')
    hFrameFD.GetYaxis().SetDecimals()
    if args.logx:
        cRatioToFONLL.cd(2).SetLogx()
    gFONLLFDUnc.Draw('2')
    lineFONLLFD.Draw('same')
    hRatioFDOverFONLL.DrawCopy('same')
    if args.syst:
        gRatioFDOverFONLL.Draw('2')
    lat.DrawLatex(0.6, 0.85, f'Non-prompt {mesonName}')
else:
    cRatioToFONLL = TCanvas('cRatioToFONLL', '', 500, 500)
    if args.logx:
        cRatioToFONLL.SetLogx()
    if args.prompt:
        hFramePrompt = cRatioToFONLL.DrawFrame(ptMinPrompt, 0., ptMaxPrompt, 5.,
                                            ';#it{p}_{T} (GeV/#it{c});Data / FONLL')
        hFramePrompt.GetYaxis().SetDecimals()
        gFONLLPromptUnc.Draw('2')
        lineFONLLPrompt.Draw('same')
        hRatioPromptOverFONLL.DrawCopy('same')
        if args.syst:
            gRatioPromptOverFONLL.Draw('2')
        lat.DrawLatex(0.7, 0.85, f'Prompt {mesonName}')
    else:
        hFrameFD = cRatioToFONLL.DrawFrame(ptMinFD, 0., ptMaxFD, 3.,
                                        ';#it{p}_{T} (GeV/#it{c});Data / FONLL')
        hFrameFD.GetYaxis().SetDecimals()
        gFONLLFDUnc.Draw('2')
        lineFONLLFD.Draw('same')
        hRatioFDOverFONLL.DrawCopy('same')
        if args.syst:
            gRatioFDOverFONLL.Draw('2')
        lat.DrawLatex(0.6, 0.85, f'Non-prompt {mesonName}')

outFileRatioName = args.outFileName.replace('.pdf', '_RatioToFONLL.pdf')
cRatioToFONLL.Update()
cRatioToFONLL.SaveAs(outFileRatioName)

input('Press enter to exit')
