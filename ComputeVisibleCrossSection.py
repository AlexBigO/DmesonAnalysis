'''
python script to compute visible cross section
run: python ComputeVisibleCrossSection.py inFilePtDiffCrossSec.root outFile.root [--Dplus] [--Ds]
'''

import sys
import argparse
import ctypes
import numpy as np
from ROOT import TFile, TH1F, TGraphAsymmErrors, TCanvas, TLegend, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kFullCircle, kOpenCircle # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import ScaleGraph

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileNamePtDiffCrossSec', metavar='text', default='inFilePtDiffCrossSec.root',
                    help='root file with pT-differential cross section')
parser.add_argument('outFileName', metavar='text', default='outFile.root', help='output root file')
parser.add_argument('--Dplus', action='store_true', default=False, help='enable calculation for D+')
parser.add_argument('--Ds', action='store_true', default=False, help='enable calculation for Ds')
parser.add_argument('--Dzero', action='store_true', default=False, help='enable calculation for D0')
parser.add_argument('--Dstar', action='store_true', default=False, help='enable calculation for D*')
parser.add_argument('--PbPb', action='store_true', default=False, help='enable calculation for PbPb collisions')
parser.add_argument('--BR', action='store_true', default=False, help='divide cross section histogram by BR')
parser.add_argument('--ptmin', default=0., type=float, help='minimum pT')
parser.add_argument('--ptmax', default=1.e10, type=float, help='maximum pT')
args = parser.parse_args()

outFileNameVisibleCrossSecPDF = args.outFileName.replace('.root', '.pdf')

argChannels = [args.Dplus, args.Ds, args.Dzero, args.Dstar]
if not any(argChannels):
    print('ERROR: you should enable the comparison for either Dplus, Ds, Dzero, or Dstar! Exit')
    sys.exit()
elif len(np.where(argChannels)[0]) > 1:
    print('ERROR: you can enable only one meson at a time! Exit')
    sys.exit()

if not args.PbPb:
    lumiUnc = 0.037 # p-Pb ; 0.021 # pp@5.02 TeV
else:
    lumiUnc = 0.

if args.Dplus:
    BR = 0.0938
    BRunc = 0.0016
    mesonName = 'D^{+}'
elif args.Ds:
    BR = 0.0224
    BRunc = 0.0008
    mesonName = 'D_{s}^{+}'
elif args.Dzero:
    BR = 0.03950
    BRunc = 0.00031
    mesonName = 'D^{0}'
elif args.Dstar:
    BR = 0.03950*0.677
    BRunc = np.sqrt(0.00031**2/0.03950**2 + 0.005*2/0.677**2) * BR
    mesonName = 'D*^{+}'


SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14, titleoffsety=1.6) #, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.65, 0.1, 0.88, 0.35)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)


# load inputs
HFsystErr = []
inFileCross = TFile.Open(args.inFileNamePtDiffCrossSec)

if not args.BR:
    hCrossSection = inFileCross.Get('hCorrYield')
else: # need to divide by BR for pp data
    hCrossSection = inFileCross.Get('hCrossSection')
    hCrossSection.Scale(1./BR)

isDataDriven = True
if not hCrossSection:
    if not args.PbPb:
        hCrossSection = inFileCross.Get('histoSigmaCorr')
        gCrossSectionFDSyst = inFileCross.Get('gSigmaCorrConservative')
        gCrossSectionFDSyst.RemovePoint(0)
        hCrossSection.Scale(1.e-6 / BR)
        ScaleGraph(gCrossSectionFDSyst, 1.e-6 / BR)
    else:
        hCrossSection = inFileCross.Get('hAAC')
        gCrossSectionFDSyst = inFileCross.Get('gaaCsystFD')
    hCrossSection.SetStats(0)
    isDataDriven = False
hCrossSection.SetDirectory(0)
if not args.PbPb:
    if not args.Dzero:
        HFsystErr.append(inFileCross.Get('AliHFSystErr'))
    else:
        HFsystErr.append(inFileCross.Get('AliHFSystErrTopol'))
        HFsystErr.append(inFileCross.Get('AliHFSystErrLowpt'))
else:
    if not args.Dzero:
        HFsystErr.append(inFileCross.Get('AliHFSystErrAA'))
    else:
        if inFileCross.Get('AliHFSystErrAA1'):
            HFsystErr.append(inFileCross.Get('AliHFSystErrAA1'))
            HFsystErr.append(inFileCross.Get('AliHFSystErrAA2'))
        else:
            HFsystErr.append(inFileCross.Get('AliHFSystErrAA'))

inFileCross.Close()

# compute visible cross section
visCrossSec, statUnc, uncorrSystUnc, corrSystUncLow, corrSystUncHigh, \
    FDSystUncLow, FDSystUncHigh, trackSystUnc, systWoTrFDBRLumi = (0 for _ in range(9))

systRawYield, systFrac, systSelEff, systPID, systPtShape, systTr = (0. for _ in range(6))

for iPt in range(hCrossSection.GetNbinsX()):

    binWidth = hCrossSection.GetBinWidth(iPt+1)
    ptCent = hCrossSection.GetBinCenter(iPt+1)

    # integrate only between ptmin and ptmax
    if ptCent-binWidth/2 < args.ptmin:
        continue
    if ptCent+binWidth/2 > args.ptmax:
        continue

    visCrossSec += hCrossSection.GetBinContent(iPt+1) * binWidth
    statUnc += (hCrossSection.GetBinError(iPt+1) * binWidth)**2

    # get correct systematic uncertainties
    if ptCent < 1 and args.Dzero:
        systErr = HFsystErr[1] # Dzero low pT
    else:
        systErr = HFsystErr[0] # topological

    # uncorrelated systematic uncertainty (yield extraction)
    uncorrSystUnc += systErr.GetRawYieldErr(ptCent)**2 * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    systRawYield += systErr.GetRawYieldErr(ptCent)**2 * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    # correlated systematic uncertainty (sel eff, PID eff, gen pT shape, tracking, FD, BR)

    systFrac += systErr.GetDataDrivenFDErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth
    systSelEff += systErr.GetCutsEffErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth
    systPID += systErr.GetPIDEffErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth
    systPtShape += systErr.GetMCPtShapeErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth
    systTr += systErr.GetTrackingEffErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth

    # trackSystUnc += systErr.GetTrackingEffErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth
    # systWoTrFDBRLumi += (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
    #     systErr.GetPIDEffErr(ptCent)**2) * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    # totCorrSystUncSqLow = (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
    #     systErr.GetPIDEffErr(ptCent)**2 + systErr.GetTrackingEffErr(ptCent)**2) \
    #         * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    # totCorrSystUncSqHigh = (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
    #     systErr.GetPIDEffErr(ptCent)**2 + systErr.GetTrackingEffErr(ptCent)**2) \
    #         * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    # if isDataDriven:
    #     FDSystUncLow += systErr.GetDataDrivenFDErr(ptCent) * binWidth
    #     FDSystUncHigh += systErr.GetDataDrivenFDErr(ptCent) * binWidth
    #     totCorrSystUncSqLow += systErr.GetDataDrivenFDErr(ptCent)**2 \
    #         * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    #     totCorrSystUncSqHigh += systErr.GetDataDrivenFDErr(ptCent)**2 \
    #         * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    # else:
    #     for iPtFD in range(gCrossSectionFDSyst.GetN()):
    #         ptCentFD, sigma = ctypes.c_double(), ctypes.c_double()
    #         gCrossSectionFDSyst.GetPoint(iPtFD, ptCentFD, sigma)
    #         if abs(ptCentFD.value-ptCent) < 0.01:
    #             break
    #     FDSystUncLow += gCrossSectionFDSyst.GetErrorYlow(iPtFD) * binWidth
    #     FDSystUncHigh += gCrossSectionFDSyst.GetErrorYlow(iPtFD) * binWidth
    #     totCorrSystUncSqLow += gCrossSectionFDSyst.GetErrorYlow(iPtFD)**2 * binWidth**2
    #     totCorrSystUncSqHigh += gCrossSectionFDSyst.GetErrorYhigh(iPtFD)**2 * binWidth**2

    # corrSystUncLow += np.sqrt(totCorrSystUncSqLow)
    # corrSystUncHigh += np.sqrt(totCorrSystUncSqHigh)

statUnc = np.sqrt(statUnc)
uncorrSystUnc = np.sqrt(uncorrSystUnc)

corrSystUnc = np.sqrt(systSelEff**2 + systPID**2 + systPtShape**2 + systTr**2 + systFrac**2)
corrSystUncLow = corrSystUnc
corrSystUncHigh = corrSystUnc

corrSystUncLowNoBRAndLumi = corrSystUncLow
corrSystUncHighNoBRAndLumi = corrSystUncHigh

corrSystUncLow = np.sqrt(corrSystUncLow**2 + (lumiUnc**2 + BRunc**2/BR**2)*visCrossSec)
corrSystUncHigh = np.sqrt(corrSystUncHigh**2 + (lumiUnc**2 + BRunc**2/BR**2)*visCrossSec)

totSystUncLowWoBRAndLumi = np.sqrt(uncorrSystUnc**2 + corrSystUncLowNoBRAndLumi**2)
totSystUncHighWoBRAndLumi = np.sqrt(uncorrSystUnc**2 + corrSystUncHighNoBRAndLumi**2)

totSystUncLow = np.sqrt(uncorrSystUnc**2 + corrSystUncLow**2)
totSystUncHigh = np.sqrt(uncorrSystUnc**2 + corrSystUncHigh**2)

# fill histos and graphs
if args.ptmin <= hCrossSection.GetBinLowEdge(1):
    ptMin = hCrossSection.GetBinLowEdge(1)
else:
    ptMin = args.ptmin
if args.ptmax >= hCrossSection.GetXaxis().GetBinUpEdge(hCrossSection.GetNbinsX()):
    ptMax = hCrossSection.GetXaxis().GetBinUpEdge(hCrossSection.GetNbinsX())
else:
    ptMax = args.ptmax

graphName = 'VisCrossSec'
utilStr = f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})'
if args.PbPb:
    graphName = 'PtVisYield'
    utilStr = f';;d#it{{N}}/d#it{{y}} ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})'

hVisibleCrossSectionStat = TH1F(f'h{graphName}', utilStr, 1, -0.96, 0.04)

gVisCrossSecUncorrSyst, gVisCrossSecCorrSyst, gVisCrossSecTotSyst, \
    gVisCrossSecSystWoBRAndLumi, gVisCrossSecSystWoTrFDBRAndLumi, \
    gVisCrossSecSystTracking, gVisCrossSecSystFD, gVisCrossSecSystLumi, gVisCrossSecSystBR = \
        (TGraphAsymmErrors(1) for _ in range(9))
gVisCrossSecUncorrSyst.SetNameTitle(f'g{graphName}UncorrSys', utilStr)
gVisCrossSecCorrSyst.SetNameTitle(f'g{graphName}PtCorrSys', utilStr)
gVisCrossSecTotSyst.SetNameTitle(f'g{graphName}TotSys', utilStr)
gVisCrossSecSystWoBRAndLumi.SetNameTitle(f'g{graphName}SysWoBRAndLumi', utilStr)
gVisCrossSecSystWoTrFDBRAndLumi.SetNameTitle(f'g{graphName}SysWoTrFDBRAndLumi', utilStr)
gVisCrossSecSystTracking.SetNameTitle('gVisCrossSecSysTracking', utilStr)
gVisCrossSecSystFD.SetNameTitle(f'g{graphName}SysFD', utilStr)
gVisCrossSecSystLumi.SetNameTitle(f'g{graphName}SysLumi', utilStr)
gVisCrossSecSystBR.SetNameTitle(f'g{graphName}SysBR', utilStr)

SetObjectStyle(hVisibleCrossSectionStat, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(gVisCrossSecUncorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecCorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecTotSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystLumi, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystBR, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystFD, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystTracking, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystWoBRAndLumi, linecolor=kRed+2, markercolor=kRed+2, fillstyle=0)
SetObjectStyle(gVisCrossSecSystWoTrFDBRAndLumi, color=kBlack, fillstyle=0)

hVisibleCrossSectionStat.SetBinContent(1, visCrossSec)
hVisibleCrossSectionStat.SetBinError(1, statUnc)
gVisCrossSecUncorrSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecCorrSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecTotSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystWoBRAndLumi.SetPoint(0, 0.04-0.5, visCrossSec)
gVisCrossSecSystWoTrFDBRAndLumi.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystTracking.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystFD.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystLumi.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystBR.SetPoint(0, 1., visCrossSec)
gVisCrossSecUncorrSyst.SetPointError(0, 0.3, 0.3, uncorrSystUnc, uncorrSystUnc)
gVisCrossSecCorrSyst.SetPointError(0, 0.3, 0.3, corrSystUncLow, corrSystUncHigh)
gVisCrossSecTotSyst.SetPointError(0, 0.3, 0.3, totSystUncLow, totSystUncHigh)
gVisCrossSecSystWoBRAndLumi.SetPointError(0, 0.3, 0.3, totSystUncLowWoBRAndLumi, totSystUncHighWoBRAndLumi)
gVisCrossSecSystWoTrFDBRAndLumi.SetPointError(0, 0.3, 0.3, systWoTrFDBRLumi, systWoTrFDBRLumi)
gVisCrossSecSystTracking.SetPointError(0, 0.3, 0.3, trackSystUnc, trackSystUnc)
gVisCrossSecSystFD.SetPointError(0, 0.3, 0.3, FDSystUncLow, FDSystUncHigh)
gVisCrossSecSystLumi.SetPointError(0, 0.3, 0.3, visCrossSec*lumiUnc, visCrossSec*lumiUnc)
gVisCrossSecSystBR.SetPointError(0, 0.3, 0.3, visCrossSec*BRunc/BR, visCrossSec*BRunc/BR)

print(f'Visible cross section: {visCrossSec:.1f} +- {statUnc:.1f} (stat) \
      +- {totSystUncLowWoBRAndLumi:.1f} (syst) +- {visCrossSec*lumiUnc:.1f} (lumi) \
      {visCrossSec*BRunc/BR:.1f} (BR)')


# Canvas
axisTitleCrossSection = '; y; d#sigma/d#it{y} (#mub)'

cVisibleCrossSection = TCanvas('cVisibleCrossSection', '', 800, 800)
cVisibleCrossSection.DrawFrame(-0.96, 0.8*visCrossSec, 0.04, 1.2*visCrossSec, axisTitleCrossSection)
SetObjectStyle(hVisibleCrossSectionStat, linecolor=kRed+2, markercolor=kRed+2,
                   markerstyle=kOpenCircle)
hVisibleCrossSectionStat.SetDirectory(0)
gVisCrossSecSystWoBRAndLumi.Draw('2')
hVisibleCrossSectionStat.Draw('same')

latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13);  # align at top
latex.SetTextFont(42)
latex.DrawLatex(0.22, 0.92, 'WORK IN PROGRESS')
latex.DrawLatex(0.22, 0.87, f'p-Pb, #sqrt{{s_{{NN}}}} = 5.02 TeV       {ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}}')
latex.SetTextSize(0.025)
latex.DrawLatex(0.22, 0.25, '#pm 1.7% BR unc. not shown')
latex.DrawLatex(0.22, 0.2, '#pm 3.7% lumi. unc. not shown')

leg.AddEntry(hVisibleCrossSectionStat, 'Non-prompt D^{+}', 'p')
leg.Draw()

# otput file
outFile = TFile.Open(args.outFileName, 'recreate')
hVisibleCrossSectionStat.Write()
gVisCrossSecUncorrSyst.Write()
gVisCrossSecCorrSyst.Write()
gVisCrossSecTotSyst.Write()
gVisCrossSecSystWoBRAndLumi.Write()
gVisCrossSecSystWoTrFDBRAndLumi.Write()
gVisCrossSecSystTracking.Write()
gVisCrossSecSystFD.Write()
gVisCrossSecSystBR.Write()
if not args.PbPb:
    gVisCrossSecSystLumi.Write()
cVisibleCrossSection.Write()
outFile.Close()

cVisibleCrossSection.SaveAs(outFileNameVisibleCrossSecPDF)

print(f'Saved output file: {args.outFileName}')
