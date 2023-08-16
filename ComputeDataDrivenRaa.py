'''
python script for the computation of the nuclear modification factor of prompt or non-prompt D mesons
run: python ComputeDataDrivenRaa.py corrYieldFile.root ppCrossSecFile.root outFile.root
                                    [--centrality] [--batch]
'''

import sys
import argparse
import numpy as np

from ROOT import gROOT, TFile, TGraphErrors, TCanvas, TLine, TLegend, TLatex # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor
from utils.AnalysisUtils import ComputeRatioDiffBins

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('corrYieldFileName', metavar='text', default='corrYieldFile.root',
                    help='root file with corrected yields')
parser.add_argument('ppCrossSecFileName', metavar='text', default='ppCrossSecFile.root',
                    help='root file with pp reference cross section')
parser.add_argument('outFileName', metavar='text', default='outFile.root', help='root output file name')
parser.add_argument('--energy', metavar='text', default='5.02', help='energy (5.02)')
parser.add_argument('--centrality', metavar='text', default='0100', help='centrality class (0100, 010, 3050, 6080)')
parser.add_argument("--batch", action='store_true', help='suppress video output', default=False)
args = parser.parse_args()

if args.energy != '5.02':
    print('ERROR: Only 5.02 TeV energy implemented! Exit')
    sys.exit()

#values from https://alice-notes.web.cern.ch/system/files/notes/public/711/2019-04-02-ALICE_public_note.pdf
if args.centrality == '010':
    taa = 23.26 # 1 / mbarn
    taaUnc = 0.168
elif args.centrality == '3050':
    taa = 3.917
    taaUnc = 0.0645
elif args.centrality == '6080':
    taa = 0.4188
    taaUnc = 0.0106
elif args.centrality == '0100':
    taa = 208
    taaUnc = 0.
    # taa = 0.09923
    # taaUnc = 0.0017
else:
    print('ERROR: Only 0-10, 30-50 and 60-80 centrality classes implemented! Exit')
    sys.exit()

# load input file
corrYieldFile = TFile.Open(args.corrYieldFileName)
hCorrYieldPbPb =  corrYieldFile.Get('hCorrYield')
hCorrYieldPbPb.SetName('hCorrYieldPbPb')
gCorrYieldPbPbSystTot = corrYieldFile.Get('gCorrYieldSystTot')
gCorrYieldPbPbSystTot.SetName('gCorrYieldSystTotpPb')
systErrPbPb = corrYieldFile.Get('AliHFSystErr')
systErrPbPb.SetName('AliHFSystErrPbPb')

ppCrossSecFile = TFile.Open(args.ppCrossSecFileName)
hCrossSectionPP = ppCrossSecFile.Get('hCrossSection')
hCrossSectionPP.SetName('hCrossSectionPP')
gCrossSectionPPSystTot = ppCrossSecFile.Get('gCrossSectionSystTot')
gCrossSectionPPSystTot.SetName('gCrossSectionSystTotPP')
gCrossSectionPPSystLumi = ppCrossSecFile.Get('gCrossSectionSystLumi')
systErrPP = ppCrossSecFile.Get('AliHFSystErr')
systErrPP.SetName('AliHFSystErrPP')

rapidityCorrectionFile = TFile.Open("/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/6_RpPb/rapidity_correction/rapidity_correction_rebinned_2_16.root")
hRapidityCorrection = rapidityCorrectionFile.Get('hYCorrection_central')

# TODO: improve protection checking the limits of the bins besides the number
# if hCorrYieldPbPb.GetNbinsX() != hCrossSectionPP.GetNbinsX():
#     print('ERROR: inconsistent number of bins in input objects! Exit')
#     sys.exit()

# hRaa = hCorrYieldPbPb.Clone('hRaa')
# hRaa.Divide(hCrossSectionPP)

# compute RpPb
hRaa = ComputeRatioDiffBins(hCorrYieldPbPb, hCrossSectionPP)

# apply rapidity correction
for iPt in range(hRaa.GetNbinsX()):
    raa = hRaa.GetBinContent(iPt+1)
    raa_statunc = hRaa.GetBinError(iPt+1)
    y_correction = hRapidityCorrection.GetBinContent(iPt+1)
    hRaa.SetBinContent(iPt+1, raa / y_correction)
    hRaa.SetBinError(iPt+1, raa_statunc / y_correction)

# hRaa = ComputeRatioDiffBins(hRaa, hRapidityCorrection)
# hRaa.Divide(hRapidityCorrection)

# name settings
hRaa.SetName('hRpPb')
hRaa.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{R}_{pPb}')
# hRaa.Scale(1.e3 / taa) #convert Taa to 1 / ub
hRaa.Scale(1. / taa)

gRaaSystTot = TGraphErrors(1)
gRaaSystTot.SetName('gRaaSystTot')
gRaaSystTot.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{R}_{AA}')
SetObjectStyle(gRaaSystTot, color=GetROOTColor('kBlack'), fillstyle=0)
gRaaSystTaa = TGraphErrors(1)
gRaaSystTaa.SetName('gRaaSystTaa')
gRaaSystTaa.SetTitle('Taa syst. unc.;;')
gRaaSystTaa.SetPoint(0, 1., 1.)
gRaaSystTaa.SetPointError(0, 0.4, taaUnc / taa)
SetObjectStyle(gRaaSystTaa, color=GetROOTColor('kBlue'), fillstyle=0)
gRaaSystNorm = TGraphErrors(1)
gRaaSystNorm.SetName('gRaaSystNorm')
gRaaSystNorm.SetTitle('Normalization syst. unc. (pp norm. + Taa);;')
gRaaSystNorm.SetPoint(0, 1., 1.)
gRaaSystNorm.SetPointError(0, 0.4, np.sqrt((taaUnc / taa)**2 + gCrossSectionPPSystLumi.GetErrorY(0)**2))
SetObjectStyle(gRaaSystNorm, color=GetROOTColor('kAzure+4'), fillstyle=0)

for iPt in range(hRaa.GetNbinsX()):
    ppValue = gCrossSectionPPSystTot.GetPointY(iPt)
    ppSystUnc = gCrossSectionPPSystTot.GetErrorY(iPt)
    pbpbValue = gCorrYieldPbPbSystTot.GetPointY(iPt)
    pbpbSystUnc = gCorrYieldPbPbSystTot.GetErrorY(iPt)

    raa = hRaa.GetBinContent(iPt+1)
    relRaaSystUnc = np.sqrt((ppSystUnc / ppValue)**2 + (pbpbSystUnc / pbpbValue)**2)
    ptCent = hRaa.GetBinCenter(iPt+1)

    gRaaSystTot.SetPoint(iPt, ptCent, raa)
    gRaaSystTot.SetPointError(iPt, 0.4, relRaaSystUnc * raa)

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14)

cRaa = TCanvas('cRaa', '', 700, 800)
hRaa.GetYaxis().SetRangeUser(0., 2.7)
hRaa.Draw('same')
gRaaSystTot.Draw('2')
line = TLine(hRaa.GetXaxis().GetXmin(), 1., hRaa.GetXaxis().GetXmax(), 1.)
line.SetLineColor(GetROOTColor('kBlack'))
line.SetLineStyle(9)
line.SetLineWidth(2)
line.Draw("same")

# description = TLatex()
# description.SetTextSize(0.045)
# description.DrawLatex(0.3, 0.9, "ALICE")

leg = TLegend(0.60, 0.1, 0.84, 0.3)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.045)
leg.AddEntry(hRaa, 'Non-prompt D^{+}', 'l')
leg.Draw()

cRaa.Update()

outFile = TFile(args.outFileName, 'recreate')
hRaa.Write()
gRaaSystTot.Write()
gRaaSystNorm.Write()
gRaaSystTaa.Write()
hCorrYieldPbPb.Write()
gCorrYieldPbPbSystTot.Write()
hCrossSectionPP.Write()
gCrossSectionPPSystTot.Write()
cRaa.Write()
systErrPbPb.Write()
systErrPP.Write()
outFile.Close()

if not args.batch:
    input('Press enter to exit')
