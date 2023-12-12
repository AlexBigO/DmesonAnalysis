'''
python script for the FONLL extrapolation of the pT-integrated cross sections of non-prompt D mesons
'''

import sys
import argparse
import numpy as np
import yaml
from ROOT import TFile, TH1D, TGraphAsymmErrors, TCanvas # pylint: disable=import-error,no-name-in-module
from ROOT import kGray, kBlack, kRed # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.ReadModel import ReadFONLL
from utils.AnalysisUtils_Fabrizio import ComputebbbarExtrapolationFactor, ComputePtExtrapolationFactor, CombineExtrapFactors

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('--Dplus', action='store_true', default=False, help='enable calculation for D+')
parser.add_argument('--Ds', action='store_true', default=False, help='enable calculation for Ds')
parser.add_argument('--Dzero', action='store_true', default=False, help='enable calculation for D0')
args = parser.parse_args()

if not args.Dplus and not args.Ds and not args.Dzero:
    print('ERROR: you should enable the comparison for either D+, Ds, or D0! Exit')
    sys.exit()
elif (args.Dplus and args.Ds) or (args.Dplus and args.Dzero) or (args.Dzero and args.Ds):
    print('ERROR: you cannot enable the comparison for more than one meson! Exit')
    sys.exit()

with open('/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/Extrapolation/config_input_files.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

if args.Dzero:
    inFileName = inputCfg['nonprompt']['Dzero']['crosssection']
    mesonName = 'Dzero'
elif args.Dplus:
    inFileName = inputCfg['nonprompt']['Dplus']['crosssection']
    mesonName = 'Dplus'
else:
    inFileName = inputCfg['nonprompt']['Ds']['crosssection']
    mesonName = 'Ds'

BRD = inputCfg['BR'][mesonName]
BRDunc = inputCfg['BRunc'][mesonName]
lumiRelUnc = inputCfg['lumiunc']
suffix = inputCfg['nonprompt']['Dplus']['suffix']

# Data
inFile = TFile.Open(inFileName)
if args.Dzero:
    hCrossSection = inFile.Get('hXSec')
    hCrossSection.SetDirectory(0)
    # from pb to mub to be consistend with D+ and Ds
    hCrossSection.Scale(1.e-6) # BR already included in original file
    hNonPromptD0SystRawYield = inFile.Get('hSys_MassFit')
    hNonPromptD0SystSelEff = inFile.Get('hSys_Selection')
    hNonPromptD0SystTrEff = inFile.Get('hSys_Tracking')
    hNonPromptD0SystPIDEff = inFile.Get('hSys_PID')
    hNonPromptD0SystPtShape = inFile.Get('hSys_pTShape')
    hNonPromptD0SystFrac = inFile.Get('hSys_Fraction')
    hNonPromptD0SystRawYield.SetDirectory(0)
    hNonPromptD0SystSelEff.SetDirectory(0)
    hNonPromptD0SystTrEff.SetDirectory(0)
    hNonPromptD0SystPIDEff.SetDirectory(0)
    hNonPromptD0SystPtShape.SetDirectory(0)
    hNonPromptD0SystFrac.SetDirectory(0)
else:
    if 'pp' in suffix:
        hCrossSection = inFile.Get('hCrossSection')
    else:
        hCrossSection = inFile.Get('hCorrYield')
    hCrossSection.Scale(1. / BRD)
    hCrossSection.SetDirectory(0)
    systErr = inFile.Get('AliHFSystErr')
inFile.Close()

nPtBins = hCrossSection.GetNbinsX()
ptMin = hCrossSection.GetBinLowEdge(1)
ptMax = hCrossSection.GetXaxis().GetBinUpEdge(nPtBins)

# compute visible cross section
visCrossSec, visCrossSecStat, systRawYield, systFrac, systSelEff, systPID, systPtShape, systTr = (0. for _ in range(8))
for iPt in range(1, nPtBins+1):
    ptCent = hCrossSection.GetBinCenter(iPt)
    ptWidth = hCrossSection.GetBinWidth(iPt)
    yCent = hCrossSection.GetBinContent(iPt)
    statUnc = hCrossSection.GetBinError(iPt)
    visCrossSec += yCent * ptWidth
    visCrossSecStat += statUnc**2 * ptWidth**2
    if args.Dzero:
        # uncorr vs pT
        systRawYield += hNonPromptD0SystRawYield.GetBinContent(iPt)**2 * yCent**2 * ptWidth**2
        systFrac += hNonPromptD0SystFrac.GetBinContent(iPt)**2 * yCent**2 * ptWidth**2
        # corr vs pT
        systSelEff += hNonPromptD0SystSelEff.GetBinContent(iPt) * yCent * ptWidth
        systPID += hNonPromptD0SystPIDEff.GetBinContent(iPt) * yCent * ptWidth
        systPtShape += hNonPromptD0SystPtShape.GetBinContent(iPt) * yCent * ptWidth
        systTr += hNonPromptD0SystTrEff.GetBinContent(iPt) * yCent * ptWidth
    else:
        # uncorr vs pT
        systRawYield += systErr.GetRawYieldErr(ptCent)**2 * yCent**2 * ptWidth**2
        # systFrac += systErr.GetDataDrivenFDErr(ptCent)**2 * yCent**2 * ptWidth**2
        # corr vs pT
        systFrac += systErr.GetDataDrivenFDErr(ptCent) * yCent * ptWidth
        systSelEff += systErr.GetCutsEffErr(ptCent) * yCent * ptWidth
        systPID += systErr.GetPIDEffErr(ptCent) * yCent * ptWidth
        systPtShape += systErr.GetMCPtShapeErr(ptCent) * yCent * ptWidth
        systTr += systErr.GetTrackingEffErr(ptCent) * yCent * ptWidth

visCrossSecStat = np.sqrt(visCrossSecStat)
systRawYield = np.sqrt(systRawYield)
# systFrac = np.sqrt(systFrac)

visCrossSecBRSyst = BRDunc/BRD * visCrossSec
visCrossSecLumiSyst = lumiRelUnc * visCrossSec
visCrossSecFDSyst = systFrac
visCrossSecTrSyst = systTr
visCrossSecDataSystWoTrFDBRLumi = np.sqrt(systRawYield**2 + systSelEff**2 + systPID**2 + systPtShape**2)
visCrossSecUncorrSyst = np.sqrt(systRawYield**2) # + systFrac**2)
visCrossSecCorrSyst = np.sqrt(systSelEff**2 + systPID**2 + systPtShape **
                              2 + systTr**2 + visCrossSecLumiSyst**2 + visCrossSecBRSyst**2 + systFrac**2)
visCrossSecDataSyst = np.sqrt(visCrossSecUncorrSyst**2 + visCrossSecCorrSyst**2)

# FONLL
if args.Dzero:
    histoNameFONLL = 'hD0Kpi'
elif args.Dplus:
    histoNameFONLL = 'hDpluskpipi'
else:
    histoNameFONLL = 'hDsPhipitoKkpi'

hFONLLPtExtrap, hFONLLbbbarExtrap = ({} for _ in range(2))
inFileFONLLPtExtrap = TFile.Open(inputCfg['extrap']['ptint']['default'])
# inFileFONLLbbbarExtrap = TFile.Open(inputCfg['extrap']['bbbar']['denominator']['default'])
for pred in ['central', 'min', 'max']:
    hFONLLPtExtrap[pred] = inFileFONLLPtExtrap.Get(f'{histoNameFONLL}fromBpred_{pred}_corr')
    # hFONLLbbbarExtrap[pred] = inFileFONLLbbbarExtrap.Get(f'{histoNameFONLL}fromBpred_{pred}_corr')

hFONLLPtExtrapAlt, hFONLLbbbarExtrapAlt = ([] for _ in range(2))
for iFONLL, fileName in enumerate(inputCfg['extrap']['ptint']['alternative']):
    hFONLLPtExtrapAlt.append({})
    inFileFONLLPtExtrapAlt = TFile.Open(fileName)
    for pred in ['central', 'min', 'max']:
        hFONLLPtExtrapAlt[iFONLL][pred] = inFileFONLLPtExtrapAlt.Get(f'{histoNameFONLL}fromBpred_{pred}_corr')
        hFONLLPtExtrapAlt[iFONLL][pred].SetName(f'{hFONLLPtExtrapAlt[iFONLL][pred].GetName()}PtAlt{iFONLL}')
        hFONLLPtExtrapAlt[iFONLL][pred].SetDirectory(0)
    inFileFONLLPtExtrapAlt.Close()
# for iFONLL, fileName in enumerate(inputCfg['extrap']['bbbar']['denominator']['alternative']):
#     hFONLLbbbarExtrapAlt.append({})
#     inFileFONLLbbbarExtrapAlt = TFile.Open(fileName)
#     for pred in ['central', 'min', 'max']:
#         hFONLLbbbarExtrapAlt[iFONLL][pred] = inFileFONLLbbbarExtrapAlt.Get(f'{histoNameFONLL}fromBpred_{pred}_corr')
#         hFONLLbbbarExtrapAlt[iFONLL][pred].SetName(f'{hFONLLbbbarExtrapAlt[iFONLL][pred].GetName()}bbbarAlt{iFONLL}')
#         hFONLLbbbarExtrapAlt[iFONLL][pred].SetDirectory(0)
#     inFileFONLLbbbarExtrapAlt.Close()

# for bbbar extrapolation
# _, FONLLbbbar = ReadFONLL(inputCfg['extrap']['bbbar']['numerator'])

# compute pT extrapolation factor
extrFactor = ComputePtExtrapolationFactor(hFONLLPtExtrap, ptMin, ptMax)
print('\n\x1b[32mpt extrapolation factor:\033[0m')
print(f'\tdefault       = {extrFactor["central"]:.4f} + {extrFactor["uncHigh"]:.4f}'
      f' - {extrFactor["uncLow"]:.4f}\033[0m')
extrFactorAlt = []
for iFONLL, _ in enumerate(hFONLLPtExtrapAlt):
    extrFactorAlt.append(ComputePtExtrapolationFactor(hFONLLPtExtrapAlt[iFONLL], ptMin, ptMax))
    print(f'\talternative {iFONLL+1} = {extrFactorAlt[iFONLL]["central"]:.4f} '
          f'+ {extrFactorAlt[iFONLL]["uncHigh"]:.4f} - {extrFactorAlt[iFONLL]["uncLow"]:.4f}')
# combine default and alternatives
print(f'\n\tchosen strategy for combination: \x1b[32m{inputCfg["extrap"]["ptint"]["combination"]}\033[0m')
extrFactorFin = CombineExtrapFactors(extrFactor, extrFactorAlt, inputCfg['extrap']['ptint']['combination'])
print(f'\t\x1b[32mfinal factor  = {extrFactorFin["central"]:.4f} + {extrFactorFin["uncHigh"]:.4f}'
      f' - {extrFactorFin["uncLow"]:.4f}\033[0m\n')

# compute bbbar extrapolation factor
# extrFactorbbbar = ComputebbbarExtrapolationFactor(FONLLbbbar, hFONLLbbbarExtrap, ptMin, ptMax, BRD)
# print('\n\x1b[32mbbbar extrapolation factor:\033[0m')
# print(f'\tdefault       = {extrFactorbbbar["central"]:.4f} + {extrFactorbbbar["uncHigh"]:.4f}'
#       f' - {extrFactorbbbar["uncLow"]:.4f}')
# extrFactorbbbarAlt = []
# for iFONLL, _ in enumerate(hFONLLbbbarExtrapAlt):
#     extrFactorbbbarAlt.append(
#         ComputebbbarExtrapolationFactor(FONLLbbbar, hFONLLbbbarExtrapAlt[iFONLL], ptMin, ptMax, BRD))
#     print(f'\talternative {iFONLL+1} = {extrFactorbbbarAlt[iFONLL]["central"]:.4f} '
#           f'+ {extrFactorbbbarAlt[iFONLL]["uncHigh"]:.4f} - {extrFactorbbbarAlt[iFONLL]["uncLow"]:.4f}\033[0m')
# # combine default and alternatives
# print(f'\n\tchosen strategy for combination: \x1b[32m{inputCfg["extrap"]["bbbar"]["combination"]}\033[0m')
# extrFactorbbbarFin = CombineExtrapFactors(extrFactorbbbar, extrFactorbbbarAlt,
#                                           inputCfg['extrap']['bbbar']['combination'])
# print(f'\t\x1b[32mfinal factor  = {extrFactorbbbarFin["central"]:.4f} + {extrFactorbbbarFin["uncHigh"]:.4f}'
#       f' - {extrFactorbbbarFin["uncLow"]:.4f}\033[0m\n')

# compute pT-integrated cross section
totCrossSec = visCrossSec * extrFactor['central']
totCrossSecStat = visCrossSecStat * extrFactor['central']
totCrossSecDataSyst = visCrossSecDataSyst * extrFactor['central']
totCrossSecCorrSyst = visCrossSecCorrSyst * extrFactor['central']
totCrossSecUncorrSyst = visCrossSecUncorrSyst * extrFactor['central']
totCrossSecDataSystWoTrFDBRLumi = visCrossSecDataSystWoTrFDBRLumi * extrFactor['central']
totCrossSecTrSyst = visCrossSecTrSyst * extrFactor['central']
totCrossSecFDSyst = visCrossSecFDSyst * extrFactor['central']
totCrossSecBRSyst = visCrossSecBRSyst * extrFactor['central']
totCrossSecLumiSyst = visCrossSecLumiSyst * extrFactor['central']
totCrossSecExtrapSystLow = visCrossSec * extrFactor['uncLow']
totCrossSecExtrapSystHigh = visCrossSec * extrFactor['uncHigh']

# # extrapolate to dsigma_bbbar / dy
# # load correction factors from Pythia and POWHEG for rapidity
# inFileRapidityPythia = TFile.Open(inputCfg['extrap']['bbbar']['rapiditycorr']['pythia']['central'])
# hCorrFactorRapidityPythia = inFileRapidityPythia.Get('hCorrFactorRapidityPythia')
# corrFactorYUncPythia = inputCfg['extrap']['bbbar']['rapiditycorr']['pythia']['uncertainty']
# if inputCfg['extrap']['bbbar']['rapiditycorr']['powheg']['apply']:
#     inFileRapidityPowheg = TFile.Open(inputCfg['extrap']['bbbar']['rapiditycorr']['powheg']['central'])
#     hCorrFactorRapidityPowheg = inFileRapidityPowheg.Get('hCorrFactorRapidityPOWHEG')
#     corrFactorYUncPowheg = inputCfg['extrap']['bbbar']['rapiditycorr']['powheg']['uncertainty']
#     corrFactorY = hCorrFactorRapidityPythia.GetBinContent(1) * hCorrFactorRapidityPowheg.GetBinContent(1)
#     corrFactorYUnc = np.sqrt(corrFactorYUncPythia**2 + corrFactorYUncPowheg**2)
# else: # do not apply powheg correction for J/psi crew
#     corrFactorY = hCorrFactorRapidityPythia.GetBinContent(1)
#     corrFactorYUnc = corrFactorYUncPythia

# totbbbarCrossSec = visCrossSec * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecStat = visCrossSecStat * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecDataSyst = visCrossSecDataSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecCorrSyst = visCrossSecCorrSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecUncorrSyst = visCrossSecUncorrSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecDataSystWoTrFDBRLumi = visCrossSecDataSystWoTrFDBRLumi * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecTrSyst = visCrossSecTrSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecFDSyst = visCrossSecFDSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecBRSyst = visCrossSecBRSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecLumiSyst = visCrossSecLumiSyst * extrFactorbbbarFin['central'] * corrFactorY
# totbbbarCrossSecExtrapSystLow = visCrossSec * extrFactorbbbarFin['uncLow'] * corrFactorY
# totbbbarCrossSecExtrapSystHigh = visCrossSec * extrFactorbbbarFin['uncHigh'] * corrFactorY
# totbbbarCrossSecExtrapRapSyst = totbbbarCrossSec * corrFactorYUnc

# draw results
SetGlobalStyle()

hPtIntCrossSecStat = TH1D('hPtIntCrossSecStat', ';;d#sigma/d#it{y} (#mub)', 1, 0.5, 1.5)
gPtIntCrossSecDataSyst, gPtIntCrossSecDataCorrSyst, gPtIntCrossSecDataUncorrSyst, \
    gPtIntCrossSecDataSystWoTrFDBRLumi, gPtIntCrossSecTrackingSyst, gPtIntCrossSecFDSyst, gPtIntCrossSecBRSyst, \
        gPtIntCrossSecLumiSyst, gPtIntCrossSecExtrapSyst = (TGraphAsymmErrors(1) for _ in range(9))
gPtIntCrossSecDataSyst.SetNameTitle('gPtIntCrossSecDataSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecDataCorrSyst.SetNameTitle('gPtIntCrossSecDataCorrSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecDataUncorrSyst.SetNameTitle('gPtIntCrossSecDataUncorrSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecDataSystWoTrFDBRLumi.SetNameTitle('gPtIntCrossSecDataSystWoTrFDBRLumi', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecTrackingSyst.SetNameTitle('gPtIntCrossSecTrackingSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecFDSyst.SetNameTitle('gPtIntCrossSecFDSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecBRSyst.SetNameTitle('gPtIntCrossSecBRSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecLumiSyst.SetNameTitle('gPtIntCrossSecLumiSyst', ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecExtrapSyst.SetNameTitle('gPtIntCrossSecExtrapSyst', ';;d#sigma/d#it{y} (#mub)')
SetObjectStyle(hPtIntCrossSecStat, color=kRed+2)
SetObjectStyle(gPtIntCrossSecDataSyst, color=kRed+2, fillstyle=0)
SetObjectStyle(gPtIntCrossSecDataCorrSyst, fillstyle=0)
SetObjectStyle(gPtIntCrossSecDataUncorrSyst, fillstyle=0)
SetObjectStyle(gPtIntCrossSecDataSystWoTrFDBRLumi, fillstyle=0)
SetObjectStyle(gPtIntCrossSecTrackingSyst, fillstyle=0)
SetObjectStyle(gPtIntCrossSecFDSyst, fillstyle=0)
SetObjectStyle(gPtIntCrossSecBRSyst, fillstyle=0)
SetObjectStyle(gPtIntCrossSecLumiSyst, fillstyle=0)
SetObjectStyle(gPtIntCrossSecExtrapSyst, fillcolor=kRed, fillalpha=0.8, fillstyle=1000)

hPtIntCrossSecStat.SetBinContent(1, totCrossSec)
hPtIntCrossSecStat.SetBinError(1, totCrossSecStat)
gPtIntCrossSecDataSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecDataSyst.SetPointError(0, 0.4, 0.4, totCrossSecDataSyst, totCrossSecDataSyst)
gPtIntCrossSecDataCorrSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecDataCorrSyst.SetPointError(0, 0.4, 0.4, totCrossSecCorrSyst, totCrossSecCorrSyst)
gPtIntCrossSecDataUncorrSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecDataUncorrSyst.SetPointError(0, 0.4, 0.4, totCrossSecUncorrSyst, totCrossSecUncorrSyst)
gPtIntCrossSecDataSystWoTrFDBRLumi.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecDataSystWoTrFDBRLumi.SetPointError(0, 0.4, 0.4,
                                                 totCrossSecDataSystWoTrFDBRLumi, totCrossSecDataSystWoTrFDBRLumi)
gPtIntCrossSecTrackingSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecTrackingSyst.SetPointError(0, 0.4, 0.4, totCrossSecTrSyst, totCrossSecTrSyst)
gPtIntCrossSecFDSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecFDSyst.SetPointError(0, 0.4, 0.4, totCrossSecFDSyst, totCrossSecFDSyst)
gPtIntCrossSecBRSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecBRSyst.SetPointError(0, 0.4, 0.4, totCrossSecBRSyst, totCrossSecBRSyst)
gPtIntCrossSecLumiSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecLumiSyst.SetPointError(0, 0.4, 0.4, totCrossSecLumiSyst, totCrossSecLumiSyst)
gPtIntCrossSecExtrapSyst.SetPoint(0, 1, totCrossSec)
gPtIntCrossSecExtrapSyst.SetPointError(0, 0.4, 0.4, totCrossSecExtrapSystLow, totCrossSecExtrapSystHigh)

gPtIntCrossSecExtrapSystPart = []
if 'correlated' in inputCfg['extrap']['ptint']['combination']: # includes both correlated and uncorrelated
    for iFONLL, _ in enumerate(extrFactorAlt):
        gPtIntCrossSecExtrapSystPart.append(TGraphAsymmErrors(1))
        gPtIntCrossSecExtrapSystPart[iFONLL].SetName(f'gPtIntCrossSecExtrapSystVar{iFONLL}')
        SetObjectStyle(gPtIntCrossSecExtrapSystPart[iFONLL], color=kGray+1, fillstyle=1000)
        gPtIntCrossSecExtrapSystPart[iFONLL].SetPoint(0, 1, totCrossSec)
        gPtIntCrossSecExtrapSystPart[iFONLL].SetPointError(0, 0.4, 0.4,
                                                           extrFactorFin[f'uncLowAlt{iFONLL}'] * visCrossSec,
                                                           extrFactorFin[f'uncHighAlt{iFONLL}'] * visCrossSec)

    gPtIntCrossSecExtrapSystPart.append(TGraphAsymmErrors(1))
    gPtIntCrossSecExtrapSystPart[-1].SetName('gPtIntCrossSecExtrapSystVarFONLL')
    SetObjectStyle(gPtIntCrossSecExtrapSystPart[-1], color=kGray+1, fillstyle=1000)
    gPtIntCrossSecExtrapSystPart[-1].SetPoint(0, 1, totCrossSec)
    gPtIntCrossSecExtrapSystPart[-1].SetPointError(0, 0.4, 0.4,
                                                   extrFactorFin['uncLowFONLL'] * visCrossSec,
                                                   extrFactorFin['uncHighFONLL'] * visCrossSec)

cPtInt = TCanvas('cPtInt', '', 500, 500)
cPtInt.DrawFrame(0.5, 0., 1.5, hPtIntCrossSecStat.GetMaximum()*1.5, ';;d#sigma/d#it{y} (#mub)')
gPtIntCrossSecExtrapSyst.Draw('2')
gPtIntCrossSecDataSyst.Draw('2')
hPtIntCrossSecStat.Draw('same')
cPtInt.Modified()
cPtInt.Update()

outFileNonPrompt = TFile(
    f'/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/Extrapolation/NonPrompt{mesonName}_PtIntegr{ptMin:.0f}-{ptMax:.0f}_{suffix}_FONLLextrap_final.root',
    'recreate')
cPtInt.Write()
hPtIntCrossSecStat.Write()
gPtIntCrossSecDataSyst.Write()
gPtIntCrossSecDataCorrSyst.Write()
gPtIntCrossSecDataUncorrSyst.Write()
gPtIntCrossSecDataSystWoTrFDBRLumi.Write()
gPtIntCrossSecTrackingSyst.Write()
gPtIntCrossSecFDSyst.Write()
gPtIntCrossSecBRSyst.Write()
gPtIntCrossSecLumiSyst.Write()
gPtIntCrossSecExtrapSyst.Write()
for graph in gPtIntCrossSecExtrapSystPart:
    graph.Write()
outFileNonPrompt.Close()

# hbbbarCrossSecStat = TH1D('hbbbarCrossSecStat', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)', 1, 0.5, 1.5)
# gbbbarCrossSecDataSyst, gbbbarCrossSecDataCorrSyst, gbbbarCrossSecDataUncorrSyst, \
#     gbbbarCrossSecDataSystWoTrFDBRLumi, gbbbarCrossSecTrackingSyst, gbbbarCrossSecFDSyst, gbbbarCrossSecBRSyst, \
#         gbbbarCrossSecLumiSyst, gbbbarCrossSecExtrapSyst, gbbbarCrossSecRapExtrapSyst = \
#             (TGraphAsymmErrors(0) for _ in range(10))
# gbbbarCrossSecDataSyst.SetNameTitle('gbbbarCrossSecDataSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecDataCorrSyst.SetNameTitle('gbbbarCrossSecDataCorrSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecDataUncorrSyst.SetNameTitle('gbbbarCrossSecDataUncorrSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecDataSystWoTrFDBRLumi.SetNameTitle('gbbbarCrossSecDataSystWoTrFDBRLumi',
#                                                 ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecTrackingSyst.SetNameTitle('gbbbarCrossSecTrackingSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecFDSyst.SetNameTitle('gbbbarCrossSecFDSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecBRSyst.SetNameTitle('gbbbarCrossSecBRSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecLumiSyst.SetNameTitle('gbbbarCrossSecLumiSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecExtrapSyst.SetNameTitle('gbbbarCrossSecExtrapSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# gbbbarCrossSecRapExtrapSyst.SetNameTitle('gbbbarCrossSecRapExtrapSyst', ';;d#sigma_{b#bar{b}}/d#it{y} (#mub)')
# SetObjectStyle(hbbbarCrossSecStat, color=kBlack)
# SetObjectStyle(gbbbarCrossSecDataSyst, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecDataCorrSyst, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecDataUncorrSyst, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecDataSystWoTrFDBRLumi, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecTrackingSyst, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecFDSyst, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecBRSyst, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecLumiSyst, color=kBlack, fillstyle=0)
# SetObjectStyle(gbbbarCrossSecExtrapSyst, color=kGray+1, fillstyle=1000)
# SetObjectStyle(gbbbarCrossSecRapExtrapSyst, color=kGray+1, fillstyle=1000)

# gbbbarCrossSecExtrapSystPart = []
# if 'correlated' in inputCfg['extrap']['bbbar']['combination']: # includes both correlated and uncorrelated
#     for iFONLL, _ in enumerate(extrFactorbbbarAlt):
#         gbbbarCrossSecExtrapSystPart.append(TGraphAsymmErrors(0))
#         gbbbarCrossSecExtrapSystPart[iFONLL].SetName(f'gbbbarCrossSecExtrapSystVar{iFONLL}')
#         SetObjectStyle(gbbbarCrossSecExtrapSystPart[iFONLL], color=kGray+1, fillstyle=1000)
#         gbbbarCrossSecExtrapSystPart[iFONLL].SetPoint(0, 1, totbbbarCrossSec)
#         gbbbarCrossSecExtrapSystPart[iFONLL].SetPointError(0, 0.4, 0.4,
#                                                            extrFactorbbbarFin[f'uncLowAlt{iFONLL}'] * visCrossSec,
#                                                            extrFactorbbbarFin[f'uncHighAlt{iFONLL}'] * visCrossSec)

#     gbbbarCrossSecExtrapSystPart.append(TGraphAsymmErrors(0))
#     gbbbarCrossSecExtrapSystPart[-1].SetName('gbbbarCrossSecExtrapSystVarFONLL')
#     SetObjectStyle(gbbbarCrossSecExtrapSystPart[-1], color=kGray+1, fillstyle=1000)
#     gbbbarCrossSecExtrapSystPart[-1].SetPoint(0, 1, totbbbarCrossSec)
#     gbbbarCrossSecExtrapSystPart[-1].SetPointError(0, 0.4, 0.4,
#                                                    extrFactorbbbarFin['uncLowFONLL'] * visCrossSec,
#                                                    extrFactorbbbarFin['uncHighFONLL'] * visCrossSec)

# hbbbarCrossSecStat.SetBinContent(1, totbbbarCrossSec)
# hbbbarCrossSecStat.SetBinError(1, totbbbarCrossSecStat)
# gbbbarCrossSecDataSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecDataSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecDataSyst, totbbbarCrossSecDataSyst)
# gbbbarCrossSecDataCorrSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecDataCorrSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecCorrSyst, totbbbarCrossSecCorrSyst)
# gbbbarCrossSecDataUncorrSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecDataUncorrSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecUncorrSyst, totbbbarCrossSecUncorrSyst)
# gbbbarCrossSecDataSystWoTrFDBRLumi.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecDataSystWoTrFDBRLumi.SetPointError(0, 0.4, 0.4, totbbbarCrossSecDataSystWoTrFDBRLumi,
#                                                  totbbbarCrossSecDataSystWoTrFDBRLumi)
# gbbbarCrossSecTrackingSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecTrackingSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecTrSyst, totbbbarCrossSecTrSyst)
# gbbbarCrossSecFDSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecFDSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecFDSyst, totbbbarCrossSecFDSyst)
# gbbbarCrossSecBRSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecBRSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecBRSyst, totbbbarCrossSecBRSyst)
# gbbbarCrossSecLumiSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecLumiSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecLumiSyst, totbbbarCrossSecLumiSyst)
# gbbbarCrossSecExtrapSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecExtrapSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecExtrapSystLow, totbbbarCrossSecExtrapSystHigh)
# gbbbarCrossSecRapExtrapSyst.SetPoint(0, 1, totbbbarCrossSec)
# gbbbarCrossSecRapExtrapSyst.SetPointError(0, 0.4, 0.4, totbbbarCrossSecExtrapRapSyst, totbbbarCrossSecExtrapRapSyst)
# gbbbarCrossSecExtrapSystTot = gbbbarCrossSecExtrapSyst.Clone('gbbbarCrossSecExtrapSystTot')
# gbbbarCrossSecExtrapSystTot.SetPointError(0, 0.4, 0.4,
#                                           np.sqrt(totbbbarCrossSecExtrapRapSyst**2 + totbbbarCrossSecExtrapSystLow**2),
#                                           np.sqrt(totbbbarCrossSecExtrapRapSyst**2 + totbbbarCrossSecExtrapSystHigh**2))

# cbbbar = TCanvas('cbbbar', '', 500, 500)
# cbbbar.DrawFrame(0.5, 0., 1.5, hbbbarCrossSecStat.GetMaximum()*1.5, ';;d#sigma/d#it{y} (#mub)')
# gbbbarCrossSecExtrapSystTot.Draw('2')
# gbbbarCrossSecDataSyst.Draw('2')
# hbbbarCrossSecStat.Draw('same')
# cbbbar.Modified()
# cbbbar.Update()


# if not inputCfg['extrap']['bbbar']['rapiditycorr']['powheg']['apply']:
#     suffix='_no_bbbar_rapidity_corr'

# outFilebbbar = TFile(f'ptintcrosssec/nonpromptD/CrossSection_bbbar_FromNonPrompt{mesonName}_FONLLextrap{suffix}.root',
#                      'recreate')
# cbbbar.Write()
# hbbbarCrossSecStat.Write()
# gbbbarCrossSecDataSyst.Write()
# gbbbarCrossSecDataCorrSyst.Write()
# gbbbarCrossSecDataUncorrSyst.Write()
# gbbbarCrossSecDataSystWoTrFDBRLumi.Write()
# gbbbarCrossSecTrackingSyst.Write()
# gbbbarCrossSecFDSyst.Write()
# gbbbarCrossSecBRSyst.Write()
# gbbbarCrossSecLumiSyst.Write()
# gbbbarCrossSecExtrapSyst.Write()
# for graph in gbbbarCrossSecExtrapSystPart:
#     graph.Write()
# gbbbarCrossSecRapExtrapSyst.Write()
# outFilebbbar.Close()

# write also txt file
# pylint: disable=anomalous-backslash-in-string,line-too-long
outFileNonPromptTxt = open(f'/home/abigot/AnalysisNonPromptDplus/Run2pPb5Tev/4_Analysis/5_CrossSection/Extrapolation/NonPrompt{mesonName}_PtIntegr{ptMin:.0f}-{ptMax:.0f}{suffix}_final.txt', 'w')
outFileNonPromptTxt.write(f'pt [{ptMin},{ptMax}], sigma = {visCrossSec:.5f} +- {visCrossSecStat:.5f} (stat) +- '
                          f'{visCrossSecDataSystWoTrFDBRLumi:.5f} (systMinusTr) +- {visCrossSecTrSyst:.5f} (systTr) '
                          f'+- {visCrossSecFDSyst:.5f} (systFD) +- {visCrossSecLumiSyst:.5f} (lumi) '
                          f'+- {visCrossSecBRSyst:.5f} (BR)\n\n')
outFileNonPromptTxt.write('------\n')
outFileNonPromptTxt.write(f'pt [{ptMin},{ptMax}], sigma = {visCrossSec:.1f} \pm {visCrossSecStat:.1f} '
                          f'(\mathrm{{stat}}) \pm {np.sqrt(visCrossSecDataSyst**2-visCrossSecLumiSyst**2-visCrossSecBRSyst**2):.1f} (\mathrm{{syst}}) \pm '
                          f'{visCrossSecLumiSyst:.1f} (\mathrm{{lumi}}) \pm '
                          f'{visCrossSecBRSyst:.1f} (\mathrm{{BR}})\n\n')
outFileNonPromptTxt.write('------\n')
outFileNonPromptTxt.write(f'Extrapolation factor to full pt = {extrFactor["central"]:.5f} +{extrFactor["uncHigh"]:.5f}'
                          f' -{extrFactor["uncLow"]:.5f}\n')
outFileNonPromptTxt.write('------\n\n')
outFileNonPromptTxt.write('------\n')
outFileNonPromptTxt.write(f'pt [0,inf], sigma = {totCrossSec:.5f} +- {totCrossSecStat:.5f} (stat) +- '
                          f'{totCrossSecDataSystWoTrFDBRLumi:.5f} (systMinusTr) +- {totCrossSecTrSyst:.5f} (systTr) '
                          f'+- {totCrossSecFDSyst:.5f} (systFD) +- {totCrossSecLumiSyst:.5f} (lumi) '
                          f'+- {totCrossSecBRSyst:.5f} (BR) + {totCrossSecExtrapSystHigh:.5f} - '
                          f'{totCrossSecExtrapSystLow:.5f} (systExtrap)\n\n')
outFileNonPromptTxt.write('------\n')
outFileNonPromptTxt.write(f'pt [0,inf], sigma = {totCrossSec:.1f} \pm {totCrossSecStat:.1f} '
                          f'(\mathrm{{stat}}) \pm {np.sqrt(totCrossSecDataSyst**2-totCrossSecBRSyst**2-totCrossSecLumiSyst**2):.1f} (\mathrm{{syst}}) \pm '
                          f'{totCrossSecLumiSyst:.1f} (\mathrm{{lumi}}) \pm '
                          f'{totCrossSecBRSyst:.1f} (\mathrm{{BR}}) ^{{+{totCrossSecExtrapSystHigh:.1f}}}'
                          f'_{{-{totCrossSecExtrapSystLow:.1f}}} (\mathrm{{extrap}})\n\n')
# outFileNonPromptTxt.write('------\n')
# outFileNonPromptTxt.write(f'Extrapolation factor to ds_bbbar/dy = {extrFactorbbbarFin["central"]:.5f} '
#                           f'+{extrFactorbbbarFin["uncHigh"]:.5f} -{extrFactorbbbarFin["uncLow"]:.5f}\n')
# outFileNonPromptTxt.write('------\n\n')
# outFileNonPromptTxt.write('------\n')
# outFileNonPromptTxt.write(f'ds/dy_bbbar = {totbbbarCrossSec:.1f} \pm {totbbbarCrossSecStat:.1f} '
#                           f'(\mathrm{{stat}}) \pm {np.sqrt(totbbbarCrossSecDataSyst**2-totbbbarCrossSecLumiSyst**2-totbbbarCrossSecBRSyst**2):.1f} (\mathrm{{syst}}) \pm '
#                           f'{totbbbarCrossSecLumiSyst:.1f} (\mathrm{{lumi}}) \pm '
#                           f'{totbbbarCrossSecBRSyst:.1f} (\mathrm{{BR}}) '
#                           f'^{{+{totbbbarCrossSecExtrapSystHigh:.1f}}}'
#                           f'_{{-{totbbbarCrossSecExtrapSystLow:.1f}}} (\mathrm{{extrap}})'
#                           f' \pm {totbbbarCrossSecExtrapRapSyst:.1f} (\mathrm{{rap. shape}})\n')
# outFileNonPromptTxt.close()

input('Press enter to exit')
