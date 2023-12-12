'''
Script with miscellanea utils methods for the analysis
'''

import ctypes
import numpy as np
from ROOT import TH1F, TF1, TList, TGraph, TGraphErrors, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module
from ROOT import kGreen, kAzure, kOrange, kMagenta, kGray # pylint: disable=import-error,no-name-in-module


def MergeHists(listOfHists):
    '''
    Method to merge histos

    Parameters
    ----------

    - listOfHists: python list of histos

    Returns
    ----------

    - hMerged: merged histo

    '''
    listMerge = TList()
    for iHist, hist in enumerate(listOfHists):
        if iHist == 0:
            hMerged = hist.Clone()
        else:
            listMerge.Add(hist)
    hMerged.Merge(listMerge)
    return hMerged


def ComputeRatioDiffBins(hNum, hDen, uncOpt=''):
    '''
    Method to compute ratio between histograms with different bins (but compatible)

    Parameters
    ----------

    - hNum: histogram for numerator
    - hDen: histogram for denominator
    - uncOpt: uncertainty option as in ROOT.TH1.Divide

    Returns
    ----------

    - hRatio: ratio histogram

    '''

    ptMinNum = hNum.GetBinLowEdge(1)
    ptMaxNum = hNum.GetXaxis().GetBinUpEdge(hNum.GetNbinsX())
    ptMinDen = hDen.GetBinLowEdge(1)
    ptMaxDen = hDen.GetXaxis().GetBinUpEdge(hDen.GetNbinsX())
    if ptMinNum < ptMinDen:
        ptMin = ptMinDen
    else:
        ptMin = ptMinNum
    if ptMaxNum > ptMaxDen:
        ptMax = ptMaxDen
    else:
        ptMax = ptMaxNum

    if hNum.GetNbinsX() < hDen.GetNbinsX():
        if np.array(hNum.GetXaxis().GetXbins(), 'd').any(): # variable binning
            ptLimsRatio = np.array(hNum.GetXaxis().GetXbins(), 'd')
        else: # constant binning
            binWidth = hNum.GetBinWidth(1)
            ptLimsRatio = np.array([ptMinDen + iBin * binWidth for iBin in range(hNum.GetNbinsX()+1)], 'd')
    else:
        if np.array(hDen.GetXaxis().GetXbins(), 'd').any(): # variable binning
            ptLimsRatio = np.array(hDen.GetXaxis().GetXbins(), 'd')
        else: # constant binning
            binWidth = hDen.GetBinWidth(1)
            ptLimsRatio = np.array([ptMinDen + iBin * binWidth for iBin in range(hDen.GetNbinsX()+1)], 'd')
    ptLimsRatio = ptLimsRatio[(ptLimsRatio >= ptMin) & (ptLimsRatio <= ptMax)]
    nPtBins = len(ptLimsRatio)-1

    hRatio = TH1F('hRatio', f';{hNum.GetXaxis().GetTitle()};ratio', nPtBins, ptLimsRatio)
    hNumReb = TH1F('hNumReb', '', nPtBins, ptLimsRatio)
    hDenReb = TH1F('hDenReb', '', nPtBins, ptLimsRatio)

    for iPtRatio in range(1, hRatio.GetNbinsX()+1):
        deltaPt = ptLimsRatio[iPtRatio]-ptLimsRatio[iPtRatio-1]
        num, numUnc, den, denUnc = (0 for _ in range(4))
        for iPtNum in range(1, hNum.GetNbinsX()+1):
            if hNum.GetBinLowEdge(iPtNum) >= ptLimsRatio[iPtRatio-1] and \
                hNum.GetXaxis().GetBinUpEdge(iPtNum) <= ptLimsRatio[iPtRatio]:
                num += hNum.GetBinContent(iPtNum) * hNum.GetBinWidth(iPtNum)
                numUnc += hNum.GetBinError(iPtNum)**2 * hNum.GetBinWidth(iPtNum)**2 # considered uncorrelated
        hNumReb.SetBinContent(iPtRatio, num/deltaPt)
        hNumReb.SetBinError(iPtRatio, np.sqrt(numUnc)/deltaPt)
        for iPtDen in range(1, hDen.GetNbinsX()+1):
            if hDen.GetBinLowEdge(iPtDen) >= ptLimsRatio[iPtRatio-1] and \
                hDen.GetXaxis().GetBinUpEdge(iPtDen) <= ptLimsRatio[iPtRatio]:
                den += hDen.GetBinContent(iPtDen) * hDen.GetBinWidth(iPtDen)
                denUnc += hDen.GetBinError(iPtDen)**2 * hDen.GetBinWidth(iPtDen)**2 # considered uncorrelated
        hDenReb.SetBinContent(iPtRatio, den/deltaPt)
        hDenReb.SetBinError(iPtRatio, np.sqrt(denUnc)/deltaPt)

    hRatio.Divide(hNumReb, hDenReb, 1., 1., uncOpt)

    return hRatio


def SumHistosDiffBins(h1, h2, uncOpt='uncorr'):
    '''
    Method to sum histograms with different bins (but compatible)

    Parameters
    ----------

    - h1: first histogram
    - h2: second histogram
    - uncOpt: uncertainty option (uncorr, corr)

    Returns
    ----------

    - hSum: sum histogram

    '''

    ptMinNum = h1.GetBinLowEdge(1)
    ptMaxNum = h1.GetXaxis().GetBinUpEdge(h1.GetNbinsX())
    ptMinDen = h2.GetBinLowEdge(1)
    ptMaxDen = h2.GetXaxis().GetBinUpEdge(h2.GetNbinsX())
    if ptMinNum < ptMinDen:
        ptMin = ptMinDen
    else:
        ptMin = ptMinNum
    if ptMaxNum > ptMaxDen:
        ptMax = ptMaxDen
    else:
        ptMax = ptMaxNum

    if h1.GetNbinsX() < h2.GetNbinsX():
        if np.array(h1.GetXaxis().GetXbins(), 'd').any(): # variable binning
            ptLimsSum = np.array(h1.GetXaxis().GetXbins(), 'd')
        else: # constant binning
            binWidth = h1.GetBinWidth(1)
            ptLimsSum = np.array([ptMinDen + iBin * binWidth for iBin in range(h1.GetNbinsX()+1)], 'd')
    else:
        if np.array(h2.GetXaxis().GetXbins(), 'd').any(): # variable binning
            ptLimsSum = np.array(h2.GetXaxis().GetXbins(), 'd')
        else: # constant binning
            binWidth = h2.GetBinWidth(1)
            ptLimsSum = np.array([ptMinDen + iBin * binWidth for iBin in range(h2.GetNbinsX()+1)], 'd')
    ptLimsSum = ptLimsSum[(ptLimsSum >= ptMin) & (ptLimsSum <= ptMax)]
    nPtBins = len(ptLimsSum)-1

    h1Reb = TH1F('h1Reb', f';{h1.GetXaxis().GetTitle()};{h1.GetYaxis().GetTitle()}', nPtBins, ptLimsSum)
    h2Reb = TH1F('h2Reb', f';{h2.GetXaxis().GetTitle()};{h2.GetYaxis().GetTitle()}', nPtBins, ptLimsSum)
    hSum = TH1F('hSum', f';{h1.GetXaxis().GetTitle()};sum', nPtBins, ptLimsSum)

    for iPtSum in range(1, hSum.GetNbinsX()+1):
        deltaPt = ptLimsSum[iPtSum]-ptLimsSum[iPtSum-1]
        y1, y1Unc, y2, y2Unc = (0 for _ in range(4))
        for iPt1 in range(1, h1.GetNbinsX()+1):
            if h1.GetBinLowEdge(iPt1) >= ptLimsSum[iPtSum-1] and \
                h1.GetXaxis().GetBinUpEdge(iPt1) <= ptLimsSum[iPtSum]:
                y1 += h1.GetBinContent(iPt1) * h1.GetBinWidth(iPt1)
                y1Unc += h1.GetBinError(iPt1)**2 * h1.GetBinWidth(iPt1)**2 # considered uncorrelated
        y1 /= deltaPt
        y1Unc = np.sqrt(y1Unc)/deltaPt
        for iPt2 in range(1, h2.GetNbinsX()+1):
            if h2.GetBinLowEdge(iPt2) >= ptLimsSum[iPtSum-1] and \
                h2.GetXaxis().GetBinUpEdge(iPt2) <= ptLimsSum[iPtSum]:
                y2 += h2.GetBinContent(iPt2) * h2.GetBinWidth(iPt2)
                y2Unc += h2.GetBinError(iPt2)**2 * h2.GetBinWidth(iPt2)**2 # considered uncorrelated
        y2 /= deltaPt
        y2Unc = np.sqrt(y2Unc)/deltaPt

        ySum = y1 + y2
        if uncOpt == 'uncorr':
            ySumUnc = np.sqrt(y1Unc**2 + y2Unc**2)
        else:
            ySumUnc = y1Unc + y2Unc

        h1Reb.SetBinContent(iPtSum, y1)
        h1Reb.SetBinError(iPtSum, y1Unc)
        h2Reb.SetBinContent(iPtSum, y2)
        h2Reb.SetBinError(iPtSum, y2Unc)
        hSum.SetBinContent(iPtSum, ySum)
        hSum.SetBinError(iPtSum, ySumUnc)

    return hSum, h1Reb, h2Reb


def ScaleGraph(graph, scaleFactor):
    '''
    Helper method to scale a TGraph

    Parameters
    ----------

    - graph: graph to scale
    - scaleFactor: scale factor
    '''
    for iPt in range(graph.GetN()):
        x, y = ctypes.c_double(), ctypes.c_double()
        graph.GetPoint(iPt, x, y)
        graph.SetPoint(iPt, x.value, y.value * scaleFactor)
        if isinstance(graph, TGraphAsymmErrors):
            yUncLow = graph.GetErrorYlow(iPt)
            yUncHigh = graph.GetErrorYhigh(iPt)
            graph.SetPointEYlow(iPt, yUncLow * scaleFactor)
            graph.SetPointEYhigh(iPt, yUncHigh * scaleFactor)
        elif isinstance(graph, TGraphErrors):
            yUncLow = graph.GetErrorYlow(iPt)
            yUncHigh = graph.GetErrorYhigh(iPt)
            graph.SetPointError(iPt, graph.GetErrorX(iPt), graph.GetErrorY(iPt) * scaleFactor)
        elif isinstance(graph, TGraph):
            continue


def ComputeRatioGraph(gNum, gDen, useDenUnc=True):
    '''
    Helper method to divide two TGraph (assuming same binning)

    Parameters
    ----------

    - gNum: graph to divide (numerator)
    - gDen: graph to divide (denominator)

    Returns
    ----------

    - gRatio: resulting graph
    '''
    if gNum.GetN() != gDen.GetN():
        print('ERROR: only graphs with same number of bins can be divided!')
        return None

    gRatio = TGraphAsymmErrors(0)
    for iPt in range(gNum.GetN()):
        x, num = ctypes.c_double(), ctypes.c_double()
        xd, den = ctypes.c_double(), ctypes.c_double()
        gNum.GetPoint(iPt, x, num)
        xUncLow = gNum.GetErrorXlow(iPt)
        xUncHigh = gNum.GetErrorXhigh(iPt)
        numUncLow = gNum.GetErrorYlow(iPt)
        numUncHigh = gNum.GetErrorYhigh(iPt)
        gDen.GetPoint(iPt, xd, den)
        denUncLow = gDen.GetErrorYlow(iPt)
        denUncHigh = gDen.GetErrorYhigh(iPt)

        ratio, ratioUncLow, ratioUncHigh = 0., 0., 0.
        if num.value != 0. and den.value != 0.:
            ratio = num.value/den.value
            if useDenUnc:
                ratioUncLow = np.sqrt((numUncLow/num.value)**2 + (denUncLow/den.value)**2) * ratio
                ratioUncHigh = np.sqrt((numUncHigh/num.value)**2 + (denUncHigh/den.value)**2) * ratio
            else:
                ratioUncLow = numUncLow / num.value * ratio
                ratioUncHigh = numUncHigh / num.value * ratio

        gRatio.SetPoint(iPt, x.value, ratio)
        gRatio.SetPointError(iPt, xUncLow, xUncHigh, ratioUncLow, ratioUncHigh)

    return gRatio


def DivideGraphByHisto(gNum, hDen, useHistoUnc=True):
    '''
    Helper method to divide a TGraph by a TH1 (assuming same binning)

    Parameters
    ----------

    - gNum: graph to divide (numerator)
    - hDen: histogram (denominator)

    Returns
    ----------

    - gRatio: resulting graph
    '''
    if gNum.GetN() != hDen.GetNbinsX():
        print('ERROR: only graphs and histos with same number of bins can be divided!')
        return None

    gRatio = TGraphAsymmErrors(0)
    for iPt in range(gNum.GetN()):
        x, num = ctypes.c_double(), ctypes.c_double()
        gNum.GetPoint(iPt, x, num)
        xUncLow = gNum.GetErrorXlow(iPt)
        xUncHigh = gNum.GetErrorXhigh(iPt)
        numUncLow = gNum.GetErrorYlow(iPt)
        numUncHigh = gNum.GetErrorYhigh(iPt)
        ptBinHisto = hDen.GetXaxis().FindBin(x.value)
        den = hDen.GetBinContent(ptBinHisto)
        if useHistoUnc:
            ratioUncLow = np.sqrt((numUncLow/num.value)**2 + (hDen.GetBinError(ptBinHisto)/den)**2) * num.value/den
            ratioUncHigh = np.sqrt((numUncHigh/num.value)**2 + (hDen.GetBinError(ptBinHisto)/den)**2) * num.value/den
        else:
            ratioUncLow = numUncLow/num.value * num.value/den
            ratioUncHigh = numUncHigh/num.value * num.value/den
        gRatio.SetPoint(iPt, x.value, num.value/den)
        gRatio.SetPointError(iPt, xUncLow, xUncHigh, ratioUncLow, ratioUncHigh)

    return gRatio


def AverageDmesonRatio(hStat, gSystPtCorr, gSystPtUncorr, gSystBR=None, gSystExtrap=None, opt='fit'):
    '''
    Helper method to fit D-meson ratio with a constant function
    considering correlated and uncorrelated syst unc properly

    Parameters
    ----------

    - hStat: histogram with statistical uncertainties
    - gSystPtCorr: graph with pt correlated systematic uncertainties
    - gSystPtUncorr: graph with pt uncorrelated systematic uncertainties
    - gSystBR: graph with BR systematic uncertainties (not mandatory)
    - gSystExtrap: graph with extrapolation systematic uncertainties (not mandatory)
    - opt: option for average (fit, weighted)

    Returns
    ----------

    - averageRatio: dictionary with
        {cent, stat, uncorrSys, corrSysLow, corrSysHigh, BRSysLow, BRSysHigh, ExtrapSysLow, ExtrapSysHigh}
    - fAverageRatio: dictionary with fit functions
        {cent, stat, uncorrSys, corrSysLow, corrSysHigh, BRSysLow, BRSysHigh, ExtrapSysLow, ExtrapSysHigh}
    '''

    if opt not in ['fit', 'weighted']:
        print('ERROR: average option should be either fit or weighted! Return None')
        return None, None

    # anyway store result in TF1
    fStat = TF1('fStat', 'pol0', 0., 100.)
    fStat.SetLineColor(kGreen+2)
    fSysPtUncorr = TF1('fSysPtUncorr', 'pol0', 0., 100.)
    fSysPtUncorr.SetLineColor(kAzure+4)
    fStatPlusSysPtUncorr = TF1('fStatPlusSysPtUncorr', 'pol0', 0., 100.)
    fStatPlusSysPtUncorrShiftUp = TF1('fStatPlusSysPtUncorrShiftUp', 'pol0', 0., 100.)
    fStatPlusSysPtUncorrShiftUp.SetLineColor(kOrange+7)
    fStatPlusSysPtUncorrShiftDown = TF1('fStatPlusSysPtUncorrShiftDown', 'pol0', 0., 100.)
    fStatPlusSysPtUncorrShiftDown.SetLineColor(kOrange+7)
    if gSystBR:
        fStatPlusSysPtUncorrBRShiftUp = TF1('fStatPlusSysPtUncorrBRShiftUp', 'pol0', 0., 100.)
        fStatPlusSysPtUncorrBRShiftUp.SetLineColor(kMagenta+1)
        fStatPlusSysPtUncorrBRShiftDown = TF1('fStatPlusSysPtUncorrBRShiftDown', 'pol0', 0., 100.)
        fStatPlusSysPtUncorrBRShiftDown.SetLineColor(kMagenta+1)
    if gSystExtrap:
        fStatPlusSysPtUncorrExtrapShiftUp = TF1('fStatPlusSysPtUncorrExtrapShiftUp', 'pol0', 0., 100.)
        fStatPlusSysPtUncorrExtrapShiftUp.SetLineColor(kGray+1)
        fStatPlusSysPtUncorrExtrapShiftDown = TF1('fStatPlusSysPtUncorrExtrapShiftDown', 'pol0', 0., 100.)
        fStatPlusSysPtUncorrExtrapShiftDown.SetLineColor(kGray+1)

    if opt == 'fit':
        gStatPlusSystPtUncorr = gSystPtUncorr.Clone('gStatPlusSystPtUncorr')
        gStatPlusSystPtUncorrShiftUp = gSystPtUncorr.Clone('gStatPlusSystPtUncorrShiftUp')
        gStatPlusSystPtUncorrShiftDown = gSystPtUncorr.Clone('gStatPlusSystPtUncorrShiftDown')
        if gSystBR: # BR uncertainty also pT corr
            gStatPlusSystPtUncorrBRShiftUp = gSystPtUncorr.Clone('gStatPlusSystPtUncorrBRShiftUp')
            gStatPlusSystPtUncorrBRShiftDown = gSystPtUncorr.Clone('gStatPlusSystPtUncorrBRShiftDown')
        if gSystExtrap:
            gStatPlusSystPtUncorrExtrapShiftUp = gSystPtUncorr.Clone('gStatPlusSystPtUncorrExtrapShiftUp')
            gStatPlusSystPtUncorrExtrapShiftDown = gSystPtUncorr.Clone('gStatPlusSystPtUncorrExtrapShiftDown')
        for iPt in range(gSystPtUncorr.GetN()):
            pt = hStat.GetBinCenter(iPt+1)
            cent = hStat.GetBinContent(iPt+1)
            statUnc = hStat.GetBinError(iPt+1)
            ptUnc = gSystPtUncorr.GetErrorXlow(iPt)
            sysUncorrLow = gSystPtUncorr.GetErrorYlow(iPt)
            sysUncorrHigh = gSystPtUncorr.GetErrorYhigh(iPt)
            sysCorrLow = gSystPtCorr.GetErrorYlow(iPt)
            sysCorrHigh = gSystPtCorr.GetErrorYhigh(iPt)
            gStatPlusSystPtUncorrShiftUp.SetPoint(iPt, pt, cent+sysCorrHigh)
            gStatPlusSystPtUncorrShiftDown.SetPoint(iPt, pt, cent-sysCorrLow)
            gStatPlusSystPtUncorr.SetPointError(iPt, ptUnc, ptUnc, np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                np.sqrt(statUnc**2 + sysUncorrHigh**2))
            gStatPlusSystPtUncorrShiftUp.SetPointError(iPt, ptUnc, ptUnc, np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                       np.sqrt(statUnc**2 + sysUncorrHigh**2))
            gStatPlusSystPtUncorrShiftDown.SetPointError(iPt, ptUnc, ptUnc, np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                         np.sqrt(statUnc**2 + sysUncorrHigh**2))
            if gSystBR:
                BRunc = gSystBR.GetErrorYlow(iPt)
                gStatPlusSystPtUncorrBRShiftUp.SetPoint(iPt, pt, cent+BRunc)
                gStatPlusSystPtUncorrBRShiftDown.SetPoint(iPt, pt, cent-BRunc)
                gStatPlusSystPtUncorrBRShiftUp.SetPointError(iPt, ptUnc, ptUnc,
                                                             np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                             np.sqrt(statUnc**2 + sysUncorrHigh**2))
                gStatPlusSystPtUncorrBRShiftDown.SetPointError(iPt, ptUnc, ptUnc,
                                                               np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                               np.sqrt(statUnc**2 + sysUncorrHigh**2))
            if gSystExtrap:
                ExtrapUncLow = gSystExtrap.GetErrorYlow(iPt)
                ExtrapUncHigh = gSystExtrap.GetErrorYhigh(iPt)
                gStatPlusSystPtUncorrExtrapShiftUp.SetPoint(iPt, pt, cent+ExtrapUncHigh)
                gStatPlusSystPtUncorrExtrapShiftDown.SetPoint(iPt, pt, cent-ExtrapUncLow)
                gStatPlusSystPtUncorrExtrapShiftUp.SetPointError(iPt, ptUnc, ptUnc,
                                                                 np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                                 np.sqrt(statUnc**2 + sysUncorrHigh**2))
                gStatPlusSystPtUncorrExtrapShiftDown.SetPointError(iPt, ptUnc, ptUnc,
                                                                   np.sqrt(statUnc**2 + sysUncorrLow**2),
                                                                   np.sqrt(statUnc**2 + sysUncorrHigh**2))
        hStat.Fit(fStat, '0Q')
        gSystPtUncorr.Fit(fSysPtUncorr, '0Q')
        gStatPlusSystPtUncorr.Fit(fStatPlusSysPtUncorr, '0Q')
        gStatPlusSystPtUncorrShiftUp.Fit(fStatPlusSysPtUncorrShiftUp, '0Q')
        gStatPlusSystPtUncorrShiftDown.Fit(fStatPlusSysPtUncorrShiftDown, '0Q')
        if gSystBR:
            gStatPlusSystPtUncorrBRShiftUp.Fit(fStatPlusSysPtUncorrBRShiftUp, '0Q')
            gStatPlusSystPtUncorrBRShiftDown.Fit(fStatPlusSysPtUncorrBRShiftDown, '0Q')
        if gSystExtrap:
            gStatPlusSystPtUncorrExtrapShiftUp.Fit(fStatPlusSysPtUncorrExtrapShiftUp, '0Q')
            gStatPlusSystPtUncorrExtrapShiftUp.Fit(fStatPlusSysPtUncorrExtrapShiftDown, '0Q')

    elif opt == 'weighted':
        avRatio, avRatioStat, avRatioUncorrSys, avRatioCorrSysLow, avRatioCorrSysHigh, \
            avRatioBRSysLow, avRatioBRSysHigh, avRatioExtrapSysLow, avRatioExtrapSysHigh, \
                sumOfWeights = (0. for _ in range(10))
        for iPt in range(hStat.GetNbinsX()):
            rat = hStat.GetBinContent(iPt+1)
            relStat = hStat.GetBinError(iPt+1) / rat
            relUncorrSys = gSystPtUncorr.GetErrorYlow(iPt) / rat
            relCorrSysLow = gSystPtCorr.GetErrorYlow(iPt) / rat
            relCorrSysHigh = gSystPtCorr.GetErrorYhigh(iPt) / rat
            weight = 1. / (relStat**2 + relUncorrSys**2)
            avRatio += rat * weight
            sumOfWeights += weight
            avRatioStat += (relStat*rat)**2 * weight**2
            avRatioUncorrSys += (relUncorrSys*rat)**2 * weight**2
            avRatioCorrSysLow += relCorrSysLow*rat * weight
            avRatioCorrSysHigh += relCorrSysHigh*rat * weight
            if gSystBR:
                BRSysLow = gSystBR.GetErrorYlow(iPt)
                BRSysHigh = gSystBR.GetErrorYhigh(iPt)
                avRatioBRSysLow += BRSysLow * weight
                avRatioBRSysHigh += BRSysHigh * weight
            if gSystExtrap:
                ExtrapSysLow = gSystBR.GetErrorYlow(iPt)
                ExtrapSysHigh = gSystBR.GetErrorYhigh(iPt)
                avRatioExtrapSysLow += ExtrapSysLow * weight
                avRatioExtrapSysHigh += ExtrapSysHigh * weight
        avRatio /= sumOfWeights
        avRatioStat = np.sqrt(avRatioStat/sumOfWeights**2)
        avRatioUncorrSys = np.sqrt(avRatioUncorrSys/sumOfWeights**2)
        avRatioCorrSysLow /= sumOfWeights
        avRatioCorrSysHigh /= sumOfWeights
        fStatPlusSysPtUncorr.SetParameter(0, avRatio)
        fStatPlusSysPtUncorr.SetParError(0, np.sqrt(avRatioStat**2 + avRatioUncorrSys**2))
        fStat.SetParameter(0, avRatio)
        fStat.SetParError(0, avRatioStat)
        fSysPtUncorr.SetParameter(0, avRatio)
        fSysPtUncorr.SetParError(0, avRatioUncorrSys)
        fStatPlusSysPtUncorrShiftDown.SetParameter(0, avRatio-avRatioCorrSysLow)
        fStatPlusSysPtUncorrShiftDown.SetParError(0, 0.)
        fStatPlusSysPtUncorrShiftUp.SetParameter(0, avRatio+avRatioCorrSysHigh)
        fStatPlusSysPtUncorrShiftUp.SetParError(0, 0.)
        if gSystBR:
            avRatioBRSysLow /= sumOfWeights
            avRatioBRSysHigh /= sumOfWeights
            fStatPlusSysPtUncorrBRShiftDown.SetParameter(0, avRatio-avRatioBRSysLow)
            fStatPlusSysPtUncorrBRShiftDown.SetParError(0, 0.)
            fStatPlusSysPtUncorrBRShiftUp.SetParameter(0, avRatio+avRatioBRSysHigh)
            fStatPlusSysPtUncorrBRShiftUp.SetParError(0, 0.)
        if gSystExtrap:
            avRatioExtrapSysLow /= sumOfWeights
            avRatioExtrapSysHigh /= sumOfWeights
            fStatPlusSysPtUncorrExtrapShiftDown.SetParameter(0, avRatio-avRatioExtrapSysLow)
            fStatPlusSysPtUncorrExtrapShiftDown.SetParError(0, 0.)
            fStatPlusSysPtUncorrExtrapShiftUp.SetParameter(0, avRatio+avRatioExtrapSysHigh)
            fStatPlusSysPtUncorrExtrapShiftUp.SetParError(0, 0.)

    averageRatio = {'cent': fStatPlusSysPtUncorr.GetParameter(0),
                    'stat': fStat.GetParError(0),
                    'uncorrSys': fSysPtUncorr.GetParError(0),
                    'corrSysLow': fStatPlusSysPtUncorr.GetParameter(0)-fStatPlusSysPtUncorrShiftDown.GetParameter(0),
                    'corrSysHigh': fStatPlusSysPtUncorrShiftUp.GetParameter(0)-fStatPlusSysPtUncorr.GetParameter(0)}

    fAverageRatio = {'cent': fStatPlusSysPtUncorr,
                     'stat': fStat,
                     'uncorrSys': fSysPtUncorr,
                     'corrSysLow': fStatPlusSysPtUncorrShiftDown,
                     'corrSysHigh': fStatPlusSysPtUncorrShiftUp}

    if gSystBR:
        averageRatio['BRSysLow'] = fStatPlusSysPtUncorr.GetParameter(0)-fStatPlusSysPtUncorrBRShiftDown.GetParameter(0)
        averageRatio['BRSysHigh'] = fStatPlusSysPtUncorrBRShiftUp.GetParameter(0)-fStatPlusSysPtUncorr.GetParameter(0)
        fAverageRatio['BRSysLow'] = fStatPlusSysPtUncorrBRShiftDown
        fAverageRatio['BRSysHigh'] = fStatPlusSysPtUncorrBRShiftUp
    if gSystExtrap:
        averageRatio['ExtrapSysLow'] = \
            fStatPlusSysPtUncorr.GetParameter(0)-fStatPlusSysPtUncorrExtrapShiftDown.GetParameter(0)
        averageRatio['ExtrapSysHigh'] = \
            fStatPlusSysPtUncorrExtrapShiftUp.GetParameter(0)-fStatPlusSysPtUncorr.GetParameter(0)
        fAverageRatio['ExtrapSysLow'] = fStatPlusSysPtUncorrExtrapShiftDown
        fAverageRatio['ExtrapSysHigh'] = fStatPlusSysPtUncorrExtrapShiftUp

    return averageRatio, fAverageRatio


def ComputePtExtrapolationFactor(hFONLLPred, ptMeasMin, ptMeasMax):
    '''
    Method for the computation of the pT extrapolation factor

    Inputs
    -------------
    - hFONLLPred: dictionary of histograms {central, min, max}
    - ptMeasMin: minimum pT of the measured cross section
    - ptMeasMax: minimum pT of the measured cross section

    Outputs
    -------------
    - extrapolation factor {central, min, max, uncLow, uncHigh}
    '''

    ptBinMinMeas = hFONLLPred['central'].GetXaxis().FindBin(ptMeasMin*1.0001)
    ptBinMaxMeas = hFONLLPred['central'].GetXaxis().FindBin(ptMeasMax*0.9999)
    # pt24 = hFONLLPred['central'].GetXaxis().FindBin(24*0.9999)
    visCrossSecTh, totCrossSecTh, extrapFactor = ({} for _ in range(3))
    # in case of variable binning,
    for prediction in hFONLLPred:
        visCrossSecTh[prediction] = hFONLLPred[prediction].Integral(ptBinMinMeas, ptBinMaxMeas, 'width')
        totCrossSecTh[prediction] = hFONLLPred[prediction].Integral('width')
        # totCrossSecTh[prediction] = hFONLLPred[prediction].Integral(1, pt24, 'width')
        extrapFactor[prediction] = totCrossSecTh[prediction] / visCrossSecTh[prediction]

    minExtrap = min([extrapFactor['central'], extrapFactor['min'], extrapFactor['max']])
    maxExtrap = max([extrapFactor['central'], extrapFactor['min'], extrapFactor['max']])
    extrapFactor['uncLow'] = extrapFactor['central']-minExtrap
    extrapFactor['uncHigh'] = maxExtrap-extrapFactor['central']

    return extrapFactor


def ComputebbbarExtrapolationFactor(FONLLPredNum, hFONLLPredDen, ptMeasMin, ptMeasMax, BRDmes):
    '''
    Method for the computation of the ds_bbbar/dy extrapolation factor

    Inputs
    -------------
    - FONLLPredNum: ds_bbbar/dy from FONLL {central, min, max}
    - hFONLLPredDen: dictionary of FONLL histograms for non-prompt D {central, min, max}
    - ptMeasMin: minimum pT of the measured cross section
    - ptMeasMax: minimum pT of the measured cross section
    - BRDmes: BR for D --> final state

    Outputs
    -------------
    - extrapolation factor {central, min, max, uncLow, uncHigh}
    '''

    ptBinMinMeas = hFONLLPredDen['central'].GetXaxis().FindBin(ptMeasMin*1.0001)
    ptBinMaxMeas = hFONLLPredDen['central'].GetXaxis().FindBin(ptMeasMax*0.9999)
    visCrossSecTh, totCrossSecTh, extrapFactor = ({} for _ in range(3))
    # in case of variable binning,
    for prediction in hFONLLPredDen:
        visCrossSecTh[prediction] = hFONLLPredDen[prediction].Integral(ptBinMinMeas, ptBinMaxMeas, 'width') / BRDmes
        totCrossSecTh[prediction] = hFONLLPredDen[prediction].Integral('width') / BRDmes
        extrapFactor[prediction] = FONLLPredNum[prediction][0] / visCrossSecTh[prediction]

    minExtrap = min([extrapFactor['central'], extrapFactor['min'], extrapFactor['max']])
    maxExtrap = max([extrapFactor['central'], extrapFactor['min'], extrapFactor['max']])
    extrapFactor['uncLow'] = extrapFactor['central']-minExtrap
    extrapFactor['uncHigh'] = maxExtrap-extrapFactor['central']

    return extrapFactor


def CombineExtrapFactors(extrFact, extrFactAlt, strategy='none'):
    '''
    Method for the combination of extrapolation factors

    Inputs
    -------------
    - extrFact: dictionary with default extrapolation factor {central, min, max, uncLow, uncHigh}
    - extrFactAlt: list of dictionaries with alternative extrapolation factors {central, min, max, uncLow, uncHigh}
    - strategy: extrapolation factor strategy (possibilities: none, envelope, correlated, uncorrelated)

    Outputs
    -------------
    - combined extrapolation factor {central, min, max, uncLow, uncHigh}
    '''

    if strategy not in ['none', 'envelope', 'correlated', 'uncorrelated']:
        print(f'\n\x1b[33mWARNING: {strategy} combination strategy for extrapolation factor not valid! '
              'Change to none\033[0m\n')
        strategy = 'none'

    extrFactFin = extrFact.copy()
    if strategy == 'none':
        return extrFactFin
    if strategy == 'envelope':
        for iAlt, _ in enumerate(extrFactAlt):
            if extrFactFin['min'] > extrFactAlt[iAlt]['min']:
                extrFactFin['min'] = extrFactAlt[iAlt]['min']
            if extrFactFin['max'] < extrFactAlt[iAlt]['max']:
                extrFactFin['max'] = extrFactAlt[iAlt]['max']
        extrFactFin['uncLow'] = extrFactFin['central']-extrFactFin['min']
        extrFactFin['uncHigh'] = extrFactFin['max']-extrFactFin['central']
    elif strategy != 'none':
        varLow, varHigh = [], []
        uncLow, uncHigh = 0., 0.
        for iAlt, _ in enumerate(extrFactAlt):
            if extrFactAlt[iAlt]['central'] > extrFact['central']:
                varLow.append(0.)
                varHigh.append(extrFactAlt[iAlt]['central']-extrFact['central'])
            else:
                varLow.append(extrFact['central']-extrFactAlt[iAlt]['central'])
                varHigh.append(0.)
            if strategy == 'correlated':
                uncLow += varLow[iAlt]
                uncHigh += varHigh[iAlt]
            elif strategy == 'uncorrelated':
                uncLow += varLow[iAlt]**2
                uncHigh += varHigh[iAlt]**2
        if strategy == 'uncorrelated':
            uncLow = np.sqrt(uncLow)
            uncHigh = np.sqrt(uncHigh)
        varLow.append(extrFact['uncLow']) # last one is FONLL
        varHigh.append(extrFact['uncHigh']) # last one is FONLL
        extrFactFin['uncLow'] = np.sqrt(uncLow**2 + varLow[-1]**2)
        extrFactFin['uncHigh'] = np.sqrt(uncHigh**2 + varLow[-1]**2)
        extrFactFin['min'] = extrFactFin['central'] - extrFactFin['uncLow']
        extrFactFin['max'] = extrFactFin['central'] + extrFactFin['uncHigh']
        for iVar, (varL, varH) in enumerate(zip(varLow, varHigh)):
            if iVar < len(varLow)-1:
                extrFactFin[f'uncLowAlt{iVar}'] = varL
                extrFactFin[f'uncHighAlt{iVar}'] = varH
            else:
                extrFactFin['uncLowFONLL'] = varL
                extrFactFin['uncHighFONLL'] = varH

    return extrFactFin
