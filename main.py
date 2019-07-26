

import os
import copy
import pickle
import logging
import ROOT
import math
import sys, getopt
import string
import os
from ROOT import TFile, TTree, TH1F, TH1D, TH2F, TH2D,  TCanvas, TLegend, TF1, gROOT, gSystem

import ROOT
import numpy as np
from array import array

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)
#gStyle.SetErrorX(1.)

gSystem.Load("/afs/cern.ch/user/l/lixia/workspace/RooUnfold/RooUnfold/libRooUnfold.so")

from ZoomInUnfold import ZoomInUnfold


def ftrue(X, par=[]):
    x = X[0]
    #return par[0]*math.exp(-x/100.0)+par[1]*ROOT.TMath.BreitWigner(x,91,2.5)+par[2]*ROOT.TMath.BreitWigner(x,125,0.005)
    return par[0]*math.exp(-x/500.0)+par[1]*ROOT.TMath.BreitWigner(x,90, 10)+par[2]*ROOT.TMath.BreitWigner(x,125, 1)

def generate_dataset(Icase=1):
    xmin = 20.
    xmax = 230.
    func = TF1('func', ftrue, xmin, xmax, 3)
    func.SetParameters(0.01,1,0.2)
    Cs = TCanvas('Cs', '', 10, 10, 800, 600)
    func.Draw()
    Cs.SaveAs('Cs_func.png')
    f = TFile('data_'+str(Icase)+'.root', 'recreate')
    tree = TTree('nominal', '')
    mtrue = array('f', [0])
    mrec = array('f', [0])
    tree.Branch('mtrue', mtrue, 'mtrue/F')
    tree.Branch('mrec', mrec, 'mrec/F')
    num = 10000
    for i in range(num):
        x0 = func.GetRandom(xmin, xmax)
        if Icase == 1:
            x = gRandom.Gaus(x0+5, 5)
        elif Icase == 2:
            x = gRandom.Gaus(x0+5, 10)
        elif Icase == 3:
            x = gRandom.Gaus(x0+5, 20)
        else:
            print('ERROR!!!')
            break
        mtrue[0]=x0
        mrec[0] = x
        tree.Fill()
    tree.Write()
    f.Save()
    return

def prepare_Axb(tree, nx=30, ny=30, wsg=0):
    nbinsx = nx
    nbinsy = ny
    xmin = 50.
    xmax = 200.
    dx = (xmax-xmin)/nbinsx
    dy = (xmax-xmin)/nbinsy

    hA = TH2D('hA', '', nbinsy, xmin, xmax, nbinsx, xmin, xmax)
    hA.Sumw2()
    #f = TFile('data.root', 'read')
    #tree = f.Get('nominal')
    tree.Draw('mtrue>>hx(%i,%.f,%.f)' % (nbinsx, xmin, xmax))
    tree.Draw('mrec>>hb(%i,%.f,%.f)' % (nbinsy, xmin, xmax))
    tree.Draw('mtrue>>hxp(%i,%.f,%.f)' % (nbinsx, xmin, xmax), 'mrec>%.f && mrec<%.f' %(xmin, xmax))
    tree.Draw('mrec>>hbp(%i,%.f,%.f)' % (nbinsx, xmin, xmax), 'mtrue>%.f && mtrue<%.f' %(xmin, xmax))
    for i in range(nbinsx): #loop true
        Ntrue = tree.GetEntries('mtrue > %.2f && mtrue < %.2f' % (xmin+i*dx, xmin+i*dx+dx))
        sumNrec = 0.
        for j in range(nbinsy): #loop rec
            if j == 0:
                #Nrec = tree.GetEntries('mtrue > %.2f && mtrue < %.2f && mrec < %.2f' % (xmin+i*dx, xmin+i*dx+dx, xmin+j*dy+dy))
                Nrec = tree.GetEntries('mtrue > %.2f && mtrue < %.2f && mrec > %.2f && mrec < %.2f' % (xmin+i*dx, xmin+i*dx+dx, xmin+j*dy, xmin+j*dy+dy))
            elif j == nbinsy-1:
                #Nrec = tree.GetEntries('mtrue > %.2f && mtrue < %.2f && mrec > %.2f' % (xmin+i*dx, xmin+i*dx+dx, xmin+j*dy))
                Nrec = tree.GetEntries('mtrue > %.2f && mtrue < %.2f && mrec > %.2f && mrec < %.2f' % (xmin+i*dx, xmin+i*dx+dx, xmin+j*dy, xmin+j*dy+dy))
            else:
                Nrec = tree.GetEntries('mtrue > %.2f && mtrue < %.2f && mrec > %.2f && mrec < %.2f' % (xmin+i*dx, xmin+i*dx+dx, xmin+j*dy, xmin+j*dy+dy))
            if wsg or 1:
                hA.SetBinContent(j+1, i+1, float(Nrec)/float(Ntrue))
            else:
                hA.SetBinContent(i+1, j+1, float(Nrec)/float(Ntrue))
            sumNrec += Nrec
        #print('point',i,sumNrec/Ntrue)
    hx = gROOT.FindObject('hx')
    hb = gROOT.FindObject('hb')
    hxp = gROOT.FindObject('hxp')
    hbp = gROOT.FindObject('hbp')
    hx.Sumw2()
    hxp.Sumw2()
    hb.Sumw2()
    hbp.Sumw2()
    if wsg:
        return hA, hx, hb, hxp, hbp
    else:
        return hA, hx, hb

def TH1Dclone(h):
    nbins = h.GetNbinsX()
    xmin = h.GetXaxis().GetXmin()
    xmax = h.GetXaxis().GetXmax()
    hnew = TH1D(h.GetName()+'_clone', '', nbins, xmin, xmax)
    for i in range(nbins):
        hnew.SetBinContent(i+1, h.GetBinContent(i+1))
        hnew.SetBinError(i+1, h.GetBinError(i+1))
    return hnew
def convert2hist(binning, L, name='tmp'):
    nbins = len(L)
    h = TH1F('h'+name, '', nbins, np.array(binning))
    for i in range(nbins):
        binwidth = h.GetBinWidth(i+1)
        h.SetBinContent(i+1, L[i][0]/binwidth)
        h.SetBinError(i+1, L[i][1]/binwidth)
    return h
def newunfold(Icase=1):
    print('Using ZoomInUnfold method....')
    xmin = 50.
    xmax = 200.
    nbins = 5
    nsplit = 20
    res = 2.
    if Icase == 1:
        nsplit = 20
        res = 2
    elif Icase == 2:
        nbins = 5
        nsplit = 25
        #nsplit = 2
        res = 10.
    dotoymc = 0 
    ntoys = 100
    sigfilepath = 'data_'+str(Icase)+'.root'
    sigtreename = 'nominal'
    datafilepath = sigfilepath
    datatreename = sigtreename
    unfold = ZoomInUnfold(sigfilepath, sigtreename, datafilepath, datatreename, 'mtrue', 'mrec', xmin, xmax, nbins, nsplit, res, dotoymc, ntoys)
    unfold.do_unfold()
    #binning, hb, hx, hout = unfold.get_output_hist()
    binning, hb, hx, hout = unfold.get_output()
    hRx, Rx = unfold.get_Rx()
    hb = convert2hist(binning, hb, 'rec')
    hx = convert2hist(binning, hx, 'true')
    hout = convert2hist(binning, hout, 'out')
    print('binning =', binning)
    if 0:
        print('binning =', binning)
        print('hb =', hb, hb.Integral(), hb.GetNbinsX())
        print('hx =', hx, hx.Integral(), hx.GetNbinsX())
        print('hout =', hout, hout.Integral(), hout.GetNbinsX())
        print('Rx =', Rx)
        print('hRx =', hRx, hRx.Integral(), hRx.GetNbinsX())

    fsc = hx.Integral() / hout.Integral()
    #hout.Scale(fsc)
    print('fsc =', fsc)
    showhist(hout)
    return hb,hx,hout, hRx

def showhist(h):
    print('showhist', h)
    nbins = h.GetNbinsX()
    print(h, nbins)
    for i in range(nbins):
        print(i, h.GetBinContent(i+1), '+/-', h.GetBinError(i+1))
    return


def unfold(hA, hx, hb, Icase=1, Imethod=1, tree=None):
    print('we are in unfold()')
    hRx = None
    if Imethod!=4:
        print('Dim(meas) =', hb.GetNbinsX(), hA.GetXaxis().GetNbins())
        print('Dim(true) =', hx.GetNbinsX(), hA.GetYaxis().GetNbins())
        response = RooUnfoldResponse(0, 0, hA)
        if Imethod == 0:
            unfold = RooUnfoldInvert(response, hb)
        elif Imethod == 1:
            unfold = RooUnfoldTUnfold(response, hb)
            unfold.FixTau(2e-3)
        elif Imethod == 2:
            unfold = RooUnfoldSvd(response, hb, 30, 10)
            #print('kterm =', unfold.GetKterm())
        elif Imethod == 3:
            unfold = RooUnfoldBayes(response, hb, 200, 0)
        hout = unfold.Hreco(2)
        Vx = unfold.Ereco(2)
        n = Vx.GetNrows()
        hRx = TH2F('hRx', '', n, 0, n, n, 0, n)
        for i in range(n):
            for j in range(n):
                vij = Vx(i, j)
                vii = Vx(i, i)
                vjj = Vx(j, j)
                rij = 0
                if vii*vjj>0:
                    rij = vij/math.sqrt(vii*vjj)
                hRx.SetBinContent(i+1,j+1, rij)
    elif Imethod == 4:
        hb, hx, hout, hRx = newunfold(Icase)
    else:
        return None
    return hb, hx, hout, hRx
def get_min_max(h):
    ymin = 9999
    ymax = -9999
    nbins = h.GetNbinsX()
    for i in range(nbins):
        y = h.GetBinContent(i+1)
        if y > ymax:
            ymax = y
        if y < ymin:
            ymin = y
    return ymin,ymax
def get_bw_min_max(h):
    ymin = 9999
    ymax = -9999
    nbins = h.GetNbinsX()
    for i in range(nbins):
        #y = h.GetBinContent(i+1)
        y = h.GetBinWidth(i+1)
        if y > ymax:
            ymax = y
        if y < ymin:
            ymin = y
    return ymin,ymax

def plot_Axb(hA, hx, hb, hout, hRx=None, uniform=1, plotname=''):
    if plotname == '':
        plotname = 'unfold'
    ymin,ymax = get_min_max(hx)
    Cs_A = TCanvas('Cs_A', '', 10, 10, 800, 600)
    Cs_A.SetRightMargin(0.17)
    hA.Draw('colz')
    hA.GetXaxis().SetTitle('Meas.')
    hA.GetYaxis().SetTitle('True')
    hA.GetZaxis().SetTitle('Prob.')
    Cs_A.SaveAs('Cs_A_'+plotname+'.png')
    Cs_A.SaveAs('Cs_A_'+plotname+'.pdf')
    if hRx != None:
        Cs_Rx = TCanvas('Cs_Rx', '', 10, 10, 800, 600)
        Cs_Rx.SetRightMargin(0.17)
        hRx.Draw('colz')
        hRx.GetXaxis().SetTitle('Unfolded')
        hRx.GetYaxis().SetTitle('Unfolded')
        hRx.GetZaxis().SetTitle('Correlation Coefficient')
        Cs_Rx.SaveAs('Cs_Rx_'+plotname+'.png')
        Cs_Rx.SaveAs('Cs_Rx_'+plotname+'.pdf')
    Cs = TCanvas('Cs', '', 10, 10, 800, 600)
    hx.Draw('hist')
    if hb.GetNbinsX() == hx.GetNbinsX() or 1:
        hb.Draw('hist,same')
    hout.Draw('PEX,same')
    hout.SetMarkerStyle(20)
    hout.SetMarkerColor(ROOT.kBlue)
    hout.SetLineColor(ROOT.kBlue)
    hb.SetLineStyle(ROOT.kDashed)
    hb.SetLineColor(ROOT.kRed)
    hx.GetXaxis().SetTitle('mass [GeV]')
    if uniform:
        hx.GetYaxis().SetTitle('Events / %.2f GeV' % hx.GetBinWidth(1))
    else:
        minbw, maxbw = get_bw_min_max(hx)
        hx.GetYaxis().SetTitle('Events / Width (%.1f-%.1f GeV)' % (minbw, maxbw))
    hx.GetYaxis().SetRangeUser(0, ymax*2.)
    leg = TLegend(0.6, 0.6, 0.95, 0.95)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hx, 'True', 'L')
    leg.AddEntry(hb, 'Meas.', 'L')
    leg.AddEntry(hout, 'Unfolded', 'PE')
    leg.Draw()
    Cs.SaveAs('Cs_xb_'+plotname+'.png')
    Cs.SaveAs('Cs_xb_'+plotname+'.pdf')
def single_unfold(hA, hx, hb, Icase=1, Imethod=1, tree=None):
    print('we are in single_unfold()')
    hb, hx, hout, hRx = unfold(hA, hx, hb, Icase, Imethod, tree)
    uniform = 1
    if Imethod == 4:
        uniform = 0
    plot_Axb(hA, hx, hb, hout, hRx, uniform, 'case'+str(Icase)+'_method'+str(Imethod))
    return
def main(Icase = 1, Imethod=1):
    filename = 'data_'+str(Icase)+'.root'
    if not os.path.isfile(filename):
        print('generate dataset ......................................................')
        generate_dataset(Icase)
    f = TFile(filename, 'read')
    tree = f.Get('nominal')
    nx = 30
    ny = 30
    if Icase == 1:
        nx = 60
        ny = nx
        if Imethod == 1:
            ny=100
    elif Icase == 2:
        nx = 50
        #nx = 10
        ny = nx
        if Imethod == 1:
            ny = 100
            #ny = 20
    elif Icase == 3:
        nx = 30
        ny = nx
        if Imethod == 1:
            ny = 60
    print('nx, ny =', nx,ny)
    hA, hx, hb = prepare_Axb(tree, nx, ny)
    if 0:
        hA, hx, hb, hxp, hbp = prepare_Axb(tree, nx, ny, 1)
        savefile = TFile('unfold_test.root', 'recreate')
        hA.Write('hA')
        hx.Write('hx')
        hb.Write('hy')
        hxp.Write('hxp')
        hbp.Write('hyp')
        savefile.Close()
    if 0:
        print('for debugging, hb = ..., hA = ...')
        nbins = hb.GetNbinsX()
        for i in range(nbins):
            print(i, hb.GetBinContent(i+1), hA.GetBinContent(1, i+1))
    single_unfold(hA, hx, hb, Icase, Imethod, tree)

for Icase in [2]:
    for Imethod in [4]:
        main(Icase, Imethod)