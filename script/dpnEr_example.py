#!/usr/bin/env python

##\file
# Extract parameters out of Diffusion-pitch-noise-Energy resolution simulation
#

import os, sys, glob,re
from subprocess import Popen
from time import sleep

from ROOT import gROOT, gSystem, gStyle, gDirectory, TFile, TCanvas, TPad, TCut, TGraph, TGraph2D, TGraphErrors, TH1D, TH2D, TF1, TF2, TProfile, TProfile2D

NCPU = 8
fpath = "../../../G4QSim/build/"
fnbase = "Xe136HPG0nbb"
cmdHead = ["../../src/pixelation", "-e", "1000"]

def getParams(fname):
    params = []
    
    tf = TFile(fname, "read")
    p1 = tf.Get("p1")

    nIonTotH = TH1D("nIonTotH", "Total number of ions", 50, 94000, 95500)
    f0 = TF1("f0", "gaus", 94000, 95500)
    f0.SetParameters(500, 94700, 131)

    nPixTotH = TH1D("nPixTotH", "Total number of pixel hits", 800, 0, 1600)
    f1 = TF1("f1", "gaus", 0, 1600)
    f1.SetParameters(100, 200, 100)

    sigTotH  = TH1D("sigTotH",  "Total signal with noise", 50, 92000, 97000)
    f2 = TF1("f2", "gaus", 92000, 97000)
    f2.SetParameters(500, 94500, 200)

    maxH     = TH1D("maxH",     "Maximum signal on a pixel", 50, 0, 10000)
# skewed Gaussian: skew-normal distribution
    f3 = TF1("f3", "[0]*exp(-0.5*(x-[1])/[2]*(x-[1])/[2])*(1+TMath::Erf([3]/sqrt(2.0)*(x-[1])/[2]))", 0, 10000)
    f3.SetParameters(20, 3000, 500, 4.0)

    p1.Draw("nIonTot>>nIonTotH", "", "goff")
    nIonTotH.Fit("f0", "MRN E0")
    params.append(f0.GetParameter(1))
    params.append(f0.GetParameter(2))

    p1.Draw("nPixTot>>nPixTotH", "", "goff")
    nPixTotH.Fit("f1", "MRN E0")
    params.append(f1.GetParameter(1))
    params.append(f1.GetParameter(2))

    p1.Draw("sigTot>>sigTotH", "", "goff")
    sigTotH.Fit("f2", "MRN E0")
    params.append(f2.GetParameter(1))
    params.append(f2.GetParameter(2))

    p1.Draw("lzH.ProjectionX().GetMaximum()>>maxH", "", "goff")
    f3.SetParameters(maxH.GetMaximum(),
                     maxH.GetMean(),
                     maxH.GetRMS(),
                     4.0)
    maxH.Fit("f3", "MRN E0")
    params.append(f3.GetParameter(1))
    params.append(f3.GetParameter(2))
    params.append(f3.GetParameter(3))

    tf.Close()
    return params

dpnL = []

def dpnL_prep():
    for D in xrange(1, 11):
        for p in xrange(1, 11):
            for n in xrange(10, 31, 10):
                dpnL.append([D, p, n])

def make_link(dpn):
    src = fpath + fnbase + ".root"
    dst = fnbase + "-%02d-%02d-%02d.root" % tuple(dpn)
    try:
        os.symlink(src, dst)
    except:
        pass
    return dst

def fname2dpn(fname):
    m = re.match(r'.*-([0-9]+)-([0-9]+)-([0-9]+)', fname)
    return [float(i) for i in m.groups()]

def main():
    of = open(sys.argv[2], "w")
    of.write("# Dt pitch noise nIonTot sigma nPixTot sigma sigTot sigma maxH sigma alpha\n")
    for fname in glob.glob(sys.argv[1]):
        dpn = fname2dpn(fname)
        for a in dpn:
            of.write(" " + str(a))
        for a in getParams(fname):
            of.write(" " + str(a))
        of.write("\n")
        of.flush()
    of.close()

if __name__ == '__main__':
    main()
