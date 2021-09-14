#!/usr/bin/env python

##\file
# Manage and run simulation processes over multiple CPUs and parameters
#

import os, sys
from subprocess import Popen
from time import sleep

NCPU = 8
fpath = "../../../G4QSim/build/"
fnbase = "Xe136HPG0nbb"
cmdHead = ["../../src/pixelation", "-e", "1000"]

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

def make_cmd(dpn, dst):
    return cmdHead + ["-t", str(dpn[0]), "-l", str(dpn[0]), "-p", str(dpn[1]), "-S", str(dpn[2]), dst]

def main():
    procs = [None for i in xrange(NCPU)]

    dpnL_prep()
    i = 0; j = 0
    while i < len(dpnL):
        while True:
            if j>=NCPU:
                j = 0
                sleep(1)
            if procs[j]:
                if procs[j].poll() != None:
                    break;
            else:
                break;
            j = j + 1

        dpn = dpnL[i]
        dst = make_link(dpn)
        cmd = make_cmd(dpn, dst)
        print cmd
        procs[j] = Popen(cmd)
        i = i + 1

if __name__ == '__main__':
    main()

