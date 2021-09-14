#!/usr/bin/env python

import math
import sys

def array2roc(ary, n):
    """
    Calculate a ROC curve using n bins between probability 0.0 and 1.0.
    """
    dp = 1.0 / float(n)
    s = [0 for i in xrange(n)]
    b = [0 for i in xrange(n)]
    ssum = 0
    bsum = 0
    for [p, t] in ary:
        i = int(math.floor(p / dp))
        if t == 0: # 0 tags signal, 1 tags background
            s[i] += 1
            ssum += 1
        else:
            b[i] += 1
            bsum += 1
    s[0] /= float(ssum)
    b[0] /= float(bsum)
    for i in xrange(1,n):
        s[i] = s[i-1] + s[i]/float(ssum)
        b[i] = b[i-1] + b[i]/float(bsum)
    return (s, b)

def file2roc(fname, n):
    """
    Read a file into an array and calculate a ROC curve
    """
    ary = []
    for line in open(fname):
        elems = line.split()
        ary.append([float(elems[0]), int(elems[1])])
    return array2roc(ary, int(n))

if __name__ == "__main__":
    (s, b) = file2roc(sys.argv[1], sys.argv[2])
    print "# %s" % sys.argv[1]
    for i in xrange(len(s)):
        print "%10g %10g" % (s[i], b[i])
    print "\n"
