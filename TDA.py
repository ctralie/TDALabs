"""
Programmer: Chris Tralie (ctralie@alumni.princeton.edu)
Purpose: To wrap around Ripser to compute persistence diagrams
"""
import subprocess
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from TDAPlotting import *

def parseCocycle(s):
    s2 = "" + s
    for c in ["]", "[", ",", "{", "}"]:
        s2 = s2.replace(c, "")
    s2 = s2.replace(":", " ")
    cocycle = [int(c) for c in s2.split()]
    cocycle = np.array(cocycle)
    cocycle = np.reshape(cocycle, [len(cocycle)/3, 3])
    return cocycle

def doRipsFiltrationDM(D, maxHomDim, thresh = -1, coeff = 2, getCocycles = False):
    """
    Wrapper around Uli Bauer's Ripser code
    :param D: An NxN pairwise distance matrix
    :param maxHomDim: The dimension up to which to compute persistent homology
    :param thresh: Threshold up to which to add edges.  If not specified, add all
        edges up to the full clique
    :param coeff: A prime to use as the field coefficients for the PH computation
    :param getCocycles: True if cocycles should be computed and returned

    :return: PDs (array of all persistence diagrams from 0D up to maxHomDim).
        Each persistence diagram is a numpy array
        OR
        tuple (PDs, Cocycles) if returning cocycles
    """
    N = D.shape[0]
    #Step 1: Extract and output distance matrix
    fout = open("D.txt", "w")
    for i in range(0, N):
        for j in range(0, N):
            fout.write("%g "%D[i, j])
        if i < N-1:
            fout.write("\n")
    fout.close()

    #Step 2: Call ripser
    callThresh = 2*np.max(D)
    if thresh > 0:
        callThresh = thresh
    if getCocycles:
        proc = subprocess.Popen(["ripser/ripser-representatives", "--format", "distance", "--dim", "%i"%maxHomDim, "--threshold", "%g"%callThresh, "--modulus", "%i"%coeff, "D.txt"], stdout=subprocess.PIPE)
    elif coeff > 2:
        proc = subprocess.Popen(["ripser/ripser-coeff", "--dim", "%i"%maxHomDim, "--threshold", "%g"%callThresh, "--modulus", "%i"%coeff, "D.txt"], stdout=subprocess.PIPE)
    else:
        proc = subprocess.Popen(["ripser/ripser", "--dim", "%i"%maxHomDim, "--threshold", "%g"%callThresh, "D.txt"], stdout=subprocess.PIPE)
    #stdout = proc.communicate()[0]
    PDs = []
    AllCocycles = []
    while True:
        output=proc.stdout.readline()
        if (output == b'' or output == '') and proc.poll() is not None:
            break
        if output:
            s = output.strip()
            if output[0:4] == b"dist":
                continue
            elif output[0:4] == b"valu":
                continue
            elif output[0:4] == b"pers":
                if len(PDs) > 0:
                    PDs[-1] = np.array(PDs[-1])
                PDs.append([])
                AllCocycles.append([])
            else:
                if getCocycles:
                    s = s.split(": ")
                    if len(s) > 1:
                        [s, s1] = s
                        c = parseCocycle(s1)
                        AllCocycles[-1].append(c)
                    else:
                        s = s[0]
                s = s.replace(b"[", b"")
                s = s.replace(b"]", b"")
                s = s.replace(b"(", b"")
                s = s.replace(b")", b"")
                s = s.replace(b" ", b"")
                fields = s.split(b",")
                b = float(fields[0])
                d = -1
                if len(fields[1]) > 0:
                    d = float(fields[1])
                PDs[-1].append([b, d])
        rc = proc.poll()
    PDs[-1] = np.array(PDs[-1])
    if getCocycles:
        return (PDs, AllCocycles)
    return PDs

def getSSM(X):
    XSqr = np.sum(X**2, 1)
    D = XSqr[:, None] + XSqr[None, :] - 2*X.dot(X.T)
    D[D < 0] = 0 #Numerical precision
    D = np.sqrt(D)
    return D

def doRipsFiltration(X, maxHomDim, thresh = -1, coeff = 2, getCocycles = False):
    """
    Run ripser assuming Euclidean distance of a point cloud X
    :param X: An N x d dimensional point cloud
    :param maxHomDim: The dimension up to which to compute persistent homology
    :param thresh: Threshold up to which to add edges.  If not specified, add all
        edges up to the full clique
    :param coeff: A prime to use as the field coefficients for the PH computation
    :param getCocycles: True if cocycles should be computed and returned

    :return: PDs (array of all persistence diagrams from 0D up to maxHomDim).
        Each persistence diagram is a numpy array
        OR
        tuple (PDs, Cocycles) if returning cocycles
    """
    D = getSSM(X)
    return doRipsFiltrationDM(D, maxHomDim, thresh, coeff, getCocycles)

if __name__ == '__main__':
    np.random.seed(10)
    X = np.random.randn(200, 2)
    X = X/np.sqrt(np.sum(X**2, 1)[:, None])
    PDs = doRipsFiltration(X, 1, coeff = 3)
    plotDGM(PDs[1])
    plt.show()
