#Programmer: Chris Tralie
#Purpose: To wrap around Rann's pipeline to compute persistence diagrams and
#Dionysus for computing bottleneck distance
import subprocess
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from SparseEdgeList import *

def plotDGM(dgm, color = 'b', sz = 20, label = 'dgm'):
    if dgm.size == 0:
        return
    # Create Lists
    # set axis values
    axMin = np.min(dgm)
    axMax = np.max(dgm)
    axRange = axMax-axMin;
    # plot points
    plt.scatter(dgm[:, 0], dgm[:, 1], sz, color,label=label)
    plt.hold(True)
    # plot line
    plt.plot([axMin-axRange/5,axMax+axRange/5], [axMin-axRange/5, axMax+axRange/5],'k');
    # adjust axis
    #plt.axis([axMin-axRange/5,axMax+axRange/5, axMin-axRange/5, axMax+axRange/5])
    # add labels
    plt.xlabel('Time of Birth')
    plt.ylabel('Time of Death')

def plot2DGMs(P1, P2, l1 = 'Diagram 1', l2 = 'Diagram 2'):
    plotDGM(P1, 'r', 10, label = l1)
    plt.hold(True)
    plt.plot(P2[:, 0], P2[:, 1], 'bx', label = l2)
    plt.legend()
    plt.xlabel("Birth Time")
    plt.ylabel("Death Time")

def savePD(filename, I):
    if os.path.exists(filename):
        os.remove(filename)
    fout = open(filename, "w")
    for i in range(I.shape[0]):
        fout.write("%g %g"%(I[i, 0], I[i, 1]))
        if i < I.shape[0]-1:
            fout.write("\n")
    fout.close()

#Wrap around Dionysus's bottleneck distance after taking the log
def getInterleavingDist(PD1, PD2):
    savePD("PD1.txt", np.log(PD1))
    savePD("PD2.txt", np.log(PD2))
    proc = subprocess.Popen(["./bottleneck", "PD1.txt", "PD2.txt"], stdout=subprocess.PIPE)
    lnd = float(proc.stdout.readline())
    return np.exp(lnd) - 1.0 #Interleaving dist is 1 + eps

def getBottleneckDist(PD1, PD2):
    savePD("PD1.txt", PD1)
    savePD("PD2.txt", PD2)
    proc = subprocess.Popen(["./bottleneck", "PD1.txt", "PD2.txt"], stdout=subprocess.PIPE)
    return float(proc.stdout.readline())

def parsePDs(filename):
    PDs = {}
    fin = open(filename)
    for l in fin.readlines():
        fs = [float(s.rstrip()) for s in l.split()]
        dim = int(fs[0])
        if not dim in PDs:
            PDs[dim] = []
        if fs[-2] == fs[-1]:
            continue #Don't display classes which die instantly
        PDs[dim].append(fs[-2:])
    fin.close()
    ret = []
    count = 0
    for i in range(len(PDs)):
        ret.append(np.array(PDs[i]))
    return ret

def getPDs(I, J, D, N, m):
    if os.path.exists("temp.dimacs"):
        os.remove("temp.dimacs")
    writeResults(I, J, D, N, "temp.dimacs")
    if os.path.exists("temp.results"):
        os.remove("temp.results")
    subprocess.call(["./phatclique", "-i", "temp.dimacs", "-m", "%i"%m, "-o", "temp.results"])
    return parsePDs("temp.results")

def doRipsFiltration(X, maxHomDim, eps = 0):
    (I, J, D) = makeComplex(X, 0)
    PDs = getPDs(I, J, D, X.shape[0], maxHomDim+2)
    return PDs
    
if __name__ == '__main__':
    X = np.random.randn(100, 2)
    X = X/np.sqrt(np.sum(X**2, 1)[:, None])
    #plt.plot(X[:, 0], X[:, 1], '.')
    #plt.show()
    PDs = doRipsFiltration(X)
    plotDGM(PDs[1])
    plt.show()
