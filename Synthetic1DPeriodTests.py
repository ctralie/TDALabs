from VideoTools import *
from TDA import *
import os
import numpy as np
import scipy.io as sio
import scipy.interpolate as interp
from sklearn.decomposition import PCA

def getSlidingWindow(x, dim, Tau, dT):
    NWindows = int(np.floor((N-dim*Tau)/dT))
    X = np.zeros((NWindows, dim))
    idx = np.arange(len(x))
    for i in range(NWindows):
        idxx = dT*i + Tau*np.arange(dim)
        start = int(np.floor(idxx[0]))
        end = int(np.ceil(idxx[-1]))
        X[i, :] = interp.spline(idx[start:end+1], x[start:end+1], idxx)
    return X

def getPulseTrain(NSamples, TMin, TMax, AmpMin, AmpMax):
    x = np.zeros(NSamples)
    x[0] = 1
    i = 0
    while i < NSamples:
        i += TMin + int(np.round(np.random.randn()*(TMax-TMin)))
        if i >= NSamples:
            break
        x[i] = AmpMin + (AmpMax-AmpMin)*np.random.randn()  
    return x

def convolveAndAddNoise(x, gaussSigma, noiseSigma):
    gaussSigma = int(np.round(gaussSigma*3))
    g = np.exp(-(np.arange(-gaussSigma, gaussSigma+1, dtype=np.float64))**2/(2*gaussSigma**2))
    x = np.convolve(x, g, 'same')
    x = x + noiseSigma*np.random.randn(len(x))
    return x

def getSyntheticPulseTrain(NSamples, T, noiseSigma, gaussSigma):
    x = np.zeros(NSamples)
    x[0::T] = 1
    x = convolveAndAddNoise(x, gaussSigma, noiseSigma)
    return x

if __name__ == '__main__':
    T = 40 #The period in number of samples
    NPeriods = 3 #How many periods to go through
    N = T*NPeriods #The total number of samples
    t = np.linspace(0, 2*np.pi*NPeriods, N+1)[0:N] #Sampling indices in time
    x = np.cos(t) #The final signal
    
    
    wins = np.linspace(2, 50, 100)
    dim = 10
    res = np.zeros(len(wins))
    for i in range(len(wins)):
        Tau = wins[i]/float(dim-1)
        dT = (N-dim*Tau)/float(N)
        XS = getSlidingWindow(x, dim, Tau, dT)

        #Mean-center and normalize sliding window
        #XS = XS - np.mean(XS, 1)[:, None]
        #XS = XS/np.sqrt(np.sum(XS**2, 1))[:, None]
        
        #XS = XS - np.mean(XS, 0)[None, :]
        #XS = XS/np.sqrt(np.sum(XS**2, 1))[:, None]
        #XS = XS/np.sqrt(XS.shape[1])
        
        PDs = doRipsFiltration(XS, 1)

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(221)
        pca = PCA(n_components = 2)
        Y = pca.fit_transform(XS)
        ax.set_title("PCA of Sliding Window Embedding")
        ax.scatter(Y[:, 0], Y[:, 1])
        ax.set_aspect('equal', 'datalim')
        plt.subplot(223)
        W = int(np.round(wins[i]))
        plt.plot(np.arange(len(x)), x, 'b')
        plt.hold(True)
        plt.plot(np.arange(W), x[0:W], 'r')
        plt.plot([wins[i], wins[i]], [np.min(x), np.max(x)])
        plt.title("Signal Chunk")
        plt.subplot(224)
        plt.plot(wins, res)
        plt.xlabel("Window Size")
        plt.ylabel("Maximum Persistence")
        
        if len(PDs) == 2:
            if PDs[1].size > 0:
                ax2 = fig.add_subplot(222)
                plotDGM(PDs[1])
                plt.title("Window = %g"%wins[i])
        plt.savefig("PD%i.png"%i, bbox_inches='tight')
        if len(PDs) < 2:
            continue
        if PDs[1].size > 0:
            res[i] = np.max(PDs[1][:, 1] - PDs[1][:, 0])
    
    sio.savemat("res.mat", {"res":res})
