from VideoTools import *
from TDA import *
import os
import numpy as np
import scipy.io as sio
import scipy.interpolate
import scipy.signal
from sklearn.decomposition import PCA

def getPCAVideo(I):
    ICov = I.dot(I.T)
    [lam, V] = linalg.eigh(ICov)
    V = V*np.sqrt(lam[None, :])
    return V

def getSlidingWindowVideo(I, dim, Tau, dT):
    N = I.shape[0] #Number of frames
    P = I.shape[1] #Number of pixels (possibly after PCA)
    pix = np.arange(P)
    NWindows = int(np.floor((N-dim*Tau)/dT))
    X = np.zeros((NWindows, dim*P))
    idx = np.arange(N)
    for i in range(NWindows):
        idxx = dT*i + Tau*np.arange(dim)
        start = int(np.floor(idxx[0]))
        end = int(np.ceil(idxx[-1]))
        f = scipy.interpolate.interp2d(pix, idx[start:end+1], I[idx[start:end+1], :], kind='linear')
        X[i, :] = f(pix, idxx).flatten()
    return X

def getTimeDerivative(I, Win):
    dw = np.floor(Win/2)
    t = np.arange(-dw, dw+1)
    sigma = 0.4*dw
    xgaussf = t*np.exp(-t**2  / (2*sigma**2))
    #Normalize by L1 norm to control for length of window
    xgaussf = xgaussf/np.sum(np.abs(xgaussf)) 
    xgaussf = xgaussf[:, None]
    IRet = scipy.signal.convolve2d(I, xgaussf, 'valid')
    validIdx = np.arange(dw, I.shape[0]-dw, dtype='int64')
    return [IRet, validIdx]

ACTIVITIES = ["boxing", "handclapping", "handwaving", "jogging", "running", "walking"]

if __name__ == '__main__':
    #Load Video
    #filename = "KTH/handwaving/person01_handwaving_d1_uncomp.avi"
    filename = "KTH/boxing/person01_boxing_d1_uncomp.avi"
    (X, FrameDims) = loadVideo(filename)
    print("Doing PCA for dimensionality reduction...")
    X = getPCAVideo(X)
    print("Finished PCA")
    [X, validIdx] = getTimeDerivative(X, 10)
    
    XSqr = np.sum(X**2, 1).flatten()
    D = XSqr[:, None] + XSqr[None, :] - 2*X.dot(X.T)
    D[D < 0] = 0
    D = np.sqrt(D)
    
    #Setup video blocks
    BlockLen = 160
    BlockHop = 80
    idxs = []
    N = X.shape[0]
    NBlocks = int(np.ceil(1 + (N - BlockLen)/BlockHop))
    print("NBlocks = ", NBlocks)
    for i in range(NBlocks):
        thisidxs = np.arange(i*BlockHop, i*BlockHop+BlockLen, dtype=np.int64)
        thisidxs = thisidxs[thisidxs < N]
        idxs.append(thisidxs)
    
    plt.imshow(D)
    plt.show()
    
    wins = np.arange(2, 50)
    dim = 40
    res = np.zeros((len(wins), NBlocks))
    for i in range(len(wins)):
        win = wins[i]
        
        #Get sliding window video in blocks
        for j in range(len(idxs)):
            idx = idxs[j]
            Tau = win/float(dim-1)
            dT = (len(idx)-dim*Tau)/float(len(idx))
            XS = getSlidingWindowVideo(X[idx, :], dim, Tau, dT)

            #Mean-center and normalize sliding window
            XS = XS - np.mean(XS, 1)[:, None]
            XS = XS/np.sqrt(np.sum(XS**2, 1))[:, None]
            
            PDs = doRipsFiltration(XS, 1)
            if len(PDs) < 2:
                continue
            if PDs[1].size > 0:
                res[i, j] = np.max(PDs[1][:, 1] - PDs[1][:, 0])
            if j == 0:
                XSqr = np.sum(XS**2, 1).flatten()
                D = XSqr[:, None] + XSqr[None, :] - 2*XS.dot(XS.T)
                D[D < 0] = 0
                D = np.sqrt(D)
            
                fig = plt.figure(figsize=(12, 12))
                plt.subplot(223)
                plt.imshow(D)
                plt.title("Self-Similarity Matrix")
                ax = fig.add_subplot(221)
                pca = PCA(n_components = 2)
                Y = pca.fit_transform(XS)
                ax.set_title("PCA of Sliding Window Embedding")
                ax.scatter(Y[:, 0], Y[:, 1])
                ax.set_aspect('equal', 'datalim')
                plt.subplot(224)
                plt.plot(wins, res[:, 0])
                plt.xlabel("Window Size")
                plt.ylabel("Maximum Persistence")
                if len(PDs) == 2:
                    if PDs[1].size > 0:
                        ax2 = fig.add_subplot(222)
                        plotDGM(PDs[1])
                        plt.title("Window = %g"%wins[i])
                
                plt.savefig("PD%i.png"%i, bbox_inches='tight')
    
    sio.savemat("res.mat", {"res":res})
