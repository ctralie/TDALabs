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
    lam[lam < 0] = 0
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
    
    filenames = []
    for a in ACTIVITIES:
        for person in range(1, 26): #Subject
            for t in range(1, 5): #Activity type
                filename = "KTH/%s/person%.2i_%s_d%i_uncomp.avi"%(a, person, a, t)
                if os.path.exists(filename):
                    filenames.append(filename)
    
    scores = []
    for i in range(len(filenames)):
        print("\n\n\n\n\n\n\n\n\n\ni = ", i)
        filename = filenames[i]
        print("Doing %s..."%filename)
        (XOrig, FrameDims) = loadVideo(filename)
        X = getPCAVideo(XOrig)
        [X, validIdx] = getTimeDerivative(X, 10)
        
        #Setup video blocks
        BlockLen = 160
        BlockHop = 80
        win = 20
        dim = 20
        idxs = []
        N = X.shape[0]
        NBlocks = int(np.ceil(1 + (N - BlockLen)/BlockHop))
        print("NBlocks = ", NBlocks)
        for k in range(NBlocks):
            thisidxs = np.arange(k*BlockHop, k*BlockHop+BlockLen, dtype=np.int64)
            thisidxs = thisidxs[thisidxs < N]
            idxs.append(thisidxs)
        
        wins = np.arange(2, 50)
        res = np.zeros(NBlocks)
        
        #Get sliding window video in blocks
        maxXS = []
        maxPD = []
        maxP = 0.0
        maxj = 0
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
                res[j] = np.max(PDs[1][:, 1] - PDs[1][:, 0])
                if res[j] > maxP:
                    maxP = res[j]
                    maxXS = np.array(XS)
                    maxPD = np.array(PDs[1])
                    maxj = j
        #Output the self-similarity matrix, the persistence diagram, and
        #PCA of the sliding window for the block with the maximum persistence 
        #score, as well as the video snippet that gave rise to that block
        
        #Self-similarity matrix
        XSqr = np.sum(maxXS**2, 1).flatten()
        D = XSqr[:, None] + XSqr[None, :] - 2*maxXS.dot(maxXS.T)
        D[D < 0] = 0
        D = np.sqrt(D)
        
        #PCA
        pca = PCA(n_components = 2)
        Y = pca.fit_transform(maxXS)
        
        #Video clip
        clip = XOrig[idxs[maxj], :]
        saveVideo(clip, FrameDims, "VideoResults/%i.ogg"%i)
        
        plt.clf()
        plt.imshow(D)
        plt.savefig("VideoResults/D%i.png"%i, bbox_inches='tight')
        
        plt.clf()
        plotDGM(maxPD)
        plt.savefig("VideoResults/PD%i.svg"%i, bbox_inches='tight')
        
        plt.clf()
        ax = plt.subplot(111)
        ax.set_title("PCA of Sliding Window Embedding")
        ax.scatter(Y[:, 0], Y[:, 1])
        ax.set_aspect('equal', 'datalim')
        plt.savefig("VideoResults/Y%i.svg"%i, bbox_inches='tight')
        
        scores.append(maxP)
    
    scores = np.array(scores)
    
    #Output results in HTML format in descending order of maximum persistence
    fout = open("VideoResults/index.html", "w")
    fout.write("<html><body><table border = '1'>")
    idx = np.argsort(-scores)
    count = 1
    for i in idx:
        fout.write("<tr><td><h2>%i</h2>%s<BR><BR>Maximum Persistence = <BR><b>%g</b></td>"%(count, filenames[i], scores[i]))
        fout.write("<td><video controls><source src=\"%i.ogg\" type=\"video/ogg\">Your browser does not support the video tag.</video>"%i)
        fout.write("<td><img src = \"PD%i.svg\"></td>"%i)
        fout.write("<td><img src = \"D%i.png\"></td>"%i)
        fout.write("<td><img src = \"Y%i.svg\"></td>"%i)
        fout.write("</tr>\n")
        count += 1
    fout.write("</table></body></html>")
    fout.close()
