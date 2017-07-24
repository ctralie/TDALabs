"""
Programmer: Chris Tralie (ctralie@alumni.princeton.edu)
Purpose: Code to compute features on audio files, including
audio novelty
"""
import numpy as np
import numpy.linalg as linalg
import scipy
from scipy.io import wavfile
from scipy.io import savemat
from scipy.fftpack import dct
import matplotlib.pyplot as plt


def Specgram(X, W, H):
    """A function to compute the spectrogram of a signal
    :parm X: N x 1 Audio Signal
    :param W: Window Size
    :param H HopSize
    :returns: S, an N x NBins spectrogram array
    """
    Q = W/H
    if Q - np.floor(Q) > 0:
        print('Warning: Window size is not integer multiple of hop size\n')
    win = np.hamming(W)
    NWin = int(np.floor((len(X) - W)/float(H)) + 1)
    S = np.zeros((NWin, W))
    for i in range(NWin):
        x = X[i*H:i*H+W]
        S[i, :] = np.abs(np.fft.fft(win*x))
    #Second half of the spectrum is redundant for real signals
    if W % 2 == 0:
        #Even Case
        S = S[:, 0:int(W/2)]
    else:
        #Odd Case
        S = S[:, 0:int((W-1)/2)+1]
    return S

def getMelFilterbank( Fs, winSize, nbands, minfreq, maxfreq ):
    #Purpose: Return a mel-spaced triangle filterbank
    #Step 1: Warp to the mel-frequency scale
    melbounds = np.array([minfreq, maxfreq])
    melbounds = 1125*np.log(1 + melbounds/700.0)
    mel = np.linspace(melbounds[0], melbounds[1], nbands)
    binfreqs = 700*(np.exp(mel/1125.0) - 1)
    binbins = np.ceil(((winSize-1)/float(Fs))*binfreqs) #Ceil to the nearest bin
    binbins = np.array(binbins, dtype = np.int64)

    #Step 2: Create mel triangular filterbank
    melfbank = np.zeros((nbands, winSize))
    for i in range(nbands):
       thisbin = binbins[i]
       lbin = thisbin
       if i > 0:
           lbin = binbins[i-1]
       rbin = thisbin + (thisbin - lbin)
       if i < nbands-1:
           rbin = binbins[i+1]
       melfbank[i, lbin:thisbin+1] = np.linspace(0, 1, 1 + (thisbin - lbin))
       melfbank[i, thisbin:rbin+1] = np.linspace(1, 0, 1 + (rbin - thisbin))
    return melfbank


def getAudioNoveltyFn(x, Fs, winSize, hopSize):
    """
    Using techniques from
    Ellis, Daniel PW. "Beat tracking by dynamic programming." 
    Journal of New Music Research 36.1 (2007): 51-60.
    """
    
    #First compute mel-spaced STFT
    S = Specgram(x, winSize, hopSize)
    S = np.abs(S)
    M = getMelFilterbank(Fs, winSize, 40, 30, 8000)
    M = M[:, 0:S.shape[1]]
    X = M.dot(S.T)
    
    novFn = X[:, 1::] - X[:, 0:-1]
    novFn[novFn < 0] = 0
    novFn = np.sum(novFn, 0)
    return (S, novFn)


if __name__ == '__main__':
    Fs, X = scipy.io.wavfile.read("journey.wav")
    X = X/(2.0**15) #Audio is loaded in as 16 bit shorts.  Convert to float
    winSize = 512
    hopSize = 256
    (S, novFn) = getAudioNoveltyFn(X, Fs, winSize, hopSize)
    
    nsamples = 500
    novFn = novFn[0:nsamples]
    t = np.arange(nsamples)*hopSize/float(Fs)
    
    plt.subplot(211)
    plt.imshow(np.log(S.T), cmap = 'afmhot', aspect = 'auto')
    plt.title("Spectrogram")
    plt.axis('off')
    plt.subplot(212)
    plt.plot(t, novFn)
    plt.title("Audio Novelty Function")
    plt.xlabel("Time (Sec)")
    plt.xlim([0, np.max(t)])
    plt.show()
