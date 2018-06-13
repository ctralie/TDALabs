"""
Author: Chris Tralie
Description: Contains methods to plot and compare persistence diagrams
              Comparison algorithms include grabbing/sorting, persistence landscapes,
              the "multiscale heat kernel" (CVPR 2015), and "persistence images" (Adams et al.)
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.misc #Used for downsampling rasterized images avoiding aliasing
import time #For timing kernel comparison
import sklearn.metrics.pairwise
import scipy.stats
from ripser import ripser, plot_dgms

##############################################################################
##########                  Plotting Functions                      ##########
##############################################################################

def plotBottleneckMatching(I1, I2, matchidx, D, labels = ['dgm1', 'dgm2']):
    plot_dgms([I1, I2], labels = labels)
    cp = np.cos(np.pi/4)
    sp = np.sin(np.pi/4)
    R = np.array([[cp, -sp], [sp, cp]])
    if I1.size == 0:
        I1 = np.array([[0, 0]])
    if I2.size == 0:
        I2 = np.array([[0, 0]])
    I1Rot = I1.dot(R)
    I2Rot = I2.dot(R)
    dists = [D[i, j] for (i, j) in matchidx]
    (i, j) = matchidx[np.argmax(dists)]
    if i >= I1.shape[0] and j >= I2.shape[0]:
        return
    if i >= I1.shape[0]:
        diagElem = np.array([I2Rot[j, 0], 0])
        diagElem = diagElem.dot(R.T)
        plt.plot([I2[j, 0], diagElem[0]], [I2[j, 1], diagElem[1]], 'g')
    elif j >= I2.shape[0]:
        diagElem = np.array([I1Rot[i, 0], 0])
        diagElem = diagElem.dot(R.T)
        plt.plot([I1[i, 0], diagElem[0]], [I1[i, 1], diagElem[1]], 'g')
    else:
        plt.plot([I1[i, 0], I2[j, 0]], [I1[i, 1], I2[j, 1]], 'g')


def plotWassersteinMatching(I1, I2, matchidx, labels = ['dgm1', 'dgm2']):
    plot_dgms([I1, I2], labels = labels)
    cp = np.cos(np.pi/4)
    sp = np.sin(np.pi/4)
    R = np.array([[cp, -sp], [sp, cp]])
    if I1.size == 0:
        I1 = np.array([[0, 0]])
    if I2.size == 0:
        I2 = np.array([[0, 0]])
    I1Rot = I1.dot(R)
    I2Rot = I2.dot(R)
    for index in matchidx:
        (i, j) = index
        if i >= I1.shape[0] and j >= I2.shape[0]:
            continue
        if i >= I1.shape[0]:
            diagElem = np.array([I2Rot[j, 0], 0])
            diagElem = diagElem.dot(R.T)
            plt.plot([I2[j, 0], diagElem[0]], [I2[j, 1], diagElem[1]], 'g')
        elif j >= I2.shape[0]:
            diagElem = np.array([I1Rot[i, 0], 0])
            diagElem = diagElem.dot(R.T)
            plt.plot([I1[i, 0], diagElem[0]], [I1[i, 1], diagElem[1]], 'g')
        else:
            plt.plot([I1[i, 0], I2[j, 0]], [I1[i, 1], I2[j, 1]], 'g')


##############################################################################
##########            Diagram Comparison Functions                  ##########
##############################################################################

def getWassersteinDist(S, T):
    """
    Perform the Wasserstein distance matching between persistence diagrams.
    Assumes first two columns of S and T are the coordinates of the persistence
    points, but allows for other coordinate columns (which are ignored in
    diagonal matching)
    :param S: Mx(>=2) array of birth/death pairs for PD 1
    :param T: Nx(>=2) array of birth/death paris for PD 2
    :returns (tuples of matched indices, total cost, (N+M)x(N+M) cross-similarity)
    """
    import hungarian #Requires having compiled the library

    # Step 1: Compute CSM between S and T, including points on diagonal
    N = S.shape[0]
    M = T.shape[0]
    #Handle the cases where there are no points in the diagrams
    if N == 0:
        S = np.array([[0, 0]])
        N = 1
    if M == 0:
        T = np.array([[0, 0]])
        M = 1
    DUL = sklearn.metrics.pairwise.pairwise_distances(S, T)

    #Put diagonal elements into the matrix
    #Rotate the diagrams to make it easy to find the straight line
    #distance to the diagonal
    cp = np.cos(np.pi/4)
    sp = np.sin(np.pi/4)
    R = np.array([[cp, -sp], [sp, cp]])
    S = S[:, 0:2].dot(R)
    T = T[:, 0:2].dot(R)
    D = np.zeros((N+M, N+M))
    D[0:N, 0:M] = DUL
    UR = np.max(D)*np.ones((N, N))
    np.fill_diagonal(UR, S[:, 1])
    D[0:N, M:M+N] = UR
    UL = np.max(D)*np.ones((M, M))
    np.fill_diagonal(UL, T[:, 1])
    D[N:M+N, 0:M] = UL
    D = D.tolist()

    # Step 2: Run the hungarian algorithm
    matchidx = hungarian.lap(D)[0]
    matchidx = [(i, matchidx[i]) for i in range(len(matchidx))]
    matchdist = 0
    for pair in matchidx:
        (i, j) = pair
        matchdist += D[i][j]

    return (matchidx, matchdist, D)

def getBottleneckDist(S, T):
    """
    Perform the Bottleneck distance matching between persistence diagrams.
    Assumes first two columns of S and T are the coordinates of the persistence
    points, but allows for other coordinate columns (which are ignored in
    diagonal matching)
    :param S: Mx(>=2) array of birth/death pairs for PD 1
    :param T: Nx(>=2) array of birth/death paris for PD 2
    :returns (tuples of matched indices, total cost, (N+M)x(N+M) cross-similarity)
    """
    from bisect import bisect_left
    from hopcroftkarp import HopcroftKarp

    N = S.shape[0]
    M = T.shape[0]
    # Step 1: Compute CSM between S and T, including points on diagonal
    # L Infinity distance
    Sb, Sd = S[:, 0], S[:, 1]
    Tb, Td = T[:, 0], T[:, 1]
    D1 = np.abs(Sb[:, None] - Tb[None, :])
    D2 = np.abs(Sd[:, None] - Td[None, :])
    DUL = np.maximum(D1, D2)
    # Put diagonal elements into the matrix, being mindful that Linfinity
    # balls meet the diagonal line at a diamond vertex
    D = np.zeros((N+M, N+M))
    D[0:N, 0:M] = DUL
    UR = np.max(D)*np.ones((N, N))
    np.fill_diagonal(UR, 0.5*(S[:, 1]-S[:, 0]))
    D[0:N, M::] = UR
    UL = np.max(D)*np.ones((M, M))
    np.fill_diagonal(UL, 0.5*(T[:, 1]-T[:, 0]))
    D[N::, 0:M] = UL

    # Step 2: Perform a binary search + Hopcroft Karp to find the
    # bottleneck distance
    N = D.shape[0]
    ds = np.unique(D.flatten())
    ds = ds[ds > 0]
    ds = np.sort(ds)
    bdist = ds[-1]
    matching = {}
    while len(ds) >= 1:
        idx = 0
        if len(ds) > 1:
            idx = bisect_left(range(ds.size), int(ds.size/2))
        d = ds[idx]
        graph = {}
        for i in range(N):
            graph['%s'%i] = {j for j in range(N) if D[i, j] <= d}
        res = HopcroftKarp(graph).maximum_matching()
        if len(res) == 2*N and d < bdist:
            bdist = d
            matching = res
            ds = ds[0:idx]
        else:
            ds = ds[idx+1::]
    matchidx = [(i, matching['%i'%i]) for i in range(N)]
    return (matchidx, bdist, D)


def sortAndGrab(dgm, NBars = 10, BirthTimes = False):
    """
    Do sorting and grabbing with the option to include birth times
    Zeropadding is also taken into consideration
    """
    dgmNP = np.array(dgm)
    if dgmNP.size == 0:
        if BirthTimes:
            ret = np.zeros(NBars*2)
        else:
            ret = np.zeros(NBars)
        return ret
    #Indices for reverse sort
    idx = np.argsort(-(dgmNP[:, 1] - dgmNP[:, 0])).flatten()
    ret = dgmNP[idx, 1] - dgmNP[idx, 0]
    ret = ret[0:min(NBars, len(ret))].flatten()
    if len(ret) < NBars:
        ret = np.append(ret, np.zeros(NBars - len(ret)))
    if BirthTimes:
        bt = dgmNP[idx, 0].flatten()
        bt = bt[0:min(NBars, len(bt))].flatten()
        if len(bt) < NBars:
            bt = np.append(bt, np.zeros(NBars - len(bt)))
        ret = np.append(ret, bt)
    return ret


def getHeatRasterized(dgm, sigma, xrange, yrange, UpFac = 10):
    """
    Get a discretized verison of the solution of the heat flow equation
    described in the CVPR 2015 paper
    """
    I = np.array(dgm)
    if I.size == 0:
        return np.zeros((yrange.size, xrange.size))
    NX = xrange.size
    NY = yrange.size
    #Rasterize on a finer grid and downsample
    NXFine = UpFac*NX
    NYFine = UpFac*NY
    xrangeup = np.linspace(xrange[0], xrange[-1], NXFine)
    yrangeup = np.linspace(yrange[0], yrange[-1], NYFine)
    X, Y = np.meshgrid(xrangeup, yrangeup)
    u = np.zeros(X.shape)
    for ii in range(I.shape[0]):
        u = u + np.exp(-( (X - I[ii, 0])**2 + (Y - I[ii, 1])**2 )/(4*sigma))
        #Now subtract mirror diagonal
        u = u - np.exp(-( (X - I[ii, 1])**2 + (Y - I[ii, 0])**2 )/(4*sigma))
    u = (1.0/(4*np.pi*sigma))*u
    u = scipy.misc.imresize(u, (NY, NX))
    return u


def evalHeatKernel(dgm1, dgm2, sigma):
    """
    Evaluate the continuous heat-based kernel between dgm1 and dgm2 (more correct
    than L2 on the discretized verison above but may be slower because can't exploit
    fast matrix multiplication when evaluating many, many kernels)
    """
    kSigma = 0
    I1 = np.array(dgm1)
    I2 = np.array(dgm2)
    for i in range(I1.shape[0]):
        p = I1[i, 0:2]
        for j in range(I2.shape[0]):
            q = I2[j, 0:2]
            qc = I2[j, 1::-1]
            kSigma += np.exp(-(np.sum((p-q)**2))/(8*sigma)) - np.exp(-(np.sum((p-qc)**2))/(8*sigma))
    return kSigma / (8*np.pi*sigma)

def evalHeatDistance(dgm1, dgm2, sigma):
    """
    Return the pseudo-metric between two diagrams based on the continuous
    heat kernel
    """
    return np.sqrt(evalHeatKernel(dgm1, dgm1, sigma) + evalHeatKernel(dgm2, dgm2, sigma) - 2*evalHeatKernel(dgm1, dgm2, sigma))

def getPersistenceImage(dgm, plims, res, weightfn = lambda b, l: l, psigma = None):
    """
    Return a persistence image (Adams et al.)
    :param dgm: Nx2 array holding persistence diagram
    :param plims: An array [birthleft, birthright, lifebottom, lifetop] \
        limits of the actual grid will be rounded based on res
    :param res: Width of each pixel
    :param weightfn(b, l): A weight function as a function of birth time\
        and life time
    :param psigma: Standard deviation of each Gaussian.  By default\
        None, which indicates it should be res/2.0
    """
    #Convert to birth time/lifetime
    I = np.array(dgm)
    I[:, 1] = I[:, 1] - I[:, 0]
    
    #Create grid
    lims = np.array([np.floor(plims[0]/res), np.ceil(plims[1]/res), np.floor(plims[2]/res), np.ceil(plims[3]/res)])
    xr = np.arange(int(lims[0]), int(lims[1])+2)*res
    yr = np.arange(int(lims[2]), int(lims[3])+2)*res
    sigma = res/2.0
    if psigma:
        sigma = psigma        
            
    #Add each integrated Gaussian
    PI = np.zeros((len(yr)-1, len(xr)-1))
    for i in range(I.shape[0]):
        [x, y] = I[i, :]
        w = weightfn(x, y)
        if w == 0:
            continue
        #CDF of 2D isotropic Gaussian is separable
        xcdf = scipy.stats.norm.cdf((xr - x)/sigma)
        ycdf = scipy.stats.norm.cdf((yr - y)/sigma)
        X = ycdf[:, None]*xcdf[None, :]
        #Integral image
        PI += weightfn(x, y)*(X[1::, 1::] - X[0:-1, 1::] - X[1::, 0:-1] + X[0:-1, 0:-1])
    return {'PI':PI, 'xr':xr[0:-1], 'yr':yr[0:-1]}

def testBottleneckWassersteinNoisyCircle():
    N = 400
    np.random.seed(N)
    t = np.linspace(0, 2*np.pi, N+1)[0:N]
    X = np.zeros((N, 2))
    X[:, 0] = np.cos(t)
    X[:, 1] = np.sin(t)
    I1 = ripser(X)['dgms'][1]
    X2 = X + 0.1*np.random.randn(N, 2)
    I2 = ripser(X2)['dgms'][1]
    tic = time.time()
    (matchidxb, bdist, bD) = getBottleneckDist(I1, I2)
    btime = time.time() - tic 
    tic = time.time()
    (matchidxw, wdist, wD) = getWassersteinDist(I1, I2)
    wtime = time.time() - tic
    print("Elapsed Time Bottleneck: %.3g\nElapsed Time Wasserstein: %.3g"%(btime, wtime))
    
    plt.figure(figsize=(12, 6))
    plt.subplot(121)
    plotBottleneckMatching(I1, I2, matchidxb, bD)
    plt.title("Bottleneck Dist: %.3g"%bdist)
    plt.subplot(122)
    plotWassersteinMatching(I1, I2, matchidxw)
    plt.title("Wasserstein Dist: %.3g"%wdist)
    plt.show()

def writePD(I, filename):
    fout = open(filename, "w")
    for i in range(I.shape[0]):
        fout.write("%g %g"%(I[i, 0], I[i, 1]))
        if i < I.shape[0]-1:
            fout.write("\n")
    fout.close()

def testBottleneckVsHera(NTrials = 10):
    import subprocess
    for trial in range(NTrials):
        x = np.random.randn(100, 2)
        y = np.random.randn(100, 2)
        I1 = ripser(x)['dgms'][1]
        I2 = ripser(y)['dgms'][1]
        (matchidxb, bdist, D) = getBottleneckDist(I1, I2)
        print(bdist)
        writePD(I1, "PD1.txt")
        writePD(I2, "PD2.txt")
        subprocess.call(["./bottleneck_dist", "PD1.txt", "PD2.txt"])
        print("\n")


if __name__ == '__main__':
    #testBottleneckVsHera()
    testBottleneckWassersteinNoisyCircle()
