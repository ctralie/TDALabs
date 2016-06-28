#Programmer: Chris Tralie
#Purpose: Return a warped edge list based on Cavanna/Jahanseir/Sheehy's approach
#to approximate rips filtrations over some metric space, using a greedy permutation.
#The goal is to return many fewer edges than all O(N^2)

from sys import argv, exit
import numpy as np

#Purpose: Naive O(N^2) algorithm to do the greedy permutation
#Inputs: D (NxN distance matrix for points)
#Returns: (permutation (N-length array of indices), 
#lambdas (N-length array of insertion radii))
def getGreedyPerm(D):
    N = D.shape[0]
    #By default, takes the first point in the list to be the
    #first point in the permutation, but could be random
    perm = np.zeros(N, dtype=np.int64)
    lambdas = np.zeros(N)
    ds = D[0, :]
    for i in range(1, N):
        idx = np.argmax(ds)
        perm[i] = idx
        lambdas[i] = ds[idx]
        ds = np.minimum(ds, D[idx, :])
    return (perm, lambdas)

#Purpose: Do a binary search to find an epsilon corresponding to
#NEdges included edges
def getEps(lambdaso, D, NEdges, eps1 = 0, eps2 = 1):
    N = D.shape[0]
    b1 = ((eps1**2+3*eps1+2)/eps1)*lambdaso
    b2 = ((eps2**2+3*eps2+2)/eps2)*lambdaso
    f1 = np.sum(D > b1)
    f2 = np.sum(D > b2)
    

#Purpose: To return the sparse edge list with the warped distances, sorted
#by weight
#Inputs: lambdaso (insertion radii for points), eps (epsilon approximation constant),
#D (NxN distance matrix, okay to modify because last time it's used)
def getSparseEdgeList(lambdaso, eps, D):
    N = D.shape[0]
    E0 = (1+eps)/eps
    E1 = (1+eps)**2/eps
    
    #Create initial sparse list candidates (Lemma 6)
    nBounds = ((eps**2+3*eps+2)/eps)*lambdaso #Search neighborhoods    
    D[D > nBounds[:, None]] = np.inf #Set all distances outside of search neighborhood to infinity
    [I, J] = np.meshgrid(np.arange(N), np.arange(N))
    idx = I > J
    I = I[(D < np.inf)*(idx == 1)]
    J = J[(D < np.inf)*(idx == 1)]
    D = D[(D < np.inf)*(idx == 1)]
    print("Checking %i of %i edges"%(len(I), N*(N-1)/2))
    
    #Prune sparse list and update warped edge lengths (Algorithm 3 pg. 14)
    minlam = np.minimum(lambdaso[I], lambdaso[J])
    maxlam = np.maximum(lambdaso[I], lambdaso[J])
    #Rule out edges between vertices whose balls stop growing before they touch
    #or where one of them would have been deleted.  M stores which of these
    #happens first
    M = np.minimum((E0 + E1)*minlam, E0*(minlam + maxlam))
    #M = E0*(minlam+maxlam)
    
    t = np.arange(len(I))
    t = t[D <= M]
    (I, J, D) = (I[t], J[t], D[t])
    minlam = minlam[t]
    maxlam = maxlam[t]
    
    #Now figure out the metric of the edges that are actually added
    t = np.ones(len(I))
    t[D <= 2*minlam*E0] = 0 #If cones haven't turned into cylinders, metric is unchanged
    #Otherwise, if they meet before the M condition above, the metric is warped
    D[t == 1] = 2.0*(D[t == 1] - minlam[t == 1]*E0) #Multiply by 2 convention
    return (I, J, D)

def makeComplex(X, eps):
    N = X.shape[0]
    #Step 1: Compute all pairwise distances
    XSqr = np.sum(X**2, 1)
    D = XSqr[:, None] + XSqr[None, :] - 2*X.dot(X.T)
    D[D < 0] = 0 #Numerical precision
    D = np.sqrt(D)
    
    I = []
    J = []
    #Step 1.5: If epsilon is zero, we're done (use all edges)
    if eps == 0:
        [I, J] = np.meshgrid(np.arange(N), np.arange(N))
        idx = (I > J)
        (I, J, D) = (I[idx == 1], J[idx == 1], D[idx == 1])
    else:
        #Step 2: Compute greedy permutation
        (perm, lambdas) = getGreedyPerm(D)
        
        #Step 3: Compute the warped edge list
        #Put the insertion radii back in the order of the array
        lambdaso = np.zeros(len(lambdas))
        lambdaso[perm] = lambdas
        (I, J, D) = getSparseEdgeList(lambdaso, eps, D)
    
    print("%i / %i possible edges added (%.3g%%)"%(len(I), N*(N-1)/2, 100*float(len(I))/(N*(N-1)/2)))
    return (I, J, D)

def writeResults(I, J, D, N, filename):
    #Step 5: Write to file
    #First sort edges by weight
    idx = np.argsort(D)
    [I, J, D] = [I[idx], J[idx], D[idx]]
    fout = open(filename, 'w')
    fout.write("p edge %i %i\n"%(N, len(I)))
    for i in range(len(I)):
        fout.write("e %i %i %g\n"%(min(I[i], J[i]), max(I[i], J[i]), D[i]))
    fout.close()

if __name__ == '__main__':
    #Step 1: Load in point cloud
    if len(argv) < 4:
        print("Usage: WarpedRipsMetric <point cloud filename> <epsilon> <outname> <DEBUG (optional)>")
        exit(0)
    #Assumes a point clouds is in a text file with each point on its own line
    #and dimensions separated by spaces    
    #X = sio.loadmat(argv[1])['X']
    #N = X.shape[0]
    fin = open(argv[1])
    X = fin.readlines()
    for i in range(len(X)):
        X[i] = [float(x) for x in X[i].split()]
    X = np.array(X)
    eps = float(argv[2])
    
    (I, J, D) = makeComplex(X, eps)
    writeResults(I, J, D, X.shape[0], argv[3])

