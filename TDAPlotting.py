"""
Plotting Tools for persistence diagrams and cocycles
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def plotDGM(dgm, color = 'b', sz = 20, label = 'dgm', axcolor = np.array([0.0, 0.0, 0.0]), marker = None):
    if dgm.size == 0:
        return
    # Create Lists
    # set axis values
    axMin = np.min(dgm)
    axMax = np.max(dgm)
    axRange = axMax-axMin
    a = max(axMin - axRange/5, 0)
    b = axMax+axRange/5
    # plot line
    plt.plot([a, b], [a, b], c = axcolor, label = 'none')
    #plt.hold(True)
    # plot points
    if marker:
        H = plt.scatter(dgm[:, 0], dgm[:, 1], sz, color, marker, label=label, edgecolor = 'none')
    else:
        H = plt.scatter(dgm[:, 0], dgm[:, 1], sz, color, label=label, edgecolor = 'none')
    # add labels
    plt.xlabel('Time of Birth')
    plt.ylabel('Time of Death')
    return H

def plotDGMAx(ax, dgm, color = 'b', sz = 20, label = 'dgm'):
    if dgm.size == 0:
        return
    axMin = np.min(dgm)
    axMax = np.max(dgm)
    axRange = axMax-axMin;
    ax.scatter(dgm[:, 0], dgm[:, 1], sz, color,label=label)
    #ax.hold(True)
    ax.plot([axMin-axRange/5,axMax+axRange/5], [axMin-axRange/5, axMax+axRange/5],'k');
    ax.set_xlabel('Time of Birth')
    ax.set_ylabel('Time of Death')

def plot2DGMs(P1, P2, l1 = 'Diagram 1', l2 = 'Diagram 2'):
    plotDGM(P1, 'r', 10, label = l1)
    #plt.hold(True)
    plt.plot(P2[:, 0], P2[:, 1], 'bx', label = l2)
    plt.legend()
    plt.xlabel("Birth Time")
    plt.ylabel("Death Time")


def plotTriangles(X, A, B, C):
    #plt.hold(True)
    ax = plt.gca()
    for i in range(len(A)):
        poly = [X[A[i], :], X[B[i], :], X[C[i], :]]
        ax.add_patch(Polygon(np.array(poly), linestyle='solid', color='#00FF00', alpha=0.05))

def drawLineColored(X, C):
    #plt.hold(True)
    for i in range(X.shape[0]-1):
        plt.plot(X[i:i+2, 0], X[i:i+2, 1], c=C[i, :], lineWidth = 3)

def plotCocycle(X, cocycle, thresh, drawTriangles = False):
    XSqr = np.sum(X**2, 1)
    D = XSqr[:, None] + XSqr[None, :] - 2*X.dot(X.T)
    D[D < 0] = 0 #Numerical precision
    D = np.sqrt(D)

    #plt.hold(True)
    ax = plt.gca()
    #Plot all edges under the threshold
    N = X.shape[0]
    t = np.linspace(0, 1, 10)
    c = plt.get_cmap('Greys')
    C = c(np.array(np.round(np.linspace(0, 255, len(t))), dtype=np.int32))
    C = C[:, 0:3]

    for i in range(N):
        for j in range(N):
            if D[i, j] <= thresh:
                Y = np.zeros((len(t), 2))
                Y[:, 0] = X[i, 0] + t*(X[j, 0] - X[i, 0])
                Y[:, 1] = X[i, 1] + t*(X[j, 1] - X[i, 1])
                drawLineColored(Y, C)
                #ax.arrow(X[i, 0], X[i, 1], X[j, 0] - X[i, 0], X[j, 1] - X[i, 1], fc="k", ec="k", head_width=0.2, head_length=0.5)
    #Plot cocycle
    for k in range(cocycle.shape[0]):
        [i, j, val] = cocycle[k, :]
        [i, j] = [min(i, j), max(i, j)]
        #plt.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], 'r', lineWidth = 3, linestyle='--')
        #ax.arrow(X[i, 0], X[i, 1], X[j, 0] - X[i, 0], X[j, 1] - X[i, 1], fc="k", ec="k", head_width=0.05, head_length=0.1)
        a = 0.5*(X[i, :] + X[j, :])
        plt.text(a[0], a[1], '%g'%val)

    #Enumerate Triangles
    if drawTriangles:
        N = X.shape[0]
        [A, B, C] = np.meshgrid(np.arange(N), np.arange(N), np.arange(N))
        [A, B, C] = [A.flatten(), B.flatten(), C.flatten()]
        tidx = np.arange(len(A), dtype=np.int32)
        tidx = tidx[(D[A, B] <= thresh)*(D[B, C] <= thresh)*(D[A, C] <= thresh)]
        [A, B, C] = [A[tidx], B[tidx], C[tidx]]
        plotTriangles(X, A, B, C)

    #Plot X
    plt.scatter(X[:, 0], X[:, 1], 100, 'k')
    for i in range(X.shape[0]):
        plt.text(X[i, 0]+0.02, X[i, 1]+0.02, '%i'%i)
