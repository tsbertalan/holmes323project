#!/usr/bin/env python
'''
Created on May 8, 2014

@author: tsbertalan
'''
import logging
from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import axes3d
import kevrekidis as kv

from leloup import LeloupParameters, LofT, leloupX0params
from leloup import diff as leloupDiff
import basisFit


class ToParameters(LeloupParameters):
    def __init__(self, forcingShape=None, *args, **kwargs):
        super(ToParameters, self).__init__(intParams=True, *args, **kwargs)
        self.forcingShape = forcingShape
    def loadDefaults(self):
        from parameters.toDefaults import paramsString
        super(ToParameters, self).loadDefaults()  # load the Leloup defaults parameters
        self.addFromList(paramsString.split('\n'))

        
def NetDiff(ALPHA, saves, parameters=None, hetNames=['vsp0'], hetVecs=[None], doSave=False, startTime=0):
    from time import time
    import copy
    if parameters is None:
        parameters = ToParameters()
        parameters.loadDefaults()
    parameters.rowsums = ALPHA.sum(axis=1)
    
    L = parameters.L
    Lbase = parameters.Lbase
        
    N = parameters.N
    plist = [copy.copy(parameters) for i in range(N)]
    for hetName, hetVec in zip(hetNames, hetVecs):
        if hetVec is None:
            hetVec = np.random.normal(loc=parameters[hetName], scale=parameters[hetName] * .1, size=(N,))
        for i in range(N):
            plist[i].addFromDict({hetName: hetVec[i]})
    if doSave:
        saves['rhoMax'] = []
    tlast = [0,time()]
    def netDiff(X, t):
        if tlast[0] != t:
            tirl = time()
            logging.debug("t=%f hr, simHrs/irlSecs=%.5f" % (t, (t-tlast[0])/(tirl - tlast[1])))
            tlast[:] = t, tirl
        assert ALPHA.shape[0] == ALPHA.shape[1] == N, "c.v. %s and %d" % (str(ALPHA.shape), N)
        out = np.empty((N * 17+1,))
        assert L > Lbase
        Lt = LofT(t, L, Lbase, shape=parameters.forcingShape)
        
        if (L - Lt) < 0.1 * L:
            RHO = [plist[i].a for i in range(N)]
        else:
            RHO = np.empty((N,))
            MP = X.reshape(N, 17)[:,0].ravel()
            for i in range(N):
                # eqn 1
                RHO[i] = parameters.a * MP[i] / (MP[i] + parameters.b)
        if doSave:
            saves['rhoMax'].append((t, max(RHO)))
        for i in range(N):
            a = i*17
            b = (i+1)*17
            out[a:b] = cellDiff(X[a:b], t, ALPHA, RHO, i, parameters=plist[i], Lt=Lt, L=L, Lbase=Lbase, forcingShape=parameters.forcingShape)
            
        return out
    return netDiff


def cellDiff(X, t, ALPHA, RHO, icell, parameters=None, Lt=0.0, L=.2, Lbase=0.0, forcingShape=None):
    
    if parameters is None:
        parameters = ToParameters()
        parameters.loadDefaults()
    p = parameters
    out = np.empty((17,))
    
    CBstar = X[16]
    
    epsilon = p.epsilon
    
    if (L - Lt) < 0.1 * L:  # we're in the dark phase
#         rho = p.a
#         assert 'rowsums' in p.__dict__
#         gamma = p.rowsums[icell] * rho / epsilon
        beta = 1.0  # This should amount to the same thing as the above three lines,
                    #  and is here to match the statement from To et al. that
                    #  "during light, the VIP release rate was modeled as constant
                    #   and sufficiently high to cause complete saturation fo VPAC2
                    #   receptors." 
    else:
#         gamma = sum([  # equivalent to the below, but slower.
#                          ALPHA[icell,j]*RHO[j]
#                          for j in range(p.N)
#                          ]) / epsilon
        gamma = (ALPHA[icell,:] * RHO).sum() / epsilon  # [2a] 
        beta = gamma / (p.KD + gamma)  # [3d]
    Ca2pCytosol = (p.v0 + p.v1 * beta + p.v2 * Lt) / p.k  # eqn 4
    vK = p.VMK * Ca2pCytosol / (p.Ka + Ca2pCytosol)  # eqn 5b
    
    dCBstardt = p.vP / p.CBT * (  # 5a
                            (vK/p.vP)*(1-CBstar) / (p.K1 + 1 - CBstar)
                            -CBstar / (p.K2 + CBstar)
                            )
    out[16] = dCBstardt

    
    l = p.CBT * CBstar / (p.KC + p.CBT * CBstar)
    vsP = p.vsP0 + p.CT * l

    out[:16] = leloupDiff(X[:16], t, parameters=parameters, L=L, Lbase=Lbase, vsP=vsP, forcingShape=forcingShape)
    
    return out


def toX0params():
    out = []
    out.extend(leloupX0params())
    p = LeloupParameters()
    p.loadDefaults()
    initial = {
               'CBstar': 0.12  # from Vasalou 2010 PLoS 
              } 
    p.addFromDict(initial)
    out.append(p.CBstar)
    return out


def toX0NetParams(N):
    X0 = np.hstack([toX0params()] * N).tolist()
    return np.array(X0)    


def netIntegrate(ALPHA, parameters, initial=None, titleInfo="", nstep=1e3, **kwargs):
    startTime = time()
    if initial is None:
        initial = toX0NetParams(parameters.N)# * np.random.uniform(low=.6, high=1.4, size=(parameters.N*17+1,))
    if 'epsilon' not in parameters.pnames:
        parameters.addFromDict({"epsilon": 1./parameters.N * sum([sum([ALPHA[i,j] for j in xrange(parameters.N)]) for i in xrange(parameters.N)])})   # [2b]
    hetName = 'vsP0'
    hetLatex = r'$v_{sP0}$'
    hetUnits = r'[nM h$^{-1}$]'
    loc = parameters[hetName] + .4
    scale = parameters[hetName] * .1
    hetVec = np.random.uniform(low=loc-scale, high=loc+scale, size=(parameters.N,))
    saves = None
    _diff = NetDiff(ALPHA, saves, parameters=parameters, hetNames=[hetName], hetVecs=[hetVec], startTime=startTime)
    X, T = kv.dynSys.integration.integrate(initial, _diff, tmin=0, giveTime=True, nstep=nstep, **kwargs)
    assert parameters.L > parameters.Lbase
    
    saves = {}
    saves['hetName'] = hetName
    saves['hetLatex'] = hetLatex
    saves['hetUnits'] = hetUnits
    saves['hetVec'] = hetVec
    saves['initial'] = initial
    saves['ALPHA'] = ALPHA
    saves['parameters'] = parameters
    saves['LofT'] = LofT(T, parameters.L, parameters.Lbase, shape=parameters.forcingShape)
    saves['fname'] = "netIntegrate%d" % int(time()) + titleInfo
    
    from os import system
    system('mkdir -p data')
    
    npzName = 'data/' + saves['fname'] + '.npz'
    logging.debug('saving %s' % npzName)
    np.savez(npzName, X=X, T=T, **saves)
    
    return saves['fname']


def mRNAPlot(X, T, N, Lt, fname, show=False):
    fig, ax = kv.fa(figsize=(10,6))
    X0p = toX0params()
    ncurves = 3
    colors = kv.plotting.colors(numColors=ncurves)
    for i in range(ncurves):
        for j in range(N):
            x = X0p[i]
            if j == 0:
                label = x.glabel()
            else:
                label = None
            a = j * 17 + i
            ax.plot(T, X[:,a], label=label, color=colors[i])

    ax2color = 'red'
    ax2 = ax.twinx()
    ax2.xaxis.label.set_color(ax2color)
    ax2.plot(T, Lt, label=r"$L$", color=ax2color, zorder=-4)
    
    ax.set_xlabel('time [hr]')
    ax.set_ylabel('mRNA concentration [nM]')

    ax2.set_ylabel('light level $L$', color=ax2color)
    [i.set_color(ax2color) for i in ax2.get_yticklabels()]

    leg1 = ax.legend(loc='lower right')
    ax2.legend(loc="upper right")
    ax2.add_artist(leg1) 
    ax.legend = None 
    
    fig.subplots_adjust(left=.05, bottom=.08, right=.9, top=.95)
    
    if show:
        kv.plotting.show()
    kv.plotting.saveFig(fig, fname)
     


def plotSaved(fname, makeMovie=True, titleInfo=None, show=False):
    
    if titleInfo is None:
        titleInfo = fname
    
    npzName = str(fname)
    if '.npz' not in npzName:
        npzName += '.npz'
    if 'data/' not in npzName:
        npzName = 'data/' + npzName
    d = kv.utils.npzLoadData(npzName)
    
    sortable = np.array(d.hetVec)
    sortable.sort()
    
    def makeFrame(data, i, fig, axes, **kwargs):
        X, T = data
        ax = axes[0]
        ax.cla()
        
        j = np.linspace(0, len(T)-1, NFrames).astype(int)[i]
        MP = X.reshape(X.shape[0], d.parameters.N, 17)[j,:,0]
        MPall = X.reshape(X.shape[0], d.parameters.N, 17)[:,:,0]
        
        order = 2
        degrees = kv.graphs.degreeVec(d.ALPHA).T

        # The actual PCE fit:
        #fitter2D = kv.dynSys.eqnFree.polyChaos.basisFit.BasisFitter(
        fitter2D = basisFit.BasisFitter(
                                                    np.vstack((d.hetVec, degrees)).T
                                                                    , MP, p=order, basisNames=['legendre', 'hermite'])
        co2D, ca2D = fitter2D.fit(np.vstack((d.hetVec, degrees)).T, MP)
        kv.dynSys.eqnFree.polyChaos.ndHermite.show2DFit(d.hetVec, degrees, MP, ca2D, p=order, show=False,
                                                        figax=(fig, ax), showScatter=True, showSurface=True,
                                                        title="", 
                                                        labels=(d.hetLatex, 'weighted degree', '$M_P$'))
        ax.set_zlim(MPall.min(), MPall.max())

        
        fig.subplots_adjust(left=0, bottom=.16, top=.9, right=.95)
        ax.set_title(r'$t=%.1f$ [hr]' % T[j])
        
    from os import system
    system('mkdir -p images')
    animFig = kv.plotting.plt.figure(figsize=(8,4.5))
    
    animAxes = [
                animFig.add_subplot(1,1,1, projection='3d')
                ]
    
    if makeMovie:
        NFrames = max(64, int(64 * (d.T.max() - d.T.min()) / 24.0))  # 64 frames/day; at least 64 frames
        logging.debug('movie with %d frames' % NFrames) 
        (animFig, animAxes), animName = kv.plotting.animateOnAxis((d.X, d.T), makeFrame,
                                          fileBaseName='images/' + fname, NFrames=NFrames,
                                          suptitle="", dpi=720/4.5, figAxes=(animFig, animAxes),
                                          verbose=True, numAxes=2, show=True, closePlot=False,
                                          )
        kv.plotting.saveFig(animFig, animName+'finalFit', allowPeriod=False)
    mRNAPlot(d.X, d.T, d.parameters.N, d.LofT, fname, show=show)
    
    return d.X, d.T


def gaussMesh(N, xsize, ysize, eps=32.0, probScale=1.0, show=False):
    Xlocs = np.arange(xsize)
    Ylocs = np.arange(ysize)
    from itertools import product
    locs = [loc for loc in product(Xlocs, Ylocs)]
    np.random.shuffle(locs)
    Xlocs = [locs[i][0] for i in range(N)]
    Ylocs = [locs[i][1] for i in range(N)]
    X, Y = np.meshgrid(Xlocs, Ylocs)
    D = np.sqrt((X - X.T) ** 2 + (Y - Y.T) ** 2)
    PROB = probScale * np.exp(-D ** 2 / 2. / eps)
    PROB[PROB > 1] = 1
    PROB[np.eye(N) == 1] = 0  # no self-connections
    ALPHA = (np.random.uniform(low=0, high=1, size=(N, N)) < PROB).astype(int)
    
    if show:
        kv.plotting.showMat(ALPHA, title=r"$\alpha$")
        fig, ax = kv.fa()
        ax.hist(kv.graphs.degreeVec(ALPHA), bins=10)
     
        sortX, sortY = np.meshgrid(np.argsort(Y[0,:]), np.argsort(Y[0,:]))
        sorter0 = sortX.T, sortY.T
         
        kv.plotting.showMat(PROB[sorter0], title=r"$p$")
        kv.plotting.show()
        kv.plotting.ubiGraphPlot(ALPHA[sorter0])
        kv.exit()
    return ALPHA


def gridMesh(N, show=False):
    rows = int(np.floor(float(N)**.5))
    done = False
    while not done:
        cols = float(N) / rows
        if cols == int(cols):
            done = True
        else:
            rows += 1
    assert rows * cols == N
    if rows / cols > 3:
        logging.warn("Long/thin geometry (%d by %d) for N=%d." % (rows, cols, N)
                     + "Try an N with better factors, such as a perfect square.")  
    X, Y = np.arange(cols), np.arange(rows)
    ALPHA = np.empty((N, N))
    from itertools import product
    for i, j in product(range(N), range(N)):
        if i==j:
            ALPHA[i,j] = 1.0
        else:
            rowa = i//cols
            rowb = j//cols
            cola = i - rowa * cols
            colb = j - rowb * cols
            ALPHA[i,j] = ((X[cola] - X[colb])**2 + (Y[rowa] - Y[rowb])**2)**-.5
    if show:
        kv.plotting.showMat(ALPHA, title=r"$\alpha_{i,j}$ (geometry is $%d \times %d = %d$)" % (rows, cols, N), show=False)
        kv.plotting.ubiGraphPlot(ALPHA > .8) #try .7 for a triangular grid, or .8 for a square grid. Smaller values might develop thickness.
    return ALPHA
    

def main(N=16, xsize=128, ysize=128, tmax=72, **kwargs):
    from os import system
    task = '(%s) N=%d, tmax=%.1f' % (__name__, N, tmax)
    start = time()
    kv.utils.notify('began %s' % task)
    p = ToParameters(forcingShape=None)
    p.loadDefaults()
    p.addFromDict({'N': N})
#     ALPHA = gaussMesh(N, xsize, ysize)  # an alternate approach--let connection strength die as a Gaussian kernel of distance.
    ALPHA = gridMesh(N, show=False)
#     ALPHA = np.ones((N,N))  # indiscriminate all-all connectivity
    out = netIntegrate(ALPHA, p, tmax=tmax, **kwargs)
    stop = time()
    kv.utils.notify('done with %s' % task, '%f min elapsed' % ((stop-start) / 60.))
    return out
    

def saveToSystem():    
    from leloupModelPrinter import saveSystem
    IC = toX0params()
    params = ToParameters()
    params.loadDefaults()
    saveSystem('to', cellDiff, IC, params)


def plotAll():
    for o in kv.utils.doCmd('ls data').split('\n'):
        try:
            d = kv.utils.npzLoadData('data/' + o)
            plotSaved(o, makeMovie=False, titleInfo=o)
        except:
            print 'failed for path "data/%s"' % o


def cmdLine():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", "--number", type=int, required=True)
    parser.add_argument("-t", "--tmax", type=float, default=8.0, required=True)
    parser.add_argument("--makeMovie", type=bool, default=True)
    parser.add_argument("-l", "--label", type=str, default="", help="extra info to put in the saved file's name")
    parser.add_argument("-f", "--file", type=str, default="", help="Display a file rather than running a simulation.")
    args = parser.parse_args()
    if args.file == "":
        result = main(N=args.number, tmax=args.tmax, titleInfo=args.label)
    else:
        result = args.file
    kv.utils.notify(result)
    plotSaved(result, makeMovie=args.makeMovie, titleInfo=args.label)
    

if __name__ == '__main__':
    from tomSims import setLogLevel
    setLogLevel('debug')  # to display progress information during time-integration

#     result = main(N=32, tmax=4, nstep=1e2)
#     plotSaved(result, makeMovie=True, show=True)

    cmdLine()
     
    kv.plotting.show()
