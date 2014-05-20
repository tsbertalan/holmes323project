'''
Created on Apr 2, 2014, 2:27:58 PM

@author: bertalan@princeton.edu
'''

import numpy as np
import kevrekidis as kv
from kevrekidis.utils import setLogLevel, logging


# from tomSims.profiler import profile_method
# @profile_method('out.prof', show_graph=True)
def basisFit(X, Y, p=3, basisNames=["hermite", "legendre"], domains=None, basis=None):
    """
    X : ndarray, shape=(npoints, nvars)
        the abscissae
    Y : ndarray, size=npoints
        the ordinates
    basisNames : __len__=nvars
        Including repetitions of 'hermite' and 'legendre', not case-sensitive.
    domains : list of 2-tuples
        Bounds on the domain for each variable. Optional.
        Passed to np.polynomial.hermite_e.HermiteE etc.
    
    >>> N = 128
    >>> np.random.seed(4)
    >>> X1 = np.random.normal(size=(N,))
    >>> X2 = np.random.uniform(size=(N,))
    >>> X = np.vstack((X1, X2)).T
    >>> Y = X1 ** 2 + 3 * X1 + X2 ** 3 - 5 * X2
    >>> X.shape
    (128, 2)
    >>> co, ca = basisFit(X, Y, p=3)
    >>> abs(np.linalg.norm(Y - ca(X1, X2)) / np.linalg.norm(Y)) < .001
    True
    """
    fitter = BasisFitter(X, Y, p=p, basisNames=basisNames, domains=domains, basis=basis)
    return fitter.fit(X, Y)

    
class BasisFitter(object):
    def __init__(self, X, Y, p=3, basisNames=["hermite", "legendre"], domains=None, basis=None):
        self.npoints = npoints = Y.size
        Y.reshape((npoints, 1))
        nvars = len(basisNames)
        if X.size == Y.size:
            X = X.reshape((npoints, 1))
        assert X.shape == (npoints, nvars), str(X.shape) + str((npoints, nvars))
        if domains is None:
            domains = []
            for i in range(nvars):
                x = X[:,i]
                domains.append((x.min(), x.max()))
        self.basis = productBasis(p, basisNames, domains)
        
    
    def fit(self, X, Y):
    #     from ndHermite import vandermonde
    #     V = vandermonde(basis, X)  # applies all the basis functions to all the points
        V = np.empty((self.npoints, len(self.basis)))
        if X.size == Y.size:
            X = X.reshape((self.npoints, 1))
        for i in range(self.npoints):
            x = tuple(X[i,:])
            for j in range(len(self.basis)):
                f = self.basis[j]
                V[i, j] = f(*x)
        
        norms = np.sqrt(np.square(V.T).sum(1))
        
        def zeroMe(coefs):
            Z = np.zeros((Y.size,))
            for i in range(V.shape[0]):
                for j in range(V.shape[1]):
                    Z[i] += coefs[j] * V[i,j]
                    
    #         Z = np.dot(V / norms, coefs.reshape((len(basis), 1))).ravel()
            res = Z - Y
            if np.isnan(res).any():
                print "ruhroh"
            return res
        
        coInit = np.zeros((len(self.basis),))
        
        co, resid, rank, sing = np.linalg.lstsq(V / norms, Y, rcond=10)
        co = (co.T/norms).ravel()  # Normalizing is important to keep the inversion here^
                                    # well-conditioned if, for example, p and X.max()/X.min()
                                    # are both "large"
        class Ca:
            def __init__(self, coefs, basis):
                self.co = coefs
                self.basis = basis
                assert len(self.co) == len(self.basis)
            def __call__(self, *variables):
                out = 0
                for co, f in zip(self.co, self.basis):
                    out += co * f(*variables)
                return out
                     
        
        ca = Ca(co, self.basis) 
        return co, ca


def hermite1Dbasis(order, domain=[-1, 1]):
    """
    >>> basis = hermite1Dbasis(3)
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> (1.0 * x**3 - 3.0*x).simplify()
    x*(1.0*x**2 - 3.0)
    >>> basis[3](x).simplify()
    x*(1.0*x**2 - 3.0)
    """
    if domain is None:
        domain = [-1, 1]
    basis = []
    for p in range(order+1):
        co = np.zeros((p+1,))
        co[-1] = 1
        
        poly = np.polynomial.hermite_e.HermiteE(co, domain=domain)
        basis.append(poly)
        
    return basis


def legendre1Dbasis(order, domain=[-1, 1]):
    """
    >>> basis = legendre1Dbasis(3)
    >>> from sympy import Symbol
    >>> x = Symbol('x')
    >>> basis[3](x).simplify()
    x*(2.5*x**2 - 1.5)
    
    >>> basis[1](x)
    1.0*x
    
    >>> basis[0](x) == 1
    True
    """
    if domain is None:
        domain = [-1, 1]
    basis = []
    for p in range(order+1):
        co = np.zeros((p+1,))
        co[-1] = 1
        
        poly = np.polynomial.legendre.Legendre(co, domain=domain)
        basis.append(poly)
        
    return basis
    

def productBasis(totalOrder, polySorts, domains=None):
    """
    >>> basis = productBasis(3, polySorts=("hermite", "legendre"))
    >>> len(basis)
    10
    """
    basis = []
    bases1D = []
    nPoly = len(polySorts)
    if domains is None:
        domains = [None] * nPoly
    for i in range(nPoly):
        p = polySorts[i]
        p = p.lower()
        if p == 'hermite':
            bases1D.append(hermite1Dbasis(totalOrder, domain=domains[i]))
        elif p == 'legendre':
            bases1D.append(legendre1Dbasis(totalOrder, domain=domains[i]))
        else:
            raise NotImplementedError, "polynomials of type '%s' not yet implemented." % str(p)
    from itertools import product
    for orders in product(*[range(totalOrder+1) for poly in range(nPoly)]):
        if sum(orders) <= totalOrder:
            
                
            class Ca:
                __doc__ = None
                def __init__(self, bases1D, orders):
                    self.bases1D = bases1D
                    self.orders = tuple(orders)
                    self.nPoly = len(orders)
                    self.nVars = len(bases1D)
                    assert self.nPoly == self.nVars
                    self.__doc__ = "Arity-%d polynomial of orders %s in those %d variables." % (self.nPoly, str(self.orders), self.nPoly)
                def __str__(self):
                    return self.__doc__
                def __repr__(self):
                    return self.__doc__
                def __call__(self, *X):
                    return self.g(*X)
                
                def g(self, *X):
                    out = 1.0
                    orders = self.orders
                    for var in range(self.nVars):
                        basisSize = len(self.bases1D[var])
                        basisFunc = self.bases1D[var][orders[var]]
                        out *= basisFunc(X[var])
                    return out
            ca = Ca(bases1D, orders)
            basis.append(ca)
    return basis
    

def plotTest():
    import doctest
    doctest.testmod()
    np.random.seed(8)
    N = 128
    X1 = np.random.normal(size=(N,1), loc=0, scale=1)
    X2 = np.random.uniform(size=(N,1))
    X = np.hstack((X1, X2))
    Y = X1 ** 2 + 3 * X1 + X2 ** 3 - 5 * X2
    p = 3
    co, ca = basisFit(X, Y, p=p, basisNames=("hermite", "legendre"))
    from ndHermite import show2DFit
    show2DFit(X1, X2, Y, ca, p, labels=("normal", "uniform"), info="hermite and legendre", showSurface=True)
    
    f, a = kv.fa(proj3d=True)
    a.scatter(X1, X2, Y)
    a.scatter(X1, X2, ca(X1, X2))

    
def testSymbolic():
    from sympy.core.symbol import Symbol
    ab = productBasis(3, ("Hermite",))
    f = ab[-1]
    x = Symbol('x')

    for g in f.bases1D[0]:
        print g(x)
    

    for f in ab:
        print f, ":", f(x)


def PDFfromData(x, nbins=10, p=4):
    from scipy.stats.stats import histogram
    y, lowrange, binsize, extrapoints = histogram(x, numbins=nbins)
    x = np.linspace(lowrange, x.max(), nbins)
    co, f = basisFit(x, y, basisNames=('legendre',), p=p, domains=((x.min(), x.max()),))
    normalizer = np.trapz(y, x, dx=x[1]-x[0])
    def ca(x):
        out = f(x) / normalizer
        out = np.array(out)
        if len(out.shape) == 0:
            out = out.ravel()
        out[out<0] = 0
        return out
    return ca


def PDFtest(p=8, mu=3, scale=.1, n=10000, axis=None):
    data = np.random.normal(size=(n,), loc=mu, scale=scale)
    ca1 = PDFfromData(data, nbins=32, p=p)
    ca2 = KernelDensity(data)
    if axis is None:
        f, axis = kv.fa(figsize=(16,9))
    a = axis
    x = np.linspace(data.min(), data.max(), 100)
    a.hist(data, normed=True, bins=32, color=[1,.5,.5,.1], label='data (%d points)' % n)
    a.plot(x, ca1(x), label='O(%d) Legendre fit to histogram' % p)
    a.plot(x, ca2(x), label='kernel density estimation')
    normal = np.exp(-(x-mu)**2 / 2 / scale ** 2)
    normal /= np.trapz(normal, x)
    a.plot(x, normal, label='actual')
    a.legend()
    
    
def KernelDensity(data, h=None):
    n = data.size
    if h is None:
        sigma = np.std(data)
        h = (4. / 3. / n * sigma**5) ** (1./5)
    normalizer = 1.0 / np.sqrt(np.pi * 2)
    def k(u):
        return np.exp(-u**2/2) * normalizer
    def ca(x):
        return sum([k((x - xi) / h) for xi in data]) / n / h
    return ca


def lolScience():
    f, A = kv.fa(numAxes=3, figsize=(20,8))
    for i in range(3):
        n = (10000, 100, 10)[i]
        a = A[i]
        a.set_title(('big data', 'little data', 'shit data')[i])
        PDFtest(n=n, axis=a)


def samplePDF(f, n, domain):
    '''
    Unfortunately, for now, I can't think of how to do this without imposing
    fininte bounds.
    You can use large finite bounds, but making them too large will introduce error.
    
    This code does not check that your PDF integrates to 1.
    '''
    # generate cdf
    x = np.linspace(domain[0], domain[1], n * 2)
    p = f(x)
    from scipy.integrate import cumtrapz
    c = cumtrapz(p, x=x)
    
    
    interpolator = LinInterp(c, x[1:])
#     xf = np.linspace(c.min(), c.max(), n*10)
#     a.scatter(interpolator(xf), xf)
    
    u_u = np.random.uniform(low=0, high=1, size=(n,))
    o_o = interpolator(u_u)
    
    return o_o
    
    
def LinInterp(xx, yy):
    '''Returns a linear interpolant function.'''
    xx = np.array(xx)
    yy = np.array(yy)
    assert xx.size == yy.size
    xmin = xx.min()
    xmax = xx.max()
    assert xmin == xx[0] and xmax == xx[-1], 'x should be sorted' 
    
    @np.vectorize
    def interp(x):
        i = kv.utils.findClosest(xx, x)
        xi = xx[i]
        
        if i == 198:
            pass
        
#         # extrapolation
#         if x <= xmin:
#             j = i + 1
#         elif x >= xmax:
#             j = i - 1

        assert x >= xmin or x <= xmax, str((x, xmin, xmax))

        # interpolation
        if x < xi:
            # direction = left
            j = i - 1
        else:
            assert x >= xi
            # direction = right
            j = i + 1
        
             
        yi = yy[i]
        m = (yi - yy[j]) / (xi - xx[j])
        y = m * (x - xi) + yi
        return y
    
    return interp
    

def testSample(n=1000):
    for pdf in (
                lambda x: np.exp(-x**2/2) / np.sqrt(2 * np.pi),
                lambda x: (np.exp(-(x-2)**2/.25)/2 + np.exp(-(x+1)**2/.2))/ (
                                    (5 + 4*np.sqrt(5)) * np.sqrt(np.pi) / 20
                                                                          )
                ):
        
        sample = samplePDF(pdf, n, (-4, 4))
        
        f, a = kv.fa()
        a.hist(sample, normed=True, color=(.5, 1, .5, .5), bins=max(10, n/10), cumulative=True)
        a.hist(sample, normed=True, color=(1, .5, .5, .5), bins=max(10, n/10))
        a.scatter(sample, [-.1]*n)


def basisDiscreteData(xx, order):
    '''
    >>> n = 100
    >>> A = kv.graphs.generators.chungLu(n, p=1.0, r=0.4, symmetric=True)
    >>> data = kv.graphs.degreeVec(A)
    >>> p = basisDiscreteData(data, 3)
    >>> len(p)
    4
    >>> assert abs(p[0](2) - 1) < .1, str(abs(p[0](2) - 1))
    >>> assert abs(p[1](2) - 2) < .1, str(abs(p[1](2) - 2))
    >>> assert abs(p[2](2) - 3) < .1, str(abs(p[2](2) - 3))
    '''
    
    N = xx.size
    p = range(order+1)
    p[0] = lambda x: 1.0
    
    h_, bin_edges = np.histogram(xx, bins=N)
    h = lambda i: h_[i-1]
    
    def inner(a, b):
        return sum([a(i) * b(i) * h(i) for i in range(1, N+1)])

    inners = range(order+1) # cache some work
    inners[0] = inner(p[0], p[0])
    
    
    def p1(i):
        ins = inner(lambda u: u * 1.0, lambda u: 1.0)
        sub = ins / inners[0]
        p0i = p[0](i)
        out = (i - sub) * p0i
        return out
    p[1] = p1
    
    inners[1] = inner(p[1], p[1])
    
    from sympy import Symbol
    x_ = Symbol('x_')
    
    for j in range(1, order):
        inners[j] = inner(p[j], p[j])
        
        sub1 = inner(lambda u: j * p[j](u), p[j]) / inners[j]
        sub2 = inners[j] / inners[j-1]
        p[j+1] = lambda i: (i - sub1) * p[j](i) - sub2 * p[j-1](i)
       
        coeffs = p[j+1](x_).as_poly().all_coeffs()
        assert len(coeffs) == j+2, str((len(coeffs), j+2, coeffs, p[j+1](x_)))
        logging.debug(str(
                        (p[j+1](x_), coeffs, j)
                        ))
        p[j+1] = lambda x: sum([a * x ** (j-i) for i, a in enumerate(coeffs)]) 
    
    def Antivectorize(g):
        def m(x):
            if not hasattr(x, '__getitem__'):
                return g(x).ravel()[0]
            else:
                return g(x)
        return m
    
    return [Antivectorize(g) for g in [np.vectorize(f) for f in p]]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    kv.exit()
    
    setLogLevel('debug')
#     PDFtest()
#     testSample()
    
    n = 100
    A = kv.graphs.generators.chungLu(n, 1, .4, symmetric=True)
    data = kv.graphs.degreeVec(A)
    p = basisDiscreteData(data, 6)
    f, a = kv.fa()
    x = np.linspace(-5, 5, n)
    for i, f in enumerate(p):
        y = f(x)
        a.plot(x, y, label=r'$p_%d(x)$' % i)
    a.legend()
    
#     plotTest()
    
#     testSymbolic()
    kv.plotting.show()
