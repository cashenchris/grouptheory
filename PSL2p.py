import numpy as np
from numty.quadResidues import legendre
import itertools

# requires package numty

# modular matrix functions
def modinverse(a,p):
    """
    Multiplicative inverse of a mod p for prime p.
    """
    if a%p==0:
        raise ZeroDivisionError
    else:
        return pow(a,p-2,p)
    
def moddot(M,N,p):
    """
    Product of matrices M and N mod p.
    """
    return M.dot(N)%p

def moddet(M,p):
    """
    Detrminant of 2x2 matrix mod p.
    """
    return (M[0][0]*M[1][1]-M[0][1]*M[1][0])%p

def modmatrixinverse(M,p):
    """
    inverse of 2x2 matrix mod p.
    """
    N=np.array([[M[1][1],-M[0][1]],[-M[1][0],M[0][0]]])
    return modinverse(moddet(M,p),p)*N%p

def modconj(M,N,p):
    """
    Return N^(-1}MN.
    """
    return moddot(moddot(modmatrixinverse(N,p),M,p),N,p)

def modtrace(M,p):
    """
    Return the trace of M mod p.
    """
    return (M[0][0]+M[1][1])%p

def modmatrixpower(M,n,p):
    """
    Raise matrix to a power mod a prime.
    """
    if n==0:
        return np.array([[1,0],[0,1]])
    elif n<0:
        return modmatrixpower(modmatrixinverse(M,p),abs(n),p)
    elif n==1:
        return M%p
    else:
        return modmatrixpower(M,n%2,p).dot(modmatrixpower(M.dot(M)%p,n//2,p))%p

def nonresidue(p):
    """
    Returns the smallest integer 0<y<p such that y is not a residue mod p.
    """
    for y in range(1,p):
        if legendre(y,p)<0:
            return y

def modsquareroot(x,p):
    """
    Return a square root of x mod p, where p is an odd prime, or raise error if x is a non-residue.
    """
    if x%p==0:
        return 0
    elif x%p==1:
        return 1
    else:
        for y in range(2,(p+1)/2):
            if x%p==(y*y)%p:
                return y
        else:
            raise ValueError(str(x)+' is not a residue mod '+str(p))

        
#----------------    
class PSLelement(object):
    """
    An element of PSL(2,p) defined by a matrix M.
    """
    def __init__(self,M, p):       
        if (M[0][0]*M[1][1]-M[0][1]*M[1][0]-1)%p==0:
            thematrix=PSLRep(M,p)
            self.M=((thematrix[0][0],thematrix[0][1]),(thematrix[1][0],thematrix[1][1]))
            self.p=p
        else:
            raise ValueError('Input does not have determinant 1.')

    def __hash__(self):
        return hash((self.M,self.p))

    def __repr__(self):
        return str(self.matrix())+' % '+str(self.p)

    def __str__(self):
        return '[['+str(self.M[0])+'],['+str(self.M[1])+']]%'+str(self.p)

    def __eq__(self,other):
        return self.M==other.M and self.p==other.p

    def matrix(self):
        return np.array([[self.M[0][0],self.M[0][1]],[self.M[1][0],self.M[1][1]]])

    def __getitem__(self,sliced):
        return self.M[sliced]
    
    def prime(self):
        return self.p
    
    def __mul__(self,other):
        assert(self.p==other.p)
        return PSLelement(self.matrix().dot(other.matrix()),self.p)

    def inverse(self):
        return PSLelement(modmatrixinverse(self.matrix(),self.p),self.p)

    def trace(self):
        return (self.M[0][0]+self.M[1][1])%self.p

    def is_trivial(self):
        return self.M==((1,0),(0,1)) or self.M==((self.p-1,0),(0,self.p-1))

    def __pow__(self,thepower):
        return PSLelement(modmatrixpower(self.matrix(),thepower,self.p),self.p)
    
class PSL2(object):
    """
    The group PSL(2,p) for an odd  prime p.
    """
    def __init__(self,p):
        self.theprime=p

    def element(self,M):
        """
        Take a matrix and return the correspoding element of PSL(2,p).
        """
        return PSLelement(M,self.theprime)

    def identityelement(self):
        return PSLelement(np.array([[1,0],[0,1]]),self.theprime)

    def is_trivial(self,M):
        try:
            return PSListrivial(M.matrix(),self.theprime)
        except AttributeError:
            return PSListrivial(M,self.theprime)
    
    def are_conjugate(self,M,N):
        try:
            Amatrix=M.matrix()
        except AttributeError:
            Amatrix=M
        try:
            Bmatrix=N.matrix()
        except AttributeError:
            Bmatrix=N
        t1=M.trace()**2
        t2=N.trace()**2
        if t1!=t2:
            return False
        elif t1!=4:
            return t1==t2
        elif self.is_trivial(M) or self.is_trivial(N):
            return self.is_trivial(M) and self.is_trivial(N)
        else:
            return self.specialconj(M)==self.specialconj(N)

    def specialconj(self,M):
        """
        Given a non-identity matrix with trace**2=4, check if it is conjugate to [[1,1],[0,1]] in PSL(2,p).
        """
        if M[0][1]%self.theprime==0:
            return legendre(-M[1][0]**(self.theprime-2),self.theprime)==1
        else:
            return legendre(M[0][1],self.theprime)==1
    
    def are_automorphic(self,M,N):
        """
        Determine if there is an automorphism that takes M to N.
        """
        # For non-identity elements there is such an automorphism if and only if M and N have the same trace^2.
        try:
            Amatrix=M.matrix()
        except AttributeError:
            Amatrix=M
        try:
            Bmatrix=N.matrix()
        except AttributeError:
            Bmatrix=N
        t1=M.trace()**2
        t2=N.trace()**2
        if t1!=t2:
            return False
        elif t1!=4:
            return True
        else:
            return self.is_trivial(M)==self.is_trivial(N)

    def conjugator(self,M,N=None):
        """
        If two matrices are given return an element P of PSLsp such that P**(-1)MP=N. If one matrix is given return an element P of PSL2p such that P**(-1)MP is the favorite representative of the conjugacy class of M.
        """
        if N is not None:
            if self.are_conjugate(M,N):
                return self.conjugator(M)*(self.conjugator(N).inverse())
            else:
                raise ValueError('elements are not conjugate')
        else:
            p=self.theprime
            try:
                Mmat=M.matrix()
            except AttributeError:
                Mmat=self.element(M).matrix()
            tr=(Mmat[0][0]+Mmat[1][1])%p
            y=Mmat[0][1]
            z=Mmat[1][0]
            x=(Mmat[0][0]-tr*modinverse(2,p))%p
            if tr==2 and x==0 and y==0 and z==0:
                return self.identityelement()
            # compute matrix [[a,b],[c,d]] with det=1 whose inverse is P
            if (tr**2-4)%p!=0: #rep is [[tr/2, 1],[tr^2/4-1,tr/2]]
                if y%p==0:
                    b=1 #arbitrary nonzero
                    d=-b*x
                    a=-(z*b**2+1)*modinverse(2*x*b,p)
                    c=(z*b**2-1)*modinverse(2*b,p)
                else:
                    for (b,d) in itertools.product(range(p),range(p)):
                        if (d**2-(tr**2*modinverse(4,p)-1)*b**2-y)%p==0:
                            break
                    # b,d arbitrary satisfying y=d^2-b^2(tr^2/4-1)
                    a=(d+b*x)*modinverse(y,p)
                    c=(b*(tr**2*modinverse(4,p)-1)+d*x)*modinverse(y,p)
            else: #rep is [[1,n],[0,1]] where n is either 1 or smallest non-residue mod p, depending on class of M
                if y%p==0:
                    if legendre(-z,p)>0:
                        n=1
                    else:
                        for n in range(p):
                            if legendre(n,p)<0:
                                break
                    # -z/n is a residue
                    for c in range(p):
                        if (c**2+z*modinverse(n,p))%p==0:
                            break
                    # c=sqrt(-z/n)
                    b=-modinverse(c,p)
                    d=0
                    a=0 #arbitrary
                else:
                    if legendre(y,p)>0:
                        n=1
                    else:
                        for n in range(p):
                            if legendre(n,p)<0:
                                break
                    # y/n is a residue
                    for d in range(p):
                        if (d**2-y*modinverse(n,p))%p==0:
                            break
                    # d=sqrt(y/n)
                    b=0 #arbitrary
                    c=x*modinverse(d*n,p)
                    a=modinverse(d,p)+b*x*modinverse(y,p)
            return self.element(np.array([[a,b],[c,d]])).inverse()
                

     

def PSL2generators(p):
    """
    Return a list of generators of PSL2p.
    """
    return [np.array([[1,1],[0,1]]),np.array([[1,0],[1,1]])]

def PGL2generators(p):
    """
    Return a list of generators of PGL2p.
    """
    if p%4==3: # -1 is not a residue
        return PSL2generators(p)+[np.array([[1,0],[0,p-1]])]
    elif p%8==5:# 2 is not a residue
        return PSL2generators(p)+[np.array([[1,0],[0,2]])]
    else:
        for a in range(p):
            if legendre(a, p)<0:
                return PSL2generators(p)+[np.array([[1,0],[0,a]])]
        else:
            raise ArithmeticError

def PSLRep(M,p):
    """
    Choose a representative matrix in the same PSL2p class as M.
    """
    Mat=np.array([[M[0][0],M[0][1]],[M[1][0],M[1][1]]])%p
    t=(Mat[0][0]+Mat[1][1])%p
    if t>0 and t<=(p-1)/2:
        return Mat
    elif t>0:
        return -Mat%p
    elif Mat[0][0]>0 and Mat[0][0]<=(p-1)/2:
        return Mat
    elif Mat[0][0]>0:
        return -Mat%p
    elif Mat[0][1]<=(p-1)/2:
        return Mat
    else:
        return -Mat%p


def PSL2p_generator(p):
    """
    Generator that yields elements of PSL2p
    """
    if p==2:
        for M in [[[1,0],[0,1]],[[0,1],[1,0]],[[1,1],[0,1]],[[1,0],[1,1]],[[1,1],[1,0]],[[0,1],[1,1]]]:
            yield PSLelement(np.array(M),p)
        return
    for b in range(1,1+(p-1)/2):
        yield PSLelement(np.array([[0,b],[-modinverse(b,p),0]]),p)
    for a in range(1,1+(p-1)/2):
        if (a*a)%p==p-1:
            for b in range(p):
                yield PSLelement(np.array([[a,b],[0,(-a)%p]]),p)
        for c in range(1,p):
            yield PSLelement(np.array([[a,(-(1+a*a)*modinverse(c,p))%p],[c,(-a)%p]]),p)
    for d in range(1,1+(p-1)/2):
        for b in range(1,p):
            yield PSLelement(np.array([[0,b],[(-modinverse(b,p))%p,d]]),p)
    for a in range(1,p):
        for d in range(p):
            if (a+d)%p>0 and (a+d)%p<=(p-1)/2:
                if (a*d)%p==1:
                    for b in range(p):
                        yield PSLelement(np.array([[a,b],[0,d]]),p)
                    for c in range(1,p):
                        yield PSLelement(np.array([[a,0],[c,d]]),p)
                else:
                    for b in range(1,p):
                        yield PSLelement(np.array([[a,b],[(-(1-a*d)*modinverse(b,p))%p,d]]),p)
 

def favoritePSLrep(theMatrix,p):
    """
    Given a matrix theMatrix in SL2p return R in PSL2p, and C in PGL2p such that theMatrix=C^-1RC in PSL2p and R is of the form [[1,0],[0,1]] or [[trace/2,1],[trace**2/4-1, trace/2]]
    """
    if modtrace(theMatrix,p)>(p-1)/2:
        M=-theMatrix
    else:
        M=theMatrix
    if PSListrivial(M,p):
        R=np.array([[1,0],[0,1]])
        C=np.array([[1,0],[0,1]])
        return R,C
    t=modtrace(M,p)
    T=(t**2)%p
    inv2=modinverse(2,p)
    x=(M[0][0]-t*inv2)%p
    y=M[0][1]%p
    z=M[1][0]%p
    #assert((x**2+y*z)%p==(T*modinverse(4,p)-1)%p)
    if 4%p==T%p:
        R=np.array([[1,1],[0,1]])%p
        if y==0:
            C=np.array([[0,1],[z,0]])%p
        else:
            yinv=modinverse(y,p)
            C=np.array([[yinv,0],[x*yinv,1]])%p
    else:
        R=np.array([[t*inv2,1],[x**2+y*z,t*inv2]])%p
        if y%p==0:
            xinv=modinverse(x,p)
            if z%p==0:
                C=np.array([[xinv,1],[1,-x]])%p
            else:
                C=np.array([[-z*xinv,1],[0,-x]])%p
        else:
            yinv=modinverse(y,p)
            C=np.array([[yinv,0],[x*yinv,1]])%p
    return R,C
        
    
def stabilizer(M,p):
    """
    Generator for elements of stabilizer of PSL2p class of matrix M by the conjugation action of PGL2p.
    """
    if PSListrivial(M,p):
        for c in range(1,p):
            for d in range(p):
                yield np.array([[0,1],[c,d]])
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    if d!=b*c%p:
                        yield np.array([[1,b],[c,d]])
    if 1!=M[0][1]%p:
        R,C=favoritePSLrep(M,p)
    else:
        R=M
        C=np.array([[1,0],[0,1]])
    t=modtrace(R,p)
    if 4%p==t**2%p:
        for b in range(p):
            yield modconj(np.array([[1,b],[0,1]]),C,p)
    elif 0==t**2%p:
        yield modconj(np.array([[1,0],[0,1]]),C,p)
        yield modconj(np.array([[0,1],[1,0]]),C,p)
        for d in range(p):
            if p-1!=d**2%p:
                yield modconj(np.array([[d,1],[-1,d]]),C,p)
                yield modconj(np.array([[-1,d],[d,1]]),C,p)
    else:
        yield modconj(np.array([[1,0],[0,1]]),C,p)
        Delta=(t**2*modinverse(4,p)-1)%p
        for d in range(p):
            if Delta%p!=d**2%p:
                yield modconj(np.array([[d,1],[Delta,d]]),C,p)


        


        
    
