from __future__ import division
import sympy as sym


def cnorm2(z):
    return sym.N(sym.simplify(z*sym.conjugate(z)))


def mobius(pairormatrix,p):
    a=pairormatrix[0]
    b=pairormatrix[1]
    return (a*p+b)/(sym.conjugate(b)*p+sym.conjugate(a))

def invmobius(pairormatrix,p):
    a=pairormatrix[0]
    b=pairormatrix[1]
    return mobius((sym.conjugate(a), -b), p)
        
def dist(p,q):
    p=complex(p)
    q=complex(q)
    if any(part==float('inf') for part in [sym.re(p),sym.im(p),sym.re(q),sym.im(q)]):
        return float('inf')
    else:
        dif=p-q
        return sym.acosh(1+2*(sym.re(dif)**2+sym.im(dif)**2)/((1 - sym.re(p)**2 - sym.im(p)**2)*(1- q.real**2 - q.imag**2)))

def direction(p):
    if p==0:
        return float('nan')
    else:
        p=complex(p)
        return p/sym.sqrt(sym.re(p)**2 + sym.im(p)**2)
    

def EuclideanNorm(d):
    """
    Euclidean Norm of point hyperbolic distance d from origin
    """
    x=(sym.cosh(d)-1)/2
    return sym.sqrt(sym.N(x/(1+x)))

def findPoint(p,q,d):
    """
    return the point hyperbolic distance d from p in the direction of q
    """
    p=complex(p)
    q=complex(q)
    if p==q or d==0:
        return p
    else:
        forward=1 if d>=0 else -1
        Q=complex(mobius((1,-p),q))
        newpoint=forward*EuclideanNorm(d)*direction(Q)
        return complex(invmobius((1,-p),newpoint))

def EuclideanBall(Hcenter, Hradius):
    if Hcenter==0:
        return (0,EuclideanNorm(Hradius))
    else:
        z1=findPoint(Hcenter,0,Hradius)
        z2=findPoint(Hcenter,0,-Hradius)
        Ecenter=(z1+z2)/2
        Eradius=abs((z1-z2)/2)
        return (Ecenter,Eradius)

def standardrep(rank):
    """
    representation of free group of rank r as hyperbolic isometries
    """
    rep=dict.fromkeys(range(1,rank+1))
    theta=sym.pi/rank
    outradius=sym.cos(theta/2)/(1+sym.sin(theta/2))
    inradius=outradius*sym.cos(theta/2)
    Hdist=2*dist(0,inradius)
    Edist=EuclideanBall(0,Hdist)[1]
    base=sym.Matrix([[1,sym.cos(theta/2)],[sym.cos(theta/2),1]])
    for i in rep:
        rot=sym.Matrix([[sym.exp(sym.I*(i-1)*theta/2),0],[0,sym.exp(-sym.I*(i-1)*theta/2)]])
        ROT=sym.Matrix([[sym.exp(-sym.I*(i-1)*theta/2),0],[0,sym.exp(sym.I*(i-1)*theta/2)]])
        newM=rot*base*ROT
        rep[i]=sym.Matrix([[sym.N(sym.simplify(newM[0])),sym.N(sym.simplify(newM[1]))],[sym.N(sym.simplify(newM[2])),sym.N(sym.simplify(newM[3]))]])
    return rep

def centerofmass(pointsandmasses):
    """
    return the center of mass of a list of points in the disc
    """
    if len(pointsandmasses)==0:
        return (complex(0) , 0)
    if len(pointsandmasses)==1:
       return(pointsandmasses[0])
    else:
        p1=centerofmass([pointsandmasses[i] for i in range(0,len(pointsandmasses),2)])
        p2=centerofmass([pointsandmasses[i] for i in range(1,len(pointsandmasses),2)])
        if (p1[1]+p2[1])==0:
            return (complex(0) , 0)
        else:
            return (findPoint(p1[0],p2[0],sym.N((p2[1]/(p1[1]+p2[1]))*dist(p1[0],p2[0]))) , p1[1]+p2[1])
        
    
            
    

