from __future__ import division
import sympy as sym
import hyperbolicdisc as hdisc


def cnorm2(z):
    return sym.N(sym.simplify(z*sym.conjugate(z)))


def mobius(M,p):
    if sym.re(p)==sym.oo or sym.im(p)==sym.oo:
        if M[2]==0:
            return sym.oo
        else:
            return M[0]/M[2]
    elif p==0:
        if M[3]==0:
            return sym.oo
        else:
            return M[1]/M[3]
    else:
        numerator=M[0]*p+M[1]
        denominator=M[2]*p+M[3]
        if denominator==0:
            return sym.oo
        else:
            return numerator/denominator

def invmobius(M,p):
    return mobius(sym.Matrix([[M[3],-M[1]],[-M[2],M[0]]]), p)
        
def dist(p,q):
    p=complex(p)
    q=complex(q)
    if sym.im(p)==0 or sym.im(q)==0 or any(part==sym.oo for part in [sym.re(p), sym.re(q), sym.im(p), sym.im(q)]):
        return sym.oo
    else:
        return sym.acosh(1+((sym.re(q)-sym.re(p))**2+(sym.im(q)-sym.im(p))**2)/(2*sym.im(p)*sym.im(q)))


def findPoint(p,q,d):
    """
    return the point hyperbolic distance d from p in the direction of q
    """
    p=complex(p)
    q=complex(q)
    if p==q or d==0:
        return p
    # Find matrix T for mobius transform that takes p to I and q to the imaginary axis above I
    else:
        Movep=sym.Matrix([[1,-p.real],[0,p.imag]])
        Q=complex(mobius(Movep,q))
        if Q.real==0 and Q.imag>1:
            Rot=sym.Matrix([[1,0],[0,1]])
        elif Q.real==0 and Q.imag<1:
            Rot=sym.Matrix([[0,1],[-1,0]])
        elif Q.real>0:
            x=Q.real
            y=Q.imag
            a=(x**2+y**2-1)/(2*x)
            r=sym.sqrt(1+a**2)
            Rot=sym.Matrix([[(a+r)/sym.sqrt(1+(a+r)**2), 1/sym.sqrt(1+(a+r)**2)],[-1/sym.sqrt(1+(a+r)**2),(a+r)/sym.sqrt(1+(a+r)**2)]])
        else:
            x=-Q.real/(Q.real**2+Q.imag**2)
            y=Q.imag/(Q.real**2+Q.imag**2)
            a=(x**2+y**2-1)/(2*x)
            r=sym.sqrt(1+a**2)
            Rot=sym.Matrix([[0,1],[-1,0]])*sym.Matrix([[(a+r)/sym.sqrt(1+(a+r)**2), 1/sym.sqrt(1+(a+r)**2)],[-1/sym.sqrt(1+(a+r)**2),(a+r)/sym.sqrt(1+(a+r)**2)]])
        T=Rot*Movep            
        return complex(invmobius(T,sym.exp(d)*sym.I))
        

def EuclideanBall(Hcenter, Hradius):
    if Hradius==0:
        return (Hcenter,0)
    elif Hradius==sym.oo:
        return (sym.oo,sym.oo)
    else:
        Movecenter=sym.Matrix([[1,-sym.re(Hcenter)],[0,sym.im(Hcenter)]])
        z1=invmobius(Movecenter,-sym.ln(Hradius))
        z2=invmobius(Movecenter,sym.ln(Hradius))
        Ecenter=(z1+z2)/2
        Eradius=abs((z1-z2)/2)
        return (Ecenter,Eradius)

def standardrep(rank):
    """
    representation of free group of rank r as hyperbolic isometries
    """
    rep=dict.fromkeys(range(1,rank+1))
    discrep=hdisc.standardrep(rank)
    for i in rep:
        newM=sym.Matrix([[sym.I,-1],[-1,sym.I]])*discrep[i]*sym.Matrix([[sym.I,1],[1,sym.I]])
        rep[i]=sym.Matrix([[sym.N(sym.simplify(newM[0])),sym.N(sym.simplify(newM[1]))],[sym.N(sym.simplify(newM[2])),sym.N(sym.simplify(newM[3]))]])
    return rep

def centerofmass(pointsandmasses):
    """
    return the center of mass of a list of points in the disc
    """
    if len(pointsandmasses)==0:
        return (sym.I , 0)
    if len(pointsandmasses)==1:
       return(pointsandmasses[0])
    else:
        p1=centerofmass([pointsandmasses[i] for i in range(0,len(pointsandmasses),2)])
        p2=centerofmass([pointsandmasses[i] for i in range(1,len(pointsandmasses),2)])
        if (p1[1]+p2[1])==0:
            return (sym.I , 0)
        else:
            return (findPoint(p1[0],p2[0],sym.N((p2[1]/(p1[1]+p2[1]))*dist(p1[0],p2[0]))) , p1[1]+p2[1])
        
    
            
    

