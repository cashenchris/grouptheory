import grouptheory.group
import networkx as nx
import grouptheory.freegroups.whiteheadgraph.orderedmultigraph as omg
import grouptheory.freegroups.whiteheadgraph.wgraph as wg
import math
import random
import grouptheory.freegroups.whiteheadgraph.draw.hyperbolicdisc as hdisc
import grouptheory.freegroups.whiteheadgraph.draw.hyperbolicupperhalfplane as hup
import sympy as sym

def Edist(p,q):
    return sym.sqrt((sym.re(p)-sym.re(q))**2+(sym.im(p)-sym.im(q))**2)

def tooclose(point, pointset, tolerance):
    """
    returns True if point is within distance tolerance of some point in pointset
    """
    return any(sym.N(Edist(point,p))<=tolerance for p in pointset)

def EuclideanVertPos(vert, rank,representation={}):
    """
    Calculuate coordinates of vertex in Euclidean plane.
    representation is dictionary containing images of the generators of free group of
    the given rank in the isometry group of the plane.
    """
    # default representation is i->(v->v+(cos((i-1)pi/rank),sin((i-1)pi/rank)))
    if type(vert)!=tuple:
        vert=(vert,)
    abelianization=dict.fromkeys(range(1,rank+1),0)
    for i in range(len(vert)-1):
        abelianization[abs(vert[i])]+=vert[i]/abs(vert[i])
    theta=sym.pi/rank
    majorposition=0
    for i in abelianization:
        majorposition+=abelianization[i]*sym.exp(sym.I*(i-1)*theta)
    minorposition=complex((vert[-1]/abs(vert[-1]))*.25*sym.exp(sym.I*theta*(abs(vert[-1])-1)))
    return complex(majorposition+minorposition)

        
def DiscVertPos(vert, rank, representation=None):
    """
    Calculuate coordinates of vertex in Poincare disc.
    representation is dictionary containing images of the generators of free group of
    the given rank in the isometry group of the plane.
    """
    if representation==None:
        representation=hdisc.standardrep(rank)
    if type(vert)!=tuple:
        vert=(vert,)
    gens=list(vert)
    last=gens.pop()
    if last<0:
        position=.33*hdisc.invmobius(representation[-last],0)
    else:
        position=.33*hdisc.mobius(representation[last],0)
    while gens!=[]:
        nextgen=gens.pop()
        if nextgen<0:
            position=hdisc.invmobius(representation[-nextgen],position)
        else:
            position=hdisc.mobius(representation[nextgen],position)
    return sym.N(position)

        

def UpperVertPos(vert, rank, representation=None):
    """
    Calculuate coordinates of vertex in upper half space model of hyperbolic plane.
    representation is dictionary containing images of the generators of free group of
    the given rank in the isometry group of the plane.
    """
    if representation==None:
        representation=hup.standardrep(rank)
    if type(vert)!=tuple:
        vert=(vert,)
    gens=list(vert)
    last=gens.pop()
    if last<0:
        position=hup.findPoint(sym.I,hup.invmobius(representation[-last],sym.I),.6)
    else:
        position=hup.findPoint(sym.I,hup.mobius(representation[last],sym.I),.6)
    while gens!=[]:
        nextgen=gens.pop()
        if nextgen<0:
            position=hup.invmobius(representation[-nextgen],position)
        else:
            position=hup.mobius(representation[nextgen],position)
    return complex(position)

def draw(G,layout=None,rank=None, representation=None,center=None):
    pos={}
    if rank==None:
        rank=G.rank
    if layout=="Euclidean":
        pos=dict((vert,EuclideanVertPos(vert, rank, representation)) for vert in G)
        # jiggle vertices to make positions unique
        positions=set([])
        for vert in pos:
            while tooclose(pos[vert], positions,.1):
                theta=2*math.pi*random.random()
                pos[vert]=sym.N(pos[vert]+.2*sym.exp(theta*sym.I))
            positions.update(set([pos[vert]]))
            pos[vert]=[complex(pos[vert]).real,complex(pos[vert]).imag]
    elif layout=="upper":
        if representation==None:
            representation=hup.standardrep(rank)
        pos=dict((vert,UpperVertPos(vert, rank, representation)) for vert in G)
        if center==None:
            vertsandmasses=[(pos[p],1) for p in pos]
            center=hup.centerofmass(vertsandmasses)[0]
        if center==sym.I:
            for vert in pos:
                pos[vert]=[complex(pos[vert]).real,complex(pos[vert]).imag]
        else:    
            for vert in pos:
                newpos=hup.mobius(sym.Matrix([[1,-center.real],[0,center.imag]]),pos[vert])
                pos[vert]=[complex(newpos).real,complex(newpos).imag]
            
    elif layout=="disc":
        if representation==None:
            representation=hdisc.standardrep(rank)
        pos=dict((vert,DiscVertPos(vert, rank, representation)) for vert in G)
        if center==None:
            vertsandmasses=[(pos[p],1) for p in pos]
            center=hdisc.centerofmass(vertsandmasses)[0]
        if center==0:
            for vert in pos:
                pos[vert]=[complex(pos[vert]).real,complex(pos[vert]).imag]
        else:
            for vert in pos:
                newpos=hdisc.invmobius((1,center),pos[vert])
                pos[vert]=[complex(newpos).real,complex(newpos).imag]
    if pos=={}:
        # if no positions have been set default to nx.spring_layout
        pos=nx.spring_layout(G)
    nx.draw(G,pos)
