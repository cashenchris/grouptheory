import numpy as np
from scipy.optimize import linprog

def sapirspakulovacondition(relator):
    """
    Checks Sapir-Spakulova condition to see if one relator group with given relator is embeddable in an ascending HNN extension of a free group (hence is residually finite). 

    Returns True, False, or None if numerical accuracy is insufficient to decide.
    Input relator should be a list of nonzero integers corresponding to a freely reduced word in the free group where positive i represents the i-th basis element and -i represents its inverse.

    >>> sapirspakulovacondition([1,2,-1,-2,-2])
    True
    >>> sapirspakulovacondition([1,2,2,-1,-2,-2,-2])
    False
    """
    tr=traces(relator)
    tracevector=tr[0][-1]
    if np.allclose(tracevector,np.zeros(tracevector.shape)):
        raise ValueError("Sapir-Spakulova condition is not valid for relators  in the commutator subgroup.")
    ambiguous=False
    for i in range(1,1+max(abs(x) for x in relator)):
        result=touching(tr[i],tracevector)
        if result == True:
            return True
        elif result == None: # numerically ambiguous if there is a touching hyperplane for this i, but continue checking other i in case one of them gives a definitive positive answer.
            ambiguous=True
    if ambiguous: # We didn't find any i for which there is definitively a touching hyperplane, but for at least one of them there was a hyperplane that we couldn't tell  if it was touching or not.
        return None
    else:
        return False
        
        

def traces(relator):
    """
    computes trace in Z^r of word 'relator' in F_r. Returns a dict tr where tr[0] is list of vectors in trace and for i>0 tr[i] is an even length list where pairs of vectors are endpoints of the i-edges in the trace.
    """
    thetraces=dict()
    rank=max(abs(x) for x in relator)
    for i in range(rank+1):
        thetraces[i]=list()
    currentvertex=np.zeros(rank)
    thetraces[0].append(currentvertex)
    for letter in relator:
        idv=np.zeros(rank)
        idv[abs(letter)-1]=1
        nextvertex=currentvertex+np.sign(letter)*idv
        thetraces[0].append(nextvertex)
        thetraces[abs(letter)].append(currentvertex)
        thetraces[abs(letter)].append(nextvertex)
        currentvertex=nextvertex
    return thetraces

def simple_edges(vertexlist):
    """
    Take an i-trace and return list of simple edges.
    """
    remaining=[(vertexlist[2*i],vertexlist[2*i+1]) for i in range(len(vertexlist)/2)]
    simple=[]
    while remaining:
        thisedge=remaining.pop()
        matching=[i for i in range(len(remaining)-1,-1,-1) if same_edge(thisedge,remaining[i])]
        if not matching:
            simple.append(thisedge)
        else:
            for i in matching:
                del remaining[i]
    return simple

def simple_duplicate_vertices(vertexlist):
    """
    Given a list of points return a list of indices of those points that appear only once and a list of indices with one entry for each point that appears multiple times.
    """
    knownsimple=[]
    knownduplicate=[]
    remainingindices=range(len(vertexlist))
    while remainingindices:
        thisindex=remainingindices.pop()
        matches=[i for i in remainingindices if np.allclose(vertexlist[i],vertexlist[thisindex])]
        if matches:
            knownduplicate.append(thisindex)
            for i in matches:
                remainingindices.remove(i)
        else:
            knownsimple.append(thisindex)
    return knownsimple,knownduplicate

def unique_vertices(vertexlist):
    """
    Take a list of vertices and return a list with one entry for each distinct element of the input list.
    """
    simp,dup=simple_duplicate_vertices(vertexlist)
    return [vertexlist[i] for i in simp+dup]
            

def same_edge(vertexpair1,vertexpair2):
    """
    Given two pairs of vertices, decide if they represent the same edge.
    """
    if (np.allclose(vertexpair1[0],vertexpair2[0]) and np.allclose(vertexpair1[1],vertexpair2[1])) or (np.allclose(vertexpair1[0],vertexpair2[1]) and np.allclose(vertexpair1[1],vertexpair2[0])):
        return True
    else:
        return False
        


def in_hull(thepoint,manypoints):
    """
    Check if a given point is cointained in the convex hull of a collection of points. Return True if point is definitely in the hull, False if it is definitely not, and None if not sure.
    """
    c = np.zeros(len(manypoints))
    A = np.r_[np.array(manypoints).T,np.ones((1,len(manypoints)))]
    b = np.r_[thepoint, np.ones(1)]
    lp = linprog(c, A_eq=A, b_eq=b)
    if lp.success==True and lp.status==0:
        return True
    elif lp.success==False and lp.status==2:
        return False
    else:
        return None


def touching(vertexlist,tracevector):
    tvnormsquared=np.linalg.norm(tracevector)**2
    simplevertexindices=simple_duplicate_vertices(vertexlist)[0]
    ambiguous=False
    # first check if there is a hyperplane touching a simple vertex
    # this requires that there is a hyperpland through simple vertex s parallel to tracevector such that all vertices are on one side and only given s is in the hyperplane
    # translate everything by -s and look for hyperplane through origin containing tracevector, not contining any other of the given points, and such that all given points are on one side.
    # collapse a dimension by taking projection to orthogonal complement of tracevector.
    # if origin is in image of vertexlist under this projection then every potential hyperplane contains more than just the simple vertex s, so no touching at this s.
    # if not, then check if origin is in the convex hull of the projected points. If no, then there is a hyperplane through origin in perp(tracevector) with all projected vertices on one side. Product of this with tracevector, then translated by s, is the desired touching hyperplane.
    for i in simplevertexindices:
        s=vertexlist[i]
        projectedverts=[v-s-(np.dot(v-s,tracevector)/tvnormsquared)*tracevector for v in vertexlist if not np.allclose(v,s)]
        if any(np.allclose(np.zeros(v.shape),v) for v in projectedverts):
            continue        
        inhull=in_hull(np.zeros(tracevector.shape),unique_vertices(projectedverts))
        if inhull==False:
            return True
        elif inhull==None:
            ambiguous=True
    simpleedges=simple_edges(vertexlist)
    # Second check if there is hyperplane touching a simple edge.
    # If we got this far then there is no hyperplane touching only a single simple vertex.
    # The only way there could be hyperplane touching only a simple edge but not one touching a simple vertex is if the direction of the simple edge matches tracevector.
    for e in simpleedges:
        s=e[0]
        t=e[1]
        pt=t-s-(np.dot(t-s,tracevector)/tvnormsquared)*tracevector
        if not np.allclose(pt,np.zeros(pt.shape)):
            continue
        projectedverts=[v-s-(np.dot(v-s,tracevector)/tvnormsquared)*tracevector for v in vertexlist if not (np.allclose(v,s) or np.allclose(v,t))]
        if any(np.allclose(np.zeros(v.shape),v) for v in projectedverts):
            continue
        inhull= in_hull(np.zeros(tracevector.shape),unique_vertices(projectedverts)) 
        if inhull==False:
            return True
        elif inhull==None:
            ambiguous=True
    if ambiguous:
        return None
    else:
        return False

    
    




# in terminal, do
# python SapirSpakulova.py
# to run doctests




if __name__ == "__main__":
    import doctest
    doctest.testmod()
