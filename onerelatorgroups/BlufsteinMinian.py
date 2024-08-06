import networkx as nx
import grouptheory.smallcancellation as sc
import itertools
from fractions import Fraction

# FutureWarning: adjacency_matrix will return a scipy.sparse array instead of a matrix in Networkx 3.0.
"""
Implements the construction of Martin A. Blufstein and Elias G. Minian, Strictly systolic angled complexes and hyperbolicity of one-relator groups, Algebr. Geom. Topol. (in press)  arXiv:1907.06738.

This implementation is only for one-relator groups.
"""


def BlufsteinMinian(relator,Cprimebound=None):
    """
    Decide if the one-relator group defined by the given relator satisfies the hyperbolicity condition of Blufstein Minian that it is Cprime(1/4) and satisfies a condition they call Tprime.
    """
    if Cprimebound is None:
        Cprimebound=sc.Cprimebound([relator])
    return Cprimebound<Fraction(1,4) and BlufsteinMinianTprime(relator)

def BlufsteinMinianTprime(relator):
    """
    Decide if the one-relator group defined with given relator satisfies the Blufstein Minian T-prime condition.

    Want to know that every interior tripod has length strictly less than half the length of the relator.

    >>> BlufsteinMinianTprime('aaaababaBAbbbAb')
    True
    >>> BlufsteinMinianTprime('DCabFEcdBAef')
    False
    >>> BlufsteinMinianTprime([-3,-3,-2,-3,-2,3,-1,-1,-2,-2,-3,1,1])
    False
    >>> BlufsteinMinianTprime('abbAbabABB')
    False
    """
    R=sc.parseinputwords([relator],asrelatorlist=True)[0]
    Rinv=sc.inverse(R)
    def longestcommonprefix(list1,list2):
        matchingprefixlength=0
        for i in range(min(len(list1),len(list2))):
            if list1[i]!=list2[i]:
                return matchingprefixlength
            else:
                matchingprefixlength+=1
        return matchingprefixlength
    def overlap(a,b):
        if a[1]==1 and b[1]==1:
            return longestcommonprefix((R+R)[a[0]+1:],(Rinv+Rinv)[len(R)-1-b[0]:])
        elif a[1]==-1 and b[1]==1:
            return longestcommonprefix((Rinv+Rinv)[a[0]+1:],(Rinv+Rinv)[len(R)-1-b[0]:])
        elif a[1]==1 and b[1]==-1:
            return longestcommonprefix((R+R)[a[0]+1:],(R+R)[len(R)-1-b[0]:])
        elif a[1]==-1 and b[1]==-1:
            return longestcommonprefix((Rinv+Rinv)[a[0]+1:],(R+R)[len(R)-1-b[0]:])
        else:
            raise InputError
    G=sc.simple_Whitehead_graph([R]) # reduced Whitehead graph with at most one edge between vertices 
    A=nx.adjacency_matrix(G)
    nodelist=list(G.nodes())
    B=A@A@A # the point of B is to cut down the number of combinations in the next line. We only look at vertices that do belong to some 3-cycle.
    threecycles={(nodelist[i],nodelist[j],nodelist[k]) for (i,j,k) in itertools.combinations([h for h in range(len(nodelist)) if B[h,h]!=0],3) if A[i,j] and A[j,k] and A[k,i]} # ordered triples of vertices forming 3-cycle in reduced Whitehead graph
    for (i,j,k) in threecycles:
        # find possible indices in R and Rinv that could make an interior tripod with whose vertex has outgoing edges labelled i,j,k
        possiblefirstindex=[(h,1) for h in range(len(R)) if (R+R)[h:h+2]==(-i,j)]+[(h,-1) for h in range(len(R)) if (Rinv+Rinv)[h:h+2]==(-i,j)] # (h,1) means R[h]=-i and R[h+1]=j, (h,-1) means Rinv[h]=-i and Rinv[h+1]=j; these are all possible turns in the relator that give an edge from i to j in the Whithead graph
        possiblesecondindex=[(h,1) for h in range(len(R)) if (R+R)[h:h+2]==(-j,k)]+[(h,-1) for h in range(len(R)) if (Rinv+Rinv)[h:h+2]==(-j,k)]
        possiblethirdindex=[(h,1) for h in range(len(R)) if (R+R)[h:h+2]==(-k,i)]+[(h,-1) for h in range(len(R)) if (Rinv+Rinv)[h:h+2]==(-k,i)]
        for (f,s,t) in itertools.product(possiblefirstindex,possiblesecondindex,possiblethirdindex):
            # for each possible interior tripod, compute its length
            tripodlength=overlap(f,s)+overlap(s,t)+overlap(t,f)
            if 2*tripodlength>=len(R):
                return False
    return True
                   
    










# in terminal, do
# python BlufsteinMinian.py
# to run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()
