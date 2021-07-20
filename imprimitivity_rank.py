import networkx as nx
import freegroup as fg
import whiteheadgraph as wg #  need this for primitivity check


    
def imprimitivityrank(theword,precomputedWsubgroups=None):
    """
    Given a list or tuple of non-zero integers interpreted as a word in a free group, returns the minimal rank of a subgroup contaiing theword as an imprimitive element. Returns float('inf') if theword is primitve.
    If precomputedWsubgroups=None then Wsubgroups are computed and their rank is returned. Otherwise, the rank of the precomputedWsubgroups is returned.

    >>> imprimitivityrank([1])
    inf
    >>> imprimitivityrank([1,1,1])
    1
    >>> imprimitivityrank('abAB')
    2
    >>> imprimitivityrank([1,1,2,2,3,3])
    3
    """
    F,w=fg.parseinputword(theword)
    if precomputedWsubgroups is None:
        graphs=constructgraphs(w.letters)
    else:
        graphs=precomputedWsubgroups
    if graphs:
        return graphrank(graphs[0])
    else:
        return float('inf')


def Wsubgroups(theword):
    """
    Given a list or tuple of non-zero integers interpreted as a word in a free group, returns a list of Stallings graphs representing the maximal subgroups of minimal rank contaiing theword as an imprimitive element.
    Returns an empty list if theword is primitive.
    """
    F,w=fg.parseinputword(theword)
    graphs=constructgraphs(w.letters)
    maximalgraphs=maximalelements(graphs)
    return maximalgraphs

def parabolic(theword):
    """
    Given a word in a free group, return Wsubgroup containing it and expression of theword as a subword of that subgroup.

    Raises an assertion error if there is not a unique Wsubgroup. Conjecturally this should not happen.

    >>> parabolic([1,2,-1,-2,3,2,1,-2,-1,-3])
    (< c, abAB >, [2, 1, -2, -1])
    """
    F,w=fg.parseinputword(theword)
    graphs=Wsubgroups(w)
    if not graphs:
        return None,None
    else:
        assert(len(graphs)==1) # Conjecturally this is always true. It would be interesting to know a counterexample.
        P=fg.FGSubgroupOfFree(F,[],graph=fg.StallingsGraph(graph=graphs[0]))
        Pword=P.restrict_word(w)
        return P,Pword

def higher_irank(theword):
    """
    Return True if the imprimitivity rank of theword is greater than 2.
    """
    # Start constructing graphs, but stop once it becomes clear that no rank 2 graph will suffice.
    F,w=fg.parseinputword(theword)
    r,p=F.max_root(w,uptoconjugacy=True,withpower=True)
    if p>1: # theword is a proper power, so irank=1
        return False
    if wg.is_primitive(F,w,guaranteenonpower=True): # theword is primitve, so irank=float('inf')
        return True
    # otherwise check if it has irank=2
    graphs=constructgraphs(w.letters,upper_rank_bound=2,guarantee_nonprimitive=True)
    if graphs: # irank <=2
        return False
    else: # 2< irank < float('inf')
        return True



def constructgraphs(theword,upper_rank_bound=None,guarantee_nonprimitive=False,notetrouble=False):
    """
    Given a list or tuple of non-zero integers interpreted as a word in a free group, returns a list of Stallings graphs representing the subgroups of minimal rank contaiing theword as an imprimitive element.

    Do not continue search branch if current candidate graph has rank exceeding upper_rank_bound.
    Returns an empty list if no suitable graph are found. This means imprimitivity rank is greater than given upper_rank_bound, or, if upper_rank_bound=None, that the word is primitive and the imprimitivity rank is float('inf').

    Set guarantee_nonprimtive=True if theword is known to be nonprimitive, to avoid repeating check.
    """
    # notetrouble is for identifying words where the construction gets "more difficult". This is for separate experments.
    if notetrouble:
        Trouble=False
    G=nx.MultiDiGraph() 
    G.add_node(0) # base vertex is named 0
    if theword:
        rank=max([abs(x) for x in theword]) # rank of the ambient free group containing the word
        if upper_rank_bound is None:
            bestrank=rank # upper bound for the rank of subgroup we are looking for, unless the theword is primitive
        else:
            bestrank=min(rank,upper_rank_bound)
        maxedges=dict([(i,[abs(x) for x in theword].count(i)/2) for i in range(1,1+rank)]) # a non-primitive word must use every edge in the Stallings graph G at least twice, so the max number of distinct edges in G labeled i is at most half the number of appearances of +-i in theword
    else: # if theword is empty return trivial graph
        if notetrouble:
            return [G],Trouble
        return [G]
    F=fg.FGFreeGroup(numgens=rank)
    r,p=F.max_root(F.word(theword),uptoconjugacy=True,withpower=True)
    if p != 1: # theword is a proper power, graph is cycle labelled r
        for i in range(len(r)):
            G.add_edge(i,(i+1)%len(r),label=r.letters[i])
        if notetrouble:
            return [G],Trouble
        return [G]
    else:
        if not guarantee_nonprimitive:
            if wg.is_primitive(F,F.word(theword),guaranteenonpower=True): # if theword is already primitive then it will be primitive in every subgroup
                if notetrouble:
                    return [],Trouble
                return []
    # if we haven't returned yet then theword is nontrivial, nonprimitive, and not a proper power
    workinggraphs=[]
    finishedgraphs=[]
    Rose=nx.MultiDiGraph()
    Rose.add_node(0)
    for i in range(1,rank+1):
        Rose.add_edge(0,0,label=i)
    finishedgraphs.append(Rose)
    workinggraphs.append((G,0,theword))
    # elements of workinggraphs are triples (graph G, active vertex, suffix), where prefix+suffix=theword and there is a path in G from  0 to activevertex labelled by prefix. The next step will be to add/follow an edge from active vertex labelled by the first letter of suffix. Such an edge can have oppositive vertex new or be one of the existing vertices in the graph, provided that adding such an edge does not create an unfolded graph, and that it hasn't exceeded maxedges for that label or bestrank. If these conditions are ok then add (new graph, new vertex, suffix[1:]) to workinggraphs, or, if suffix[1:]=[], to finishedgraphs. 
    while workinggraphs:
        oldg=workinggraphs.pop()
        if graphrank(oldg[0])>bestrank:
            continue
        currentvertex=oldg[1]
        nextlabel=oldg[2][0]
        nextsuffix=oldg[2][1:]
        nextvert=vertexhaslabel(oldg[0],currentvertex,nextlabel,returnopvert=True)
        if nextvert is not None: # there is already an incident edge with correct label. Follow it.
            if not nextsuffix: # we have exhausted theword
                if nextvert==0: # if we have returned to the baseponit theword belongs to this subgroup. Otherwise throw this graph away.
                    thisrank=graphrank(oldg[0])
                    if thisrank<=bestrank:# if graph has rank at most current best then it might be a candidate. If not, discard.
                        if everyedgetwice(oldg[0],theword):# if theword traverses every edge of the graph twice then it might be a candidate. If not, theword is either contained in a proper free factor or is primitive, so this graph is not a candidate for minimal rank subgroup containing word as imprimitive; discard.
                            if nonprimitive(oldg[0],theword):# this is definitive check that theword defines imprimitive element 
                                finishedgraphs.append(oldg[0].copy())
                                bestrank=thisrank
                            else:#quickcheck for imprimitivity passed, but real check failed. Thsese words might be trouble.
                                if notetrouble:
                                    Trouble=True
            else: # theword is not exhausted. iterate.
                workinggraphs.append((oldg[0].copy(),nextvert,nextsuffix))
        else: #there is not already an incident edge with the correct label. We must add a new edge. 
            # adding a new edge will mean that any completed successor of this graph will have strictly higher rank.
            # if that rank would be too high then we can stop now
            if graphrank(oldg[0])==bestrank or graphrank(oldg[0])>=rank-1: 
                continue
            if len([e for e in oldg[0].edges(keys=True,data=True) if abs(e[3]['label'])==abs(nextlabel)])<maxedges[abs(nextlabel)]: #we haven't yet exceeded max number of edges with this label, so can try to add a new one
                if not nextsuffix: #if we are out of letters then the only way to make our immersed cycle is to now connect back to the basepoint
                    if vertexhaslabel(oldg[0],0,-nextlabel):
                        pass # base vertex already has an edge with the desired label, so adding another would make unfolded graph, discard
                    else:# basepoint does not already have an incident edge with this lable, so ok to make one
                        newg=nx.MultiDiGraph(oldg[0])
                        newg.add_edge(currentvertex,0,label=nextlabel)
                        newrank=graphrank(newg)
                        if newrank<=bestrank:
                            if everyedgetwice(newg,theword):
                                if nonprimitive(newg,theword):
                                    finishedgraphs.append(newg)
                                    bestrank=newrank
                                else:
                                    if notetrouble:
                                        Trouble=True
                else: # we are not out of leffers, so can add a new edge going to any available vertex, or to a new vertex
                    # first we add an edge with a new opp vertex
                    nextvertex=len(oldg[0]) # label for a new vertex is number of vertices in the current graph
                    newg=nx.MultiDiGraph(oldg[0])
                    newg.add_edge(currentvertex,nextvertex,label=nextlabel)
                    if graphrank(newg)<=bestrank:
                        workinggraphs.append((newg,nextvertex,nextsuffix))
                    # next we try adding an edge connecting to one of the existing vertices, but only in places such that the graph remains folded.
                    for nextvertex in oldg[0]:
                        if vertexhaslabel(oldg[0],nextvertex,-nextlabel):
                            pass # this vertex already has an edge with the desired label, skip it
                        else:
                            newg=nx.MultiDiGraph(oldg[0])
                            newg.add_edge(currentvertex,nextvertex,label=nextlabel)
                            if graphrank(newg)<=bestrank:
                                workinggraphs.append((newg,nextvertex,nextsuffix))
    if notetrouble:
            return [G for G in finishedgraphs if graphrank(G)<=bestrank],Trouble
    return [G for G in finishedgraphs if graphrank(G)<=bestrank]


##########---------auxiliary functions

def vertexhaslabel(thegraph,thevertex,thelabel,returnopvert=False):
    """
    Check if thevertex already has an outgoing edge labelled with thelabel.
    returnopvert=False then return bool.
    returnopvert=True then return the oppositve vertex on the edge labelled with thelabel, if such an edge exists, or None if no such edge.
    """
    for e in thegraph.out_edges(thevertex,data=True,keys=True):
        if e[3]['label']==thelabel:
            if returnopvert:
                return e[1]
            else:
                return True
    else:
        for e in thegraph.in_edges(thevertex,data=True,keys=True):
            if e[3]['label']==-thelabel:
                if returnopvert:
                    return e[0]
                else:
                    return True
        else:
            if returnopvert:
                return None
            else:
                return False

def graphrank(thegraph):
    """
    Returns the rank of a connected graph.
    """
    return len(thegraph.edges())-len(thegraph)+1

def spanningtree(G):
    """
    Return a list of edges of a graph G that give a spanning tree.
    """
    basepoint=list(G)[0]
    seen=set([basepoint])
    newvertices=[basepoint]
    treeedges=[]
    while newvertices:
        thisvertex=newvertices.pop()
        for e in G.out_edges(thisvertex,keys=True):
            if e[1] in seen:
                pass
            else:
                seen.add(e[1])
                treeedges.append(e)
                newvertices.append(e[1])
        for e in G.in_edges(thisvertex,keys=True):
            if e[0] in seen:
                pass
            else:
                seen.add(e[0])
                treeedges.append(e)
                newvertices.append(e[0])
    return treeedges

def freebasis(G):
    """
    The complement of a spanning tree.
    """
    T=spanningtree(G)
    return [e for e in G.edges(keys=True) if e not in T]

def wordexpressedinfreebasis(thegraph,thebasepoint,theword,thefreebasis):
    """
    Given a labelled, based graph,  a list of edges that form the complement of a spanning tree and a word that is the concatenation of labels read on a loop based at thebasepoint, return the expression of that loop in terms of the given freebasis for the fundamental group of thegraph with respect to thebasepoint.
    Returns a list of non-zero integers where i corresponds to the edge thefreebasis[i-1] and -i corresponds to the edge thefreebasis[i-1] traversed backwards.
    """
    newexpression=[]
    currentvertex=thebasepoint
    for theletter in theword:
        for e in thegraph.out_edges(currentvertex,keys=True,data=True):
            if e[3]['label']==theletter:
                if e[:3] in thefreebasis:
                    newexpression.append(1+thefreebasis.index(e[:3]))
                currentvertex=e[1]
                break
        else:
            for e in thegraph.in_edges(currentvertex,keys=True,data=True):
                if e[3]['label']==-theletter:
                    if e[:3] in thefreebasis:
                        newexpression.append(-1-thefreebasis.index(e[:3]))
                    currentvertex=e[0]
                    break
            else:
                raise KeyError
    if currentvertex==thebasepoint:
        return newexpression
    else:
        raise KeyError

def everyedgetwice(thegraph,theword):
    """
    Checks if theword defines a loop in the graph based at 0 that traverses every edge at least twice.
    """
    # if not, the element defined by theword in the subgroup defined by thegraph is either primitive or containedin a proper free factor
    edgecount=dict()
    for e in thegraph.edges(keys=True):
        edgecount[e[2]]=0
    currentvertex=0
    for theletter in theword:
        for e in thegraph.out_edges(currentvertex,keys=True,data=True):
            if e[3]['label']==theletter:
                edgecount[e[2]]+=1
                currentvertex=e[1]
                break
        else:
            for e in thegraph.in_edges(currentvertex,keys=True,data=True):
                if e[3]['label']==-theletter:
                    edgecount[e[2]]+=1
                    currentvertex=e[0]
                    break
            else:
                raise KeyError
    assert(currentvertex==0)
    return all(edgecount[x]>=2 for x in edgecount)
                
    

def nonprimitive(thegraph,theword):
    """
    thegraph is a labelled graph with basepoint 0.
    theword is a word in the lables of the graph that labells a loop based at 0.
    Decides whether theword defines a nonprimitive element of the freegroup defined by the fundamental group of the graph.
    This function assumes that theword is already guaranteed to be not a proper power. 
    """
    fb=freebasis(thegraph)
    theletters=wordexpressedinfreebasis(thegraph,0,theword,fb)
    F=fg.FGFreeGroup(numgens=len(fb))
    w=F.word(theletters)
    return not wg.is_primitive(F,w,guaranteenonpower=True)
                    

    
def contains(G,Gbase,H,Hbase):
    """
    Return true if the subgroup defined by G contains the subgroup defined by H, where G and H are folded Stallings graphs.
    """
    vertmap=dict()
    edgemap=dict()
    vertmap[Hbase]=Gbase
    newvertices=[Hbase]
    while newvertices:
        currentvertex=newvertices.pop()
        for e in H.out_edges(currentvertex,data=True,keys=True):
            if tuple(e[:3]) in edgemap:
                continue
            if e[1] not in vertmap:
                newvertices.append(e[1])
            for f in G.out_edges(vertmap[currentvertex],data=True,keys=True):
                if f[3]['label']==e[3]['label']:
                    if e[1] in vertmap:
                        if vertmap[e[1]]!=f[1]:
                            return False
                    edgemap[tuple(e[:3])]=tuple(f[:3])
                    vertmap[e[1]]=f[1]
                    break
            else:
                for f in G.in_edges(vertmap[currentvertex],data=True,keys=True):
                    if f[3]['label']==-e[3]['label']:
                        if e[1] in vertmap:
                            if vertmap[e[1]]!=f[0]:
                                return False
                        edgemap[tuple(e[:3])]=tuple(f[:3])
                        vertmap[e[1]]=f[0]
                        break
                else:
                    return False
    return True
                

def maximalelements(graphs):
    """
    Given a list of Stallings graphs, return the ones representing subgroups that are maximal with respect to inclusion among elements in the input list.
    """
    maximalelements=[]
    for i in range(len(graphs)):
        for j in range(len(graphs)):
            if i==j:
                continue
            if contains(graphs[j],0,graphs[i],0):
                break
        else:
            maximalelements.append(graphs[i])
    return maximalelements


















if __name__ == "__main__":
    import doctest
    doctest.testmod()
