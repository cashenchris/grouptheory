import networkx as nx
from smallcancellation import *
import itertools

# relator as a string with capitalization as inverse
# assume relator is freely and cyclically reduced and not a proper power and every generator appears at least twice, including possibly as inverse, otherwise the relator is primitive.
# segment=(v,e,l) signifies starting vertex v, direction e +1 or -1, and length l of a cyclic subword of w or its inverse W

def inverse(w):
    return w[::-1].swapcase()

def subword(w,segment):
    (v,e,l)=segment
    if e==1:
        return (w+w)[v:v+l]
    if e==-1:
        W=inverse(w)
        return (W+W)[-v%len(w):(-v%len(w))+l]

def pieces(w):
    if any(w[i]==w[(i+1)%len(w)].swapcase() for i in range(len(w))):
        raise ValueError("Word is not freely and cyclically reduced.")
    theletters=dict()
    for c in w:
        theletters[c.lower()]=theletters.setdefault(c.lower(),0)+1
    if any(theletters[c]<2 for c in theletters):
        raise ValueError("Word is primitive.")
    if any(w==w[:len(w)//p]*p for p in range(2,len(w)+1)):
        raise ValueError("Word is a proper power.")
    thepieces=dict()
    for segment in itertools.product(range(len(w)),[1,-1],range(1,len(w))):
        thesubword=subword(w,segment)
        thepieces.setdefault(thesubword,set()).add(segment)
    for k in [k for k in thepieces]:
        if len(thepieces[k])<2:
            del thepieces[k]
    return thepieces
    
def piecesegments(w,precomputed_pieces=None):
    if precomputed_pieces is None:
        thepieces=pieces(w)
    else:
        thepieces=precomputed_pieces
    return set().union(*[thepieces[k] for k in thepieces])

def linkgraph(w,precomputed_piecesegments=None):
    if precomputed_piecesegments is None:
        thepiecesegments=piecesegments(w)
    else:
        thepiecesegments=precomputed_piecesegments
    G=nx.DiGraph()
    nodes={(v,e) for (v,e,l) in thepiecesegments}
    G.add_nodes_from(nodes)
    for (u,e_u),(v,e_v) in itertools.product(nodes,repeat=2):
        if u==v:
            continue
        if subword(w,(u,e_u,1))==subword(w,(v,-e_v,1)):
            G.add_edge((u,e_u),(v,e_v))
    return G

def tauprime(w,strict=False,verbose=False):
    """
    >>> tauprime('aaabbbbaaabbbbAAABBBBAAABBBB')
    True
    >>> tauprime('aaabbbbaaabbbbAAABBBBAAABBBB',strict=True)
    False
    """
    thepieces=pieces(w)
    thepiecesegments=piecesegments(w,precomputed_pieces=thepieces)
    def maxpiece(cornerwithorientation):
        c,o=cornerwithorientation
        maxpiecelength=max((l for (v,e,l) in thepiecesegments if v==c and e==o))
        return subword(w,(c,o,maxpiecelength))
    G=linkgraph(w,precomputed_piecesegments=thepiecesegments)
    def specialcycles():
        corners=sorted(list({v for (v,e) in G}))
        for v in corners:
            thiscycle=[(v,1),]
            nextcycle=extendcycle(thiscycle)
            while len(nextcycle)>1:
                if nextcycle[0]==nextcycle[-1]: # its a closed cycle now
                    if len(nextcycle)>3: # if nextcycle has length three it only has two distinct entries. We don't want valence 2 vertices.
                        yield nextcycle[:-1]
                    nextcycle=incrementcycle(nextcycle,len(nextcycle)-1)
                else:
                    nextcycle=extendcycle(nextcycle)
    def extendcycle(thiscycle):
        firstvertex=thiscycle[0]
        lastvertex=thiscycle[-1]
        if firstvertex in {e[1] for e in G.out_edges(lastvertex)}: # if possible to close the cycle, do it
            return thiscycle+[firstvertex,]
        potentialnextvertices=sorted([e[1] for e in G.out_edges(lastvertex) if e[1][0]>firstvertex[0] and e[1][0] not in {c[0] for c in thiscycle}])
        candidateindex=0
        while candidateindex<len(potentialnextvertices):
            nextvertex=potentialnextvertices[candidateindex]
            candidateindex+=1
            nextvertexpiece=maxpiece(nextvertex)
            #originalpiece=maxpiece((firstvertex[0],-1))
            #if common_prefix_length(originalpiece,nextvertexpiece):
            # the new thing we are adding results in an outgoing edge in the link with the same initial label as the frst edge. Usually we can exclude such edges via an inductive argument.
            #    if len(originalpiece)+len(nextvertexpiece)-2*common_prefix_length(originalpiece,nextvertexpiece)<=len(w):
            #        continue
            if len(thiscycle)>1:
                interiorpieces=[]
                for i in range(len(thiscycle)-1):
                    leftpiece=maxpiece(thiscycle[i])
                    rightpiece=maxpiece((thiscycle[i+1][0],-thiscycle[i+1][1]))
                    interiorpieces.append(leftpiece[:common_prefix_length(leftpiece,rightpiece)])
                if any(common_prefix_length(interiorpiece,nextvertexpiece)>0 and len(interiorpiece)+len(nextvertexpiece)-2*common_prefix_length(interiorpiece,nextvertexpiece)<=len(w) for interiorpiece in interiorpieces):
                    continue
            # if we got here we didn't continue, so no surgery simplification of the link
            return thiscycle+[nextvertex,]
        # no valid way to extend
        if len(thiscycle)==1:
            return thiscycle
        else:
            return incrementcycle(thiscycle,len(thiscycle)-1)
    def incrementcycle(thiscycle,index):
        if index==0:
            return thiscycle[:1]
        firstvertex=thiscycle[0]
        lastvertex=thiscycle[index-1]
        potentialnextvertices=sorted([e[1] for e in G.out_edges(lastvertex) if e[1]>thiscycle[index] and e[1][0] not in {c[0] for c in thiscycle[:index]}])
        candidateindex=0
        while candidateindex<len(potentialnextvertices):
            nextvertex=potentialnextvertices[candidateindex]
            candidateindex+=1
            nextvertexpiece=maxpiece(nextvertex)
            #originalpiece=maxpiece((firstvertex[0],-1))
            #if common_prefix_length(originalpiece,nextvertexpiece):
            # the new thing we are adding results in an outgoing edge in the link with the same initial label as the frst edge. Usually we can exclude such edges via an inductive argument.
            #    if len(originalpiece)+len(nextvertexpiece)-2*common_prefix_length(originalpiece,nextvertexpiece)<=len(w):
            #        continue
            if len(thiscycle)>1:
                interiorpieces=[]
                for i in range(len(thiscycle)-1):
                    leftpiece=maxpiece(thiscycle[i])
                    rightpiece=maxpiece((thiscycle[i+1][0],-thiscycle[i+1][1]))
                    interiorpieces.append(leftpiece[:common_prefix_length(leftpiece,rightpiece)])
                if any(common_prefix_length(interiorpiece,nextvertexpiece) and len(interiorpiece)+len(nextvertexpiece)-2*common_prefix_length(interiorpiece,nextvertexpiece)<=len(w) for interiorpiece in interiorpieces):
                    continue
            return thiscycle[:index]+[nextvertex,]
        #no valid way to increment this index
        return incrementcycle(thiscycle,index-1)

    for cycle in specialcycles(): 
        if verbose:
            print(cycle)
        thesepieces=[]
        for i in range(len(cycle)):
            leftpiece=maxpiece(cycle[i-1])
            rightpiece=maxpiece((cycle[i][0],-cycle[i][1]))
            thesepieces.append(leftpiece[:common_prefix_length(leftpiece,rightpiece)])
        if strict:
            if 2*sum(len(p) for p in thesepieces) >=(len(cycle)-2)*len(w):
                if verbose:
                    print(thesepieces)
                return False
        else:
            if 2*sum(len(p) for p in thesepieces) >(len(cycle)-2)*len(w):
                if verbose:
                    print(thesepieces)
                return False
    return True



    
    
        
def common_prefix_length(stringone,stringtwo):
    if len(stringone)==0 or len(stringtwo)==0:
        return 0
    for L in range(min(len(stringone),len(stringtwo))):
        if stringone[L]!=stringtwo[L]:
            return L
    return min(len(stringone),len(stringtwo))

            






# in terminal, do
# python BlufsteinMinianCosta.py
# to run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()

        
    
    
            
    
