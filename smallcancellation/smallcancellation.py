import networkx as nx
import math
from fractions import Fraction
import itertools

# relatorlist returned by parseinputwords as list of lists of nonzero integers 
# segment=(r,v,e,l) means take the reltor at index r in relatorlist, subword starting at vertex v, direction e +1 or -1, and length l 




def CT(relator_list,quit_at=float('inf'),piece_up_to_automorphism=True,precomputed_piecedict=None,noparse=False):
    """
    Given a list of relators that are freely and cyclically reduced, and independent up to inversion and cyclic permutation, return (p,q) such that given relators define a C(p)-T(q) presentation. 
    >>> CT([[1,2,-1,-2]])
    (4, 4)
    >>> CT(['abABcdCD'])
    (8, 8)
    >>> CT(['abABcdCDefEF'])
    (12, 12)
    >>> CT(['aabbccddeeffgg'])
    (14, 14)
    >>> CT(['ababbabbbabbbbabbbbb'])
    (5, 4)
    >>> CT(['abAB'*8])
    (32, 4)
    >>> CT(['abab'])
    (inf, inf)
    >>> CT(['aa','bbb','cccc','dddd','abcd'])
    (2, 8)
    """
    if noparse:
        rels=relator_list
    else:
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels,piece_up_to_automorphism=piece_up_to_automorphism)
    else:
        thepiecedict=precomputed_piecedict
    return C(rels,precomputed_piecedict=thepiecedict),T(rels,precomputed_piecedict=thepiecedict)

def smallcancellation(relatorlist,theCprimebound=None,noparse=False):
    """
    Check if the relatorlist satisfy any of several small cancellation conditions that guarantee hyperbolicity.

    If theCprimebound has already been computed for this relatorlist can input it to shortcircuit redundant computation.
    
    By Gersten-Short, a Cp-Tq presentation defines a hyperbolic group when (p,q) are (3,7), (4,5), (5,4), or (>=7,3).

    >>> smallcancellation(['abABcdCD']) # C'(1/6)
    True
    >>> smallcancellation(['aCbcABCac']) # C7 and not C'(1/6)
    True
    >>> smallcancellation([[-2, -2, -3, -1, -2, -3, -3, -2, 1, 2, 2, 3, 3]]) # C8 and not C'(1/6)
    True
    >>> smallcancellation(['ababccbAbaBCCB']) # C5-T4 
    True
    >>> smallcancellation(['cacbcbcbcabacbcaba']) # C4-T6
    True
    >>> smallcancellation([[1,2,-1,-2]]) # C4-T4
    False
    >>> smallcancellation([[1,2,1,2,2,1,2,2,2]]) # C3-Tinf
    True
    """
    if noparse:
        rels=relatorlist
    else:
        rels=parseinputwords(relatorlist,asrelatorlist=True)
    if theCprimebound is None:
        theCprimebound=Cprimebound(rels,noparse=True)
    if theCprimebound<Fraction(1,6):
        return True
    (theC,theT)=CT(rels,noparse=True)
    if theC>=7 or (theC>=5 and theT>=4) or  (theC>=4 and theT>=5) or (theC>=3 and theT>=7):
        return True
    else:
        return False

    
def T(relator_list, precomputed_piecedict=None, noparse=False, precomputed_link_graph=None,piece_up_to_automorphism=True):
    """
    Return the minimum degree of an essential interior vertex in a reduced van Kampen diagram. 

    >>> T(['aB','bC','cA'])
    3
    >>> T([[1,2,1,2,2,1,2,2,2]]) 
    4
    >>> T(['aabbcc'])
    6
    >>> T(['aaa','bbb','ccc','abc'])
    6
    >>> T(['aaa','bbb','ccc']) # no pieces
    inf
    >>> T(['abABcdCD'])
    8
    >>> T(['c'*6, 'AbabAB', 'AbabA'+'C'*5+'bAbaCCCCBABccc']) # this and following examples from Hill-Pride-Vella
    4
    >>> T(['caCBCacb'])
    4
    >>> T(['c'*6,'AbabAB','abaBA'+'C'*5+'ABabACCCCAbaccc'])
    4
    >>> T(['ababcac','ababcACBABABabacaCBAB'])
    4
    >>> T(['a'*10,'c'*6,'aaCCCbcccccAbccAAAAAAAB'])
    4
    >>> T(['ac'*4,'acacACACACAbacaCAB'])
    4
    """
    if precomputed_link_graph is None:
        if noparse:
            rels=relator_list
        else:
            rels=parseinputwords(relator_list)
        if precomputed_piecedict is None:
            thepiecedict=piecedict(rels,piece_up_to_automorphism=piece_up_to_automorphism)
        else:
            thepiecedict=precomputed_piecedict
        thepiecesegments=piecesegments(rels,precomputed_piecedict=thepiecedict)
        the_link_graph=corner_graph(rels,unit_piecesegments(thepiecesegments),piece_up_to_automorphism=piece_up_to_automorphism)
    else:
        the_link_graph=precomputed_link_graph
    cyclelength= shortest_cycle_length(the_link_graph,immersed=True)
    return cyclelength

def Cprimebound(relator_list, precomputed_piecedict=None, noparse=False,piece_up_to_automorphism=True):
    return Cprime_bound(relator_list, precomputed_piecedict=None, noparse=False,piece_up_to_automorphism=True)
def Cprime_bound(relator_list, precomputed_piecedict=None, noparse=False,piece_up_to_automorphism=True):
    """
    The largest ratio of piece length to length of relator containing it.

    The group is C'(1/Lambda) for all Lambda such that 1/Lambda > Cprime_bound.

    Stop and return 1 if we find any such ratio >= 1/Lambda.

    >>> Cprime_bound(['abABcdCD'])
    Fraction(1, 8)
    >>> Cprime_bound(['aCbcABCac'])
    Fraction(2, 9)
    >>> Cprime_bound(['abcabcabc','dededede']) # no pieces
    0
    >>> Cprime_bound([[-2, -2, -3, -1, -2, -3, -3, -2, 1, 2, 2, 3, 3]])
    Fraction(3, 13)
    """
    if noparse:
        rels=relator_list
    else:
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels,piece_up_to_automorphism=piece_up_to_automorphism)
    else:
        thepiecedict=precomputed_piecedict
    thepiecesegments=piecesegments(rels,precomputed_piecedict=thepiecedict)
    return max((Fraction(l,len(rels[r])) for (r,v,e,l) in  thepiecesegments),default=0)

def C(relator_list,quit_at=float('inf'),piece_up_to_automorphism=True,precomputed_piecedict=None,noparse=False):
    """
    Find the minimum number p such that there exists some cyclic permutation of some relator that can be expressed as a freely reduced product of p pieces.

    If quit_at=q is specified the algorithm will stop and return q once it is determined that p>=q.

    If piece_up_to_automorphism=True then a word only counts as a pieces if it occurs in places in the relator list that are distinct up to automorphism. This means that copies of the root word in a relation that is a proper power do not yield pieces. 

    >>> C(['axB','bxC','cxD','dxE','exF','fxG','gxA'])
    3
    >>> C(['abABcdCD'])
    8
    >>> C([[1,2,-1,-2]])
    4
    >>> C(['abababc']) # 'c' is not a piece, so no relator can be written as a product of pieces.
    inf
    >>> C([[-2, -2, -3, -1, -2, -3, -3, -2, 1, 2, 2, 3, 3]])
    8
    """
    if noparse:
        rels=relator_list
    else:
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels,piece_up_to_automorphism=piece_up_to_automorphism)
    else:
        thepiecedict=precomputed_piecedict
    minrelatorpiecelength=float('inf')
    for relator_index in range(len(rels)):
        thisrelatorpiecelength=relator_piece_length(rels,relator_index,thepiecedict)
        minrelatorpiecelength=min(minrelatorpiecelength,thisrelatorpiecelength)
    return minrelatorpiecelength

def segment_piece_graph(rels, thesegment, thepiecedict):
    (relator_index,startvertex,direction,segmentlength)=thesegment
    relatorlength=len(rels[relator_index])
    endvertex=(thesegment[1]+direction*segmentlength)
    if direction==1 and endvertex>=relatorlength:
        wrap=True
    elif direction==-1 and endvertex<0:
        wrap=True
    else:
        wrap=False
    G=nx.DiGraph()
    G.add_nodes_from(n for n in range(startvertex,endvertex+direction,direction))
    thisrelatorsegments=(piece for piece in piecesegments(rels,precomputed_piecedict=thepiecedict) if piece[0]==relator_index and piece[2]==direction)
    this_segment_pieces=set()
    for piece in thisrelatorsegments:
        piecestart=piece[1]
        piecelength=piece[3]
        pieceend=piecestart+direction*piecelength
        if direction==1 and pieceend>=relatorlength:
            piecewrap=True
        elif direction==-1 and pieceend<0:
            piecewrap=True
        else:
            piecewrap=False
        if direction==1:
            if wrap == piecewrap:
                if startvertex<=piecestart and pieceend<=endvertex:
                    G.add_edge(piecestart,pieceend)
            elif wrap and not piecewrap:
                if startvertex<=piecestart:
                    G.add_edge(piecestart,pieceend)
                elif pieceend<=endvertex-relatorlength:
                    G.add_edge(piecestart+relatorlength,pieceend+relatorlength)
        else: #direction == -1
            if wrap == piecewrap:
                if piecestart<=startvertex and endvertex<=pieceend:
                    G.add_edge(piecestart,pieceend)
            elif wrap and not piecewrap:
                if piecestart<=startvertex:
                    G.add_edge(piecestart,pieceend)
                elif pieceend>=endvertex+relatorlength:
                    G.add_edge(piecestart-relatorlength,pieceend-relatorlength)
    return G

def segment_piece_length(rels,thesegment,thepiecedict):
    if thesegment[3]==0:
        return 0
    (relator_index,startvertex,direction,segmentlength)=thesegment
    relatorlength=len(rels[relator_index])
    endvertex=(thesegment[1]+direction*segmentlength)
    G=segment_piece_graph(rels, thesegment, thepiecedict)
    try:
        p=nx.shortest_path_length(G,startvertex,endvertex)
    except nx.NetworkXNoPath:
        p=float('inf')
    return p

def relator_piece_length(rels,relator_index,thepiecedict):
    relator_length=len(rels[relator_index])
    bestpiecelength=float('inf')
    for startvertex in range(relator_length):
        thispermutationpiecelength=segment_piece_length(rels,(relator_index,startvertex,1,relator_length),thepiecedict)
        bestpiecelength=min(bestpiecelength,thispermutationpiecelength)
    return bestpiecelength
            
def corner_remainder_segment(rels,corner,thepiecedict):
    """
    Given a corner, return the segment that is the remainder of the relator containing the corner after removing the two legs of the corner.
    """
    relator_index=corner[0][0]
    in_leg_1=reverse_segment(rels,corner[0])
    out_leg_2=corner[1]
    direction=out_leg_2[2]
    startvertex=(out_leg_2[1]+direction*out_leg_2[3])%len(rels[relator_index])
    length=len(rels[relator_index])-in_leg_1[3]-out_leg_2[3]
    return (relator_index,startvertex,direction,length)

def corner_remainder_piece_length(rels,corner,thepiecedict):
    remainder_segment=corner_remainder_segment(rels,corner,thepiecedict)
    return segment_piece_length(rels,remainder_segment,thepiecedict)

def subword(rels,segment):
    """
    Given the relatorlist and segment, return the subword of the relator defined by the segment.
    """
    (r,v,e,l)=segment
    w=rels[r]
    if e==1:
        return tuple((w+w)[v:v+l])
    if e==-1:
        W=inverse(w)
        return tuple((W+W)[-v%len(w):(-v%len(w))+l])

def reverse_segment(rels,segment):
    """
    The segment that defines the reverse traversal of the same path defined by the input segment. 
    """
    (r,v,e,l)=segment
    endvertex=(v+e*l)%len(rels[r])
    return (r,endvertex,-e,l)

def segments(rels):
    """
    Generator of all segments.
    """
    for r in range(len(rels)):
        w=rels[r]
        for v in range(len(w)):
            for e in [-1,1]:
                for l in range(1,len(w)+1):
                    yield (r,v,e,l)

def piecedict(rels,piece_up_to_automorphism=True):
    """
    Return a dictionary whose keys are pieces and whose values are sets of segments that specify that piece as their subword.
    """
    if piece_up_to_automorphism:
        roots,powers=zip(*[maxroot(relator) for relator in rels])
    thepiecesandwheretheycomefrom=dict()
    for segment in segments(rels):
        w=subword(rels,segment)
        thepiecesandwheretheycomefrom.setdefault(w,set([])).add(segment)
    if piece_up_to_automorphism:
        nonpieces=set()
        for k in thepiecesandwheretheycomefrom:
            origin_words={segment[0] for segment in thepiecesandwheretheycomefrom[k]}
            if len(origin_words)==1:
                origin_word_index=origin_words.pop()
                if len(thepiecesandwheretheycomefrom[k])==powers[origin_word_index]:
                    nonpieces.add(k)
    else:
        nonpieces={k for k in thepiecesandwheretheycomefrom if len(thepiecesandwheretheycomefrom[k])<2}
    for k in nonpieces:
        del thepiecesandwheretheycomefrom[k]
    return thepiecesandwheretheycomefrom
        
def piecesegments(rels,precomputed_piecedict=None):
    """
    Returns the set of all segments that correspond to pieces. 
    """
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    return set().union(*[thepiecedict[k] for k in thepiecedict])

def maximal_piecesegments(piecesegments):
    mpc=dict()
    for segment in piecesegments:
        mpc.setdefault(segment[:3],set([])).add(segment[3])
    return {k+(max(mpc[k]),) for k in mpc}

def unit_piecesegments(piecesegments):
    return {segment[:3]+(1,) for segment in piecesegments}
    
def successor_pieces(rels,thepiecesegments,thesegment):
    """
    Given a segment, yield all piecesegments that begin where the given segment ends. 
    """
    (r,v,e,l)=thesegment
    endvertex=(v+e*l)%len(rels[r])
    for nextl in range(1,len(rels[r])-l+1):
        nextsegment=(r,endvertex,e,nextl)
        if nextsegment in thepiecesegments:
            yield nextsegment

def interior_corners(rels,thepiecesegments):
    """
    Yield interior corners. They are pairs of piecesegments based at the same vertex of a relator and pointing away from the common vertex. 
    """
    for firstsegment in thepiecesegments:
        for secondsegment in successor_pieces(rels,thepiecesegments,firstsegment):
            yield (reverse_segment(rels,firstsegment),secondsegment)

def normalizedcornerangle(rels,corner,measure,precomputed_piecedict=None):
    """
    Give the angle of the corner divided by 2pi.
    """
    if measure=='metric':
        return Fraction(1,2)-Fraction(corner[0][3]+corner[1][3],2*len(rels[corner[0][0]]))
    elif measure=='nonmetric':
        if precomputed_piecedict==None:
            thepiecedict=piecedict(rels)
        else:
            thepiecedict=precomputed_piecedict
        crpl=corner_remainder_piece_length(rels,corner,thepiecedict)
        assert(crpl<float('inf'))
        return Fraction(1,2)-Fraction(1,2+crpl)
    else:
        raise ValueError("Argument 'measure' should be either 'metric' or 'nonmetric'.")

def corner_graph(rels,somepiecesegments,piece_up_to_automorphism=True, measure='metric',precomputed_piecedict=None):
    """
    Returns a weighted graph whose vertices are corners such that corner1 and corner2 are connected by a directed edge if the second leg of corner1 and the first leg of corner2 are distinct segments that define the same piece. A directed cycle in this graph corresponds to the arc neighborhood of an interior vertex in a reduced van Kampen diagram. 
    Vertices have attribte weight that is corner angle divided by 2pi; edges have attribute wieght that is average of its two vertex weights. Thus, the weight of a directed cycle represents the total angle around an interior vertex in a reduced van Kampen diagram, which is a positive curvature vertex if weight<1, flat if weight=1, and hyperbolic if weight>1.
    """
    if precomputed_piecedict==None:
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    if piece_up_to_automorphism:
        roots,powers=zip(*[maxroot(relator) for relator in rels])
    G=nx.DiGraph()
    for corner in interior_corners(rels,somepiecesegments):
        G.add_node(corner, weight=normalizedcornerangle(rels,corner,measure=measure,precomputed_piecedict=thepiecedict))
    for u in G:
        for v in G:
            if v[0]==u[1]: # second leg of corner u is same segment as first leg of corner v
                continue
            elif piece_up_to_automorphism and v[0][0]==u[1][0] and v[0][1]%len(roots[v[0][0]])==u[1][1]%len(roots[u[1][0]]): # second leg of corner u differs from  first leg of corner v by rotation of the relator
                continue
            elif subword(rels,v[0])!=subword(rels,u[1]): # second leg of corner u and first leg of corner v do not define same subword
                continue
            else:
                G.add_edge(u,v,weight=(G.nodes[u].get("weight")+G.nodes[v].get("weight"))/2)
    return G

def simple_light_loops_at(G,startvertex,maxweight=1):
    """
    Generator that yields simple loops based at given startvertex whose weight does not exceed maxweight. 
    """
    currentpath=[startvertex,]
    def advance(path):
        lastvertex=path[-1]
        potentialnext=sorted(G[lastvertex])
        for i in range(len(potentialnext)):
            if potentialnext[i] in path[1:]:
                continue
            elif potentialnext[i]==path[0]:
                return path+[potentialnext[i],]
            elif G.nodes[potentialnext[i]].get('weight')+sum(G.nodes[v].get('weight') for v in path)<=maxweight:
                return path+[potentialnext[i],]
        #backtrack
        for indextobacktrackto in range(len(path)-2,0,-1):
            prevvertex=path[indextobacktrackto]
            potentialnext=sorted(G[prevvertex])
            currentindex=potentialnext.index(path[indextobacktrackto+1])
            for newbranchindex in range(currentindex+1,len(potentialnext)):
                if potentialnext[newbranchindex] in path[1:]:
                    continue
                elif potentialnext[newbranchindex]==path[0]:
                    path[:indextobacktrackto+1]+[potentialnext[newbranchindex],]
                elif G.nodes[potentialnext[newbranchindex]].get('weight')+sum(G.nodes[v].get('weight') for v in path[:1+indextobacktrackto])<=maxweight:
                    return path[:indextobacktrackto+1]+[potentialnext[newbranchindex],]
        return []
    while currentpath:
        if len(currentpath)==1:
            currentpath=advance(currentpath)
        else:
            if currentpath[0]==currentpath[-1]:
                yield currentpath
            currentpath=advance(currentpath)
            
def simple_light_loops(G,maxweight=1):
    """
    Generator that yields simple loops in graph G whose weight does not exceed maxweight.
    Simple loop is given as a list of vertices such that first and last are the same, and no other repetitions, such that successive vertices have a directed edge between them in the graph. 
    """
    workingG=G.copy()
    nodes=sorted(G, key=lambda v: (G.nodes[v].get('weight'),v))
    for i in range(len(nodes)):
        for path in simple_light_loops_at(G.subgraph(nodes[i:]),nodes[i],maxweight):
            if len(path)>3:
                yield path
        workingG.remove_node(nodes[i]) # in later steps only consider paths that do not go through node[i]
    
def worst_vertex_curvature(relator_list,measure,sharing_factor=2,precomputed_piecedict=None,precomputed_corner_graph=None, relative_to_incoming_arc=None,relative_to_incoming_corners=None,noparse=False,verbose=False):
    if noparse:
        rels=relator_list
    else:
        if verbose:
            print("Parsing input words.")
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is not None:
        if verbose:
            print("Computing pieces.")
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    thepiecesegments=piecesegments(rels,thepiecedict)
    if precomputed_corner_graph is None:
        if verbose:
            print("Constructing  corner graph.")
        G=corner_graph(rels,thepiecesegments,measure=measure)
    else:
        G=precomputed_corner_graph
    if relative_to_incoming_arc is None and relative_to_incoming_corners is None:
        if verbose:
            print("Searching for short loops.")
        ssloops=simple_light_loops(G,maxweight=1)
        shortestlength=float('inf')
        shortestloop=None
        for loop in ssloops:
            looplength=sum(G.nodes[v].get('weight') for v in loop[1:])
            if looplength<shortestlength:
                shortestlength=looplength
                shortestloop=loop
        if verbose:
            if shortestloop is None:
                print("All vertices have negative  curvature.")
            else:
                print(shortestlength,shortestloop,[subword(rels,c[0]) for c in shortestloop[:-1]])
        return shortestlength
    elif relative_to_incoming_arc is not None and measure=='metric': # should require C(4) so that the different neighbors cannot interact with one another
        incomingsegement1=relative_to_incoming_arc[1]
        incomingsegement2=relative_to_incoming_arc[0]
        possible_starting_corners={v for v in G if v[0]==reverse_segment(rels,incomingsegement1)}
        possible_ending_corners={v for v in G if v[1]==reverse_segment(rels,incomingsegement2)}
        thisneightborminimalcontribution=min(Fraction(sharing_factor-1,sharing_factor*len(path))*(-1+G[v][u]['weight']+sum(G[path[j]][path[j+1]]['weight'] for j in range(-1+len(path)))) for u,v in itertools.product(possible_starting_corners,possible_ending_corners) for path in nx.all_shortest_paths(G,u,v,weight='weight') )
        return max(0,thisneightborminimalcontribution)
    elif relative_to_incoming_corners is not None and measure=='nonmetric':
        u,v=relative_to_incoming_corners
        thisneightborminimalcontribution=min(Fraction(1,sharing_factor*len(path)) * (-1+G[u][v]['weight']+sum(G[path[j]][path[j+1]]['weight'] for j in range(-1+len(path)))) for path in nx.all_shortest_paths(G,u,v,weight='weight'))
        return max(0,thisneightborminimalcontribution)
    else:
        assert(0)##########
                
def negative_metric_curvature(relator_list,precomputed_piecedict=None,noparse=False):
    return worst_vertex_curvature(relator_list,measure='metric',precomputed_piecedict=precomputed_piecedict,noparse=noparse) > 1

def negative_mean_metric_curvature(relator_list, sharing_factor=2, noparse=False,precomputed_piecedict=None,preC=None,verbose=False):
    assert(sharing_factor>1)
    if noparse:
        rels=relator_list
    else:
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        if verbose:
            print("Computing pieces.")
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    if preC is None:
        if verbose:
            print("Checking presentation is C(4).")
        theC=C(rels,noparse=True,precomputed_piecedict=thepiecedict)
    else:
        theC=preC
    assert(theC>=4) # Algorithm only implmented for C(4) presentaitons.
    thepiecesegments=piecesegments(rels,thepiecedict)
    if verbose:
        print("Constructing full corner graph.")
    G=corner_graph(rels,thepiecesegments,measure='metric')
    if verbose:
        print("Searching for short loops.")
    ssloops=simple_light_loops(G,maxweight=1)
    shortestlength=float('inf')
    shortestloop=None
    for loop in ssloops: # for light loops, equiv, for interior vertices with total angle <= 2pi, check if contributions from neighbor vertices will average out to give total angle >2pi.
        looplength=sum(G.nodes[v].get('weight') for v in loop[1:])
        if verbose:
            print("Found loop with value "+str(looplength)+", checking for neighbor support.")
        neighborcontributions=0
        for i in range(len(loop)-1):
            neighborcontributions+=worst_vertex_curvature(rels,measure='metric',sharing_factor=sharing_factor,precomputed_piecedict=thepiecedict,precomputed_corner_graph=G,relative_to_incoming_arc=(loop[i][1],loop[i+1][0]),noparse=True,verbose=False)
        if verbose:
            print("Neighbors contribute "+str(neighborcontributions)+".")
        if looplength+neighborcontributions<shortestlength:
            shortestlength= looplength+neighborcontributions
            shortestloop=loop
        if shortestlength<=1:
            if verbose:
                print(shortestlength,shortestloop,[subword(rels,c[0]) for c in shortestloop[:-1]])
            return False
    if verbose:
        if shortestloop is None:
            print("All interior vertices have negative curvature.")
        else:
            print(shortestlength,shortestloop,[subword(rels,c[0]) for c in shortestloop[:-1]])
    return True

def negative_nonmetric_curvature(relator_list, precomputed_piecedict=None,noparse=False):
    return worst_vertex_curvature(relator_list,measure='nonmetric',precomputed_piecedict=precomputed_piecedict,noparse=noparse) > 1

def negative_mean_nonmetric_curvature(relator_list, sharing_factor=2,noparse=False,precomputed_piecedict=None,preC=None,verbose=False):
    if noparse:
        rels=relator_list
    else:
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        if verbose:
            print("Computing pieces.")
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    if preC is None:
        if verbose:
            print("Checking presentation is C(4).")
        theC=C(rels,noparse=True,precomputed_piecedict=thepiecedict)
    else:
        theC=preC
    assert(theC>=4) # Algorithm only implmented for C(4) presentaitons.
    thepiecesegments=piecesegments(rels,thepiecedict)
    if verbose:
        print("Constructing full corner graph.")
    G=corner_graph(rels,thepiecesegments,measure='nonmetric')
    if verbose:
        print("Searching for short links.")
    ssloops=simple_light_loops(G,maxweight=1)
    shortestlength=float('inf')
    diagram_vertex_2_neighborhood=None
    def four_arc_corners(corner):
        reversedsecondsegment=corner[0]
        thirdsegment=corner[1]
        direction=thirdsegment[2]
        relator_index=thirdsegment[0]
        for reversedfirstsegment in successor_pieces(rels,thepiecesegments,reversedsecondsegment):
            firstsegment=reverse_segment(rels,reversedfirstsegment)
            secondsegment=reverse_segment(rels,reversedsecondsegment)
            for fourthsegment in successor_pieces(rels,thepiecesegments,thirdsegment):
                thesefoursegments=(firstsegment, secondsegment, thirdsegment,fourthsegment)
                assert(all(seg[2]==direction for seg in thesefoursegments))
                assert(all(seg[0]==relator_index for seg in thesefoursegments))
                yield thesefoursegments
    for loop in ssloops: # for light loops, equiv, for interior vertices with total angle <= 2pi, check if contributions from neighbor vertices will average out to give total angle >2pi.
        if verbose:
            print("Found small link, checking for neighbor support.")
        for fac in itertools.product(*[four_arc_corners(corner) for corner in loop[1:]]):
            new_central_vertex_angle=0
            neighborcontributions=0
            for x in fac:
                (firstsegment,secondsegment,thirdsegment,fourthsegment)=x
                direction=firstsegment[2]
                relator_index=firstsegment[0]
                totallength=firstsegment[3]+secondsegment[3]+thirdsegment[3]+fourthsegment[3]
                endvertex=fourthsegment[1]+direction*fourthsegment[3]
                cornerremaindersegment=(relator_index,endvertex,direction,len(rels[relator_index])-totallength)
                remainderarclength=segment_piece_length(rels,cornerremaindersegment,thepiecedict)
                thisnormalizedcornerangle=Fraction(1,2)-Fraction(1,4+remainderarclength)
                new_central_vertex_angle+=thisnormalizedcornerangle
            for i in range(len(fac)):
                thesefour=fac[i]
                nextfour=fac[(i+1)%len(fac)]
                neighborcorners=((reverse_segment(rels,nextfour[0]),nextfour[1]),(reverse_segment(rels,thesefour[2]),thesefour[3]))
                neighborworstcontribution=worst_vertex_curvature(rels,measure='nonmetric',sharing_factor=sharing_factor,precomputed_piecedict=thepiecedict,precomputed_corner_graph=G,relative_to_incoming_corners=neighborcorners,noparse=True)
                neighborcontributions+=neighborworstcontribution
            this_diagram_averaged_central_vertex_angle=new_central_vertex_angle+neighborcontributions
            if this_diagram_averaged_central_vertex_angle<shortestlength:
                shortestlength=this_diagram_averaged_central_vertex_angle
                diagram_vertex_2_neighborhood=fac
                if shortestlength<=1:
                    if verbose:
                        print(this_diagram_averaged_central_vertex_angle,diagram_vertex_2_neighborhood)
                    return False
    if verbose:
        if diagram_vertex_2_neighborhood is None:
            print("All interior vertices have negative curvature.")
        else:
            print(shortestlength,diagram_vertex_2_neighborhood)
    return True
            
                
            

def piece_length(theword,thepieces,quit_at=float('inf')):
    """
    Calculate minimal length of the given input string theword as a concatenation of thepieces. 

    Returns float('inf') if the input string cannot be written as a concatentation of thepieces. 
    """
    shortest_expression=shortest_piece_expression(theword,thepieces,quit_at)
    return len(shortest_expression) if shortest_expression is not None else float('inf')


def shortest_piece_expression(theword,thepieces,quit_at=float('inf'),as_cyclic_word=True):
    """
    Recursive determination of shortest expression of string theword as concatentation of thepieces.
  
    Return None if no such expression exists.

    Returned expression is guaranteed to be shortest possible, but is not necessarily unique piece expression of this length. 

    Default as_cyclic_word=True then theword is treated as cyclic word, so allow the possibility that a minimal expression is actually an expression of a cyclic permutation of theword. 
    """
    shortest_expression=None
    if as_cyclic_word:
        the_root,the_power=maxroot(theword)
        rotations=len(the_root)
    else:
        rotations=1
    currentword=theword
    currentbest=quit_at
    while rotations:
        for x in shorter_string_piece_expressions(currentword,thepieces,quit_at=currentbest):
             shortest_expression=x
        currentbest=len(shortest_expression) if shortest_expression is not None else currentbest
        currentword=currentword[1:]+currentword[0:1]
        rotations-=1
    return shortest_expression

def shorter_string_piece_expressions(theword,thepieces,quit_at=float('inf')):
    """
    Yield lists of decreasing length bounded above by quit_at whose elements are in thepieces and whose concatenation is theword.
    """
    if not theword:
        yield list([])
        return
    shortest_so_far=None
    if quit_at>0:
        for p in (p for p in thepieces if p==theword[:len(p)]):
            if shortest_so_far is None:
                completion_generator=shorter_string_piece_expressions(theword[len(p):],thepieces,quit_at-1)
            else:
                completion_generator=shorter_string_piece_expressions(theword[len(p):],thepieces,shortest_so_far-2)
            new_upper_bound=None
            while True:
                try:
                    therest=completion_generator.send(new_upper_bound)
                except StopIteration:
                    break
                if new_upper_bound is None: # first time through or stepping with next
                    new_upper_bound=yield [p,]+therest
                    shortest_so_far=min(1+len(therest),float('inf') if new_upper_bound is None else new_upper_bound,float('inf') if shortest_so_far is None else shortest_so_far)
                    new_upper_bound=min(len(therest),float('inf') if new_upper_bound is None else new_upper_bound-1,float('inf') if shortest_so_far is None else shortest_so_far-1)
                else: 
                    if len(therest)<new_upper_bound:
                        new_upper_bound=yield [p,]+therest
                        shortest_so_far=min(1+len(therest),float('inf') if new_upper_bound is None else new_upper_bound,float('inf') if shortest_so_far is None else shortest_so_far)
                        new_upper_bound=min(len(therest),float('inf') if new_upper_bound is None else new_upper_bound-1,float('inf') if shortest_so_far is None else shortest_so_far-1)
               
def all_piece_expressions(theword,thepieces, as_cyclic_word=True, quit_at=float('inf')):
    """
    Recursively yield lists of words of length at most quit_at whose elements are in thepieces and whose concatenation is either theword, or, when default as_cyclic_word=True, a cyclic permutation of theword. 
    When as_cyclic_word=True the yielded lists are normalized so that the first element of the list has a (possibly non-proper) suffix that agrees with a prefix of relator. 

    >>> [ expr for expr in all_piece_expressions('aba',['aa','b'])]
    [['aa', 'b']]
    >>> sorted([ expr for expr in all_piece_expressions('aaba',['a','aa','b'])])
    [['a', 'a', 'b', 'a'], ['aa', 'a', 'b'], ['aa', 'b', 'a']]
    >>> sorted([ expr for expr in all_piece_expressions('aaba',['a','aa','b'],quit_at=3)]) # only expressions of length <=3
    [['aa', 'a', 'b'], ['aa', 'b', 'a']]
    """
    r=theword
    if not r or quit_at==0:
        return list([])
    if as_cyclic_word:
        for p in (p for p in thepieces if len(p)<=len(r)):
            possiblestartingindices=[] # for given p there may be different possible choices of y
            for startingindex in range(len(r)-len(p)+1,len(r)+1):
                if p==(r+r)[startingindex:startingindex+len(p)]:
                    possiblestartingindices.append(startingindex)
            if not possiblestartingindices:
                continue
            for startingindex in possiblestartingindices:
                # found a way to fit p into r spanning the beginning of r. 
                whatwevegot=[p,]
                whatsleft=(r+r)[startingindex+len(p):startingindex+len(r)]
                if whatsleft:
                     for therest in all_piece_expressions(whatsleft,thepieces,quit_at=quit_at-1,as_cyclic_word=False):
                        yield whatwevegot+therest
                else:
                    yield whatwevegot
    else:
        for p in (p for p in thepieces if r[:len(p)]==p):
            whatwevegot=[p,]
            whatsleft=r[len(p):]
            if whatsleft:
                for therest in all_piece_expressions(whatsleft,thepieces,quit_at=quit_at-1,as_cyclic_word=False):
                    yield whatwevegot+therest
            else:
                yield whatwevegot
        
        
def pieces(relatorlist,piece_up_to_automorphism=True,noparse=False,asstring=False):
    """
    Given input container of relators, return set of pieces, which are subwords occuring more than once in relators or their inverses, as cyclic words.

    If piece_up_to_automorphism=True then do not count as a piece a subword of a relator that occurs only periodically.

    >>> pieces(['aa'],asstring=True,piece_up_to_automorphism=True)
    set()
    >>> sorted(pieces(['aa'],piece_up_to_automorphism=False,asstring=True))
    ['A', 'AA', 'a', 'aa']
    >>> sorted(pieces(['abABcdCD'],asstring=True))
    ['A', 'B', 'C', 'D', 'a', 'b', 'c', 'd']
    >>> sorted(pieces(['abbccbABCCB'],asstring=True))
    ['A', 'B', 'BC', 'BCC', 'BCCB', 'C', 'CB', 'CC', 'CCB', 'a', 'b', 'bc', 'bcc', 'bccb', 'c', 'cb', 'cc', 'ccb']
    >>> sorted(pieces([[1,2,3,1,2],[2,2,1,1,2,3]]))
    [(-3,), (-3, -2), (-3, -2, -1), (-2,), (-2, -1), (-1,), (-1, -2), (1,), (1, 2), (1, 2, 3), (2,), (2, 1), (2, 3), (3,)]
    """
    if noparse:
        rels=relatorlist
    else:
        rels=parseinputwords(relatorlist,asrelatorlist=True)
    thepieces=piecedict(rels,piece_up_to_automorphism=piece_up_to_automorphism).keys()
    if asstring:
        return {intlisttostring(piece) for piece in thepieces}
    else:
        return thepieces


 
            

def shortest_cycle_length(inputgraph,at_vertex=None, first_edge=None, immersed=False,weight=None):
    """
    Return the length of the shortest (directed) cycle in the graph, or float('inf') if none exists. 
    If weight=string (typically string='weight')  is given then use edge attribte string as edge lengths. If weigh=None then all edges have length 1.

    If one of at_vertex = v or first_edge=e then give the length of the shortest cycle that either starts at v, or has first edge e, respectively. 

    If immersed=True do not allow cycles that have edge (u,v) followed by edge (v,u). If, in addition, at_vertex = v or first_edge=e, then it is possible that the shortest cycle is not simple. Eg, if e is the bar of the barbell graph then the shortest immersed cycle starting with e is of the form e + loop + e backwards + loop. 
    If immersed=1 and at_vertex = v or first_edge=e then reuqire loops to be immersed except possibly at the initial vertex.
    """
    G=inputgraph.copy()
    if first_edge is not None:
        theedges=[first_edge,]
    elif at_vertex is None:
        theedges=[e for e in G.edges()]
    else:
        theedges=[e for e in G.edges(at_vertex)]
    shortestcyclelength=float('inf')
    for e in theedges:
        if weight is None:
            first_edge_weight=1
        else:
            first_edge_weight=G[e[0]][e[1]][weight]
        G.remove_edge(*e)
        if immersed:
            reverse_edge_is_present=False
            if (e[1],e[0]) in G.edges():
                reverse_edge_is_present=True
                if weight is None:
                    reverse_edge_weight=1
                else:
                    reverse_edge_weight=G[e[1]][e[0]][weight]
                G.remove_edge(e[1],e[0])
        try:
            shortest_e_path=nx.shortest_path_length(G,e[1],e[0],weight=weight) # compute distance between endpoints of e in G-e then add 1
        except nx.NetworkXNoPath:
            shortest_e_path=float('inf')
        shortest_simple_cycle_using_e=shortest_e_path+first_edge_weight
        shortestcycleusing_e=shortest_simple_cycle_using_e
        if immersed and reverse_edge_is_present and (first_edge is not None or at_vertex is not None): # in this case it is possible that there is no simple cycle using e as first edge, but there may be an immersed cycle.
            shortest_1_loop=shortest_cycle_length(G,at_vertex=e[1],immersed=1,weight=weight)
            if immersed is True:
                shortest_0_loop=shortest_cycle_length(G,at_vertex=e[0],immersed=1,weight=weight)
            else:
                shortest_0_loop=0
            shortestcycleusing_ebar=first_edge_weight+shortest_1_loop+reverse_edge_weight+shortest_0_loop
            shortestcycleusing_e=min(shortest_simple_cycle_using_e,shortestcycleusing_ebar)
        if reverse_edge_is_present:
            if weight is None:
                G.add_edge(e[1],e[0])
            else:
                G.add_edge(e[1],e[0],weight=reverse_edge_weight)
        shortestcyclelength=min(shortestcyclelength,shortestcycleusing_e)
    return shortestcyclelength


   
def common_prefix_length(stringone,stringtwo):
    """
    Given two strings, return the length of their longest common prefix. 
    """
    if len(stringone)==0 or len(stringtwo)==0:
        return 0
    for L in range(min(len(stringone),len(stringtwo))):
        if stringone[L]!=stringtwo[L]:
            return L
    return min(len(stringone),len(stringtwo))



def parseinputwords(inputwords,asrelatorlist=True):
    """
    Take as input an iterator of things representing words in a free group. Return representation of same words given as tuples of nonzero integers where positive i indicates the i-th generator of the free group and -i indicates its inverse. 

    If input consists of alphabetic strings then it is convertex to tuples of integers with a -> 1, A -> -1, b -> 2, B-> -2, etc.

    If asrelatorlist=True, raise an exception if some word in the list is not freely or cyclically reduced, or if some word in the list has length less than 2, or if some pair of words are conjugates or inverse conjugates.  
    """
    if all(type(w)==str and (w=='' or w.isalpha()) for w in inputwords):
        rels=tuple([tuple(stringtointlist(w)) for w in inputwords])
    elif all(hasattr(w,"letters") for w in inputwords):
        rels=tuple([tuple(w.letters) for w in inputwords])
    elif all(all(type(x)==int and x!=0 for x in w) for w in inputwords):
        rels=tuple([tuple(w) for w in inputwords])
    else:
        raise ValueError('Input must be either a list of lists of nonzero integers or a list of alphabetic strings.')
    if not asrelatorlist:
        return rels
    if any(len(w)<2 for w in rels):
        raise ValueError("Some input word has length less than 2.")
    if not all(freely_reduced(w) and cyclically_reduced(w) for w in rels):
        raise ValueError('Some input word not freely or cyclically reduced.')
    for i in range(len(rels)):
        for j in range(i+1,len(rels)):
            if are_conjugate(rels[i],rels[j]) or are_conjugate(rels[i],inverse(rels[j])):
                raise ValueError("Input words are not unique up to inversion and conjugation.")
    return rels
        
# Some utility function for words in free groups.
def are_conjugate(u,v):
    """
    Decide if two input lists of nonzero integers represent conjugate elements of the free group.

    >>> are_conjugate([],[1,2,-2,-1])
    True
    >>> are_conjugate([1,2,-1,-2],[-2,1,2,-1])
    True
    >>> are_conjugate([1,2,3],[3,2,1])
    False
    >>> are_conjugate([1,2,1,-1,2,1,-2,-1],[1,2])
    True
    >>> are_conjugate([1,2,-1],[2,1,-2])
    False
    """
    if Abelianization(u)!=Abelianization(v):
        return False
    x=cyclicreduce(freelyreduce(u))
    y=cyclicreduce(freelyreduce(v))
    if len(x)==0 and len(y)==0:
        return True
    elif len(x)==0 or len(y)==0:
        return False
    for i in range(len(y)):
        if x==y[i:]+y[:i]:
            return True
    return False

def Abelianization(intlist):
    """
    >>> Abelianization([])
    []
    >>> Abelianization([1,1,-1,-1])
    []
    >>> Abelianization([1,2,-1,2,2,2,2,-1,-2,-2,4])
    [-1, 2, 2, 2, 4]
    """
    rank=max((abs(x) for x in intlist),default=0)
    abelianized=[]
    vectorform=dict()
    for x in intlist:
        if x>0:
            vectorform[x]=vectorform.setdefault(x,0)+1
        if x<0:
            vectorform[-x]=vectorform.setdefault(-x,0)-1
    for i in range(1,rank+1):
        if vectorform.setdefault(i,0)>=0:
            abelianized+=vectorform[i]*[i,]
        else:
            abelianized+=(-vectorform[i])*[-i,]
    return abelianized
        

def maxroot(thestring):
    """
    Given an input string, return the shortest string of which the input string is a positive multiple, and the multiple.

    >>> maxroot('abcabcabc')
    ('abc', 3)
    >>> maxroot('abababababa')
    ('abababababa', 1)
    >>> maxroot('')
    ('', 1)
    >>> maxroot([1,2,1,2])
    ([1, 2], 2)
    """
    if len(thestring)<=1:
        return thestring,1
    for the_power in (n for n in range(len(thestring),0,-1) if len(thestring)%n==0):
        if thestring==the_power*thestring[:len(thestring)//the_power]:
            the_root=thestring[:len(thestring)//the_power]
            return the_root,the_power
    


def freely_reduced(intlist):
    """
    >>> freely_reduced([])
    True
    >>> freely_reduced([1,2,3])
    True
    >>> freely_reduced([1,2,1,-1,2])
    False
    """
    return all(intlist[i]!=-intlist[i+1] for i in range(len(intlist)-1))

def cyclically_reduced(intlist):
    """
    >>> cyclically_reduced([])
    True
    >>> cyclically_reduced([1,-1,2,-2])
    False
    >>> cyclically_reduced([1,2,-1,-2,1])
    True
    >>> cyclically_reduced([1,2,1,2,-1])
    False
    """
    if len(intlist)==0:
        return True
    else:
        return freely_reduced(intlist) and intlist[0]!=-intlist[-1]

def freelyreduce(intlist):
    """
    >>> freelyreduce([])
    ()
    >>> freelyreduce([1,2,3])
    (1, 2, 3)
    >>> freelyreduce([1,2,1,-1,2])
    (1, 2, 2)
    >>> freelyreduce([1,2,-2,3,1,2,-2,-1,-3,-1])
    ()
    """
    reduced=[x for x in intlist]
    if len(reduced)<2:
        return tuple(reduced)
    currentindex=0
    while currentindex<len(reduced)-1:
        if reduced[currentindex]==-reduced[currentindex+1]:
            try:
                reduced=reduced[:currentindex]+reduced[currentindex+2:]
                currentindex=max(0,currentindex-1)
            except IndexError:
                reduced=reduced[:currentindex]
        else:
            currentindex+=1
    return tuple(reduced)
        

def cyclicreduce(intlist):
    """
    >>> cyclicreduce([])
    ()
    >>> cyclicreduce([1,-1])
    ()
    >>> cyclicreduce([2,1,-1,2,3,-2,-2])
    (3,)
    >>> cyclicreduce([1,2,1,2,-1,-2,-2,-1])
    (1, 2, -1, -2)
    """
    theword=freelyreduce(intlist)
    if len(theword)<2:
        return tuple(theword)
    conjugatorindex=0
    while conjugatorindex<=len(theword)//2 and theword[conjugatorindex]==-theword[-conjugatorindex-1]:
        conjugatorindex+=1
    return tuple(theword[conjugatorindex:len(theword)-conjugatorindex])

def inverse(input_word):
    """
    Return the inverse of the given input_word in the same form as given.

    >>> inverse([])
    ()
    >>> inverse([1])
    (-1,)
    >>> inverse([-2])
    (2,)
    >>> inverse([1,2,-3])
    (3, -2, -1)
    >>> inverse('')
    ''
    >>> inverse('a')
    'A'
    >>> inverse('B')
    'b'
    >>> inverse('abC')
    'cBA'
    """
    if (type(input_word)==list or type(input_word)==tuple) and all(type(x)==int for x in input_word):
        return tuple(map(lambda x:-1*x, input_word[::-1]))
    elif type(input_word)==str and (input_word=='' or input_word.isalpha()):
        return input_word[::-1].swapcase()
    else:
        raise ValueError("Input should be either a list of non-zero integers or an alphabetic string.")

def simple_Whitehead_graph(rels):
    rank=max(max(abs(x) for x in w) for w in rels)
    G=nx.Graph()
    for i in range(1,rank+1):
        G.add_node(i)
        G.add_node(-i)
    for w in rels:
        if len(w)==1:
            G.add_edge(w[0],-w[0])
        if len(w)>1:
            for i in range(len(w)):
               G.add_edge(-w[i-1],w[i])
    return G


def stringtointlist(thestring):
    intlist=[]
    for c in thestring:
        if c.islower():
            intlist.append(1+'abcdefghijklmnopqrstuvwxyz'.index(c))
        else:
            intlist.append(-1*(1+'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.index(c)))
    return intlist

def intlisttostring(intlist):
    thestring=''
    if any(x==0 for x in intlist) or any(abs(x)>26 for x in intlist):
        raise ValueError("Input integer out of alphabet range.")
    for x in intlist:
        if x>0:
            thestring+='abcdefghijklmnopqrstuvwxyz'[x-1]
        else:
            thestring+='ABCDEFGHIJKLMNOPQRSTUVWXYZ'[-x-1]
    return thestring
        

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
            
    
