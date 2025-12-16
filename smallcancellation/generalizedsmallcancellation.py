import networkx as nx
import math
from fractions import Fraction
import itertools
import functools


# Input relators as either list of alphabetic strings or list of lists of nonzero integers. For alphabetic strings, lower case letters a-z represent generators 1-26 of a free group, and change of case denotes inversion, A=a^-1, B=b^-1 etc. For list of integers, positive integers represent generators and negative integers represent their inverses. For example, [[1,2,-1,-2]] and ['abAB'] both represent lists of relators with one element that is the commutator of the first two generators.

# relatorlist returned by parseinputwords as list of lists of nonzero integers

# segment=(r,v,e,l) means take the reltor at index r in relatorlist, subword starting at vertex v, direction e +1 or -1, and length l.

# Vertex v means the place between index v-1 and index v. If direction is +1 then the letter at vertex v is contained at index v. If direction is -1 then the letter at vertex v is the inverse of the letter at index v-1.

# the subword of a segment is the tuple of nonzero integers that is the subword of a relator corresponding to the given segment. Example: relator_list=[(1,2,3)] here are a few examples of segment -> subword: (0,0,1,2) -> (1,2), (0,2,-1,2)->(-2,-1), (0,2,1,2)->(3,1)

# a piece is a word that occurs as the subword corresponding to distinct segments

# a piece-segment is a segment whose subword is a piece

# a corner is a pair of piece-segments with the same relator, same vertex, opposite directions, and whose lengths sum to at most the length of the relator. 



def CT(relator_list,quit_at=float('inf'),piece_up_to_automorphism=True,precomputed_piecedict=None,noparse=False, precomputed_corner_graph=None):
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
    return C(rels,precomputed_piecedict=thepiecedict),T(rels,precomputed_piecedict=thepiecedict,precomputed_corner_graph=precomputed_corner_graph)

    
def T(relator_list, precomputed_piecedict=None, noparse=False, precomputed_corner_graph=None,piece_up_to_automorphism=True):
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
    if precomputed_corner_graph is None:
        if noparse:
            rels=relator_list
        else:
            rels=parseinputwords(relator_list)
        if precomputed_piecedict is None:
            thepiecedict=piecedict(rels,piece_up_to_automorphism=piece_up_to_automorphism)
        else:
            thepiecedict=precomputed_piecedict
        thepiecesegments=piecesegments(rels,precomputed_piecedict=thepiecedict)
        G=unweighted_corner_graph(rels,unit_piecesegments(thepiecesegments),piece_up_to_automorphism=piece_up_to_automorphism)
    else:
        G=precomputed_corner_graph
    return shortest_cycle_length(G,nobigon=True)

def Cprime_bound(relator_list, precomputed_piecedict=None, noparse=False,piece_up_to_automorphism=True):
    """
    The largest ratio of piece length to length of relator containing it.

    The group is C'(1/Lambda) for all Lambda such that 1/Lambda > Cprime_bound.


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

def C(relator_list,quit_at=float('inf'),piece_up_to_automorphism=True,precomputed_piecedict=None,noparse=False, verbose=False):
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
    relator_piece_decompositions=[relator_piece_decomposition(rels,relator_index,thepiecedict)  for relator_index in range(len(rels))]
    relator_piece_lengths=[]
    for D in relator_piece_decompositions:
        if D is None:
            relator_piece_lengths.append(float('inf'))
        else:
            relator_piece_lengths.append(len(D))
    min_decomp_length=min(relator_piece_lengths)
    if verbose:
        min_decomp_index=relator_piece_lengths.index(min_decomp_length)
        decompsegments=relator_piece_decompositions[min_decomp_index]
        outputstring=intlisttostring(rels[min_decomp_index])+'~'
        for i in range(len(decompsegments)-1):
            outputstring+=intlisttostring(subword(rels,decompsegments[i]))+'+'
        outputstring+=intlisttostring(subword(rels,decompsegments[-1]))
        return min_decomp_length,outputstring
    else:
        return min_decomp_length

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



def segment_piece_length(rels,segment,precomputed_piecedict=None, piece_segment_weights=None):
    """
    Among piece decompositions of the subword of the given segment, return the minimum of their weighted lengths. 
    """
    if segment[3]==0:
        return 0
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    if piece_segment_weights is None:
        thepsweights=lambda x:1
    else:
        thepsweights=piece_segment_weights
    (relator_index,startvertex,direction,segmentlength)=segment
    relatorlength=len(rels[relator_index])
    endvertex=(segment[1]+direction*segmentlength)
    G=segment_piece_graph(rels, segment, thepiecedict,thepsweights)
    try:
        p=nx.shortest_path_length(G,startvertex,endvertex,weight="weight")
    except nx.NetworkXNoPath:
        p=float('inf')
    return p

def segment_piece_decomposition(rels,segment,precomputed_piecedict=None,piece_segment_weights=None):
    """
    Return a shortest decomposition of the given segment as a concatenation of pieces.
    """
    if segment[3]==0:
        return []
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    if piece_segment_weights is None:
        thepsweights=lambda x:1
    else:
        thepsweights=piece_segment_weights
    (relator_index,startvertex,direction,segmentlength)=segment
    relatorlength=len(rels[relator_index])
    endvertex=(segment[1]+direction*segmentlength)
    G=segment_piece_graph(rels, segment, thepiecedict, thepsweights)
    try:
        p=nx.shortest_path(G,startvertex,endvertex)
    except nx.NetworkXNoPath:
        return None
    return [subword(rels,G[p[i]][p[i+1]]['label']) for i in range(len(p)-1)]


    
def segment_piece_graph(rels, thesegment, thepiecedict, piece_segment_weights=None):
    """
    Auxiliary function for segment_piece_length, segment_piece_decomposition. Constructs a digraph whose vertices are piece-segments and there is a directed edge denotes successor. If piece_segment_weights is given it is used to weight the edges. 
    """
    if piece_segment_weights is None:
        thepsweights=lambda x:1
    else:
        thepsweights=piece_segment_weights
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
    thisrelatorsegments=(ps for ps in piecesegments(rels,precomputed_piecedict=thepiecedict) if ps[0]==relator_index and ps[2]==direction)
    this_segment_pieces=set()
    for ps in thisrelatorsegments:
        piecestart=ps[1]
        piecelength=ps[3]
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
                    G.add_edge(piecestart,pieceend,label=ps,weight=thepsweights(ps))
            elif wrap and not piecewrap:
                if startvertex<=piecestart:
                    G.add_edge(piecestart,pieceend,label=ps,weight=thepsweights(ps))
                elif pieceend<=endvertex-relatorlength:
                    G.add_edge(piecestart+relatorlength,pieceend+relatorlength,label=ps,weight=thepsweights(ps))
        else: #direction == -1
            if wrap == piecewrap:
                if piecestart<=startvertex and endvertex<=pieceend:
                    G.add_edge(piecestart,pieceend,label=ps,weight=thepsweights(ps))
            elif wrap and not piecewrap:
                if piecestart<=startvertex:
                    G.add_edge(piecestart,pieceend,label=ps,weight=thepsweights(ps))
                elif pieceend>=endvertex+relatorlength:
                    G.add_edge(piecestart-relatorlength,pieceend-relatorlength,label=ps,weight=thepsweights(ps))
    return G

def relator_piece_length(rels,relator_index,thepiecedict):
    """
    Find a shortest expression of a cyclic permutation of the given relator as a concatentation of pieces. 
    """
    relator_length=len(rels[relator_index])
    bestpiecelength=float('inf')
    for startvertex in range(relator_length):
        thispermutationpiecelength=segment_piece_length(rels,(relator_index,startvertex,1,relator_length),thepiecedict)
        bestpiecelength=min(bestpiecelength,thispermutationpiecelength)
    return bestpiecelength

def relator_piece_decomposition(rels,relator_index,thepiecedict):
    """
    Return a shortest decomposition of (a cyclic conjugate of) the given relator as a concatentation of pieces.

    Output is list of piece-segments whose concatenation is a cyclic permutation of the relator. 
    """
    relator_length=len(rels[relator_index])
    piece_decomposition_by_starting_index=[]
    piece_decomposition_by_starting_index_lengths=[]
    for startvertex in range(relator_length):
        D=segment_piece_decomposition(rels,(relator_index,startvertex,1,relator_length),thepiecedict)
        piece_decomposition_by_starting_index.append(D)
        if D is None:
            piece_decomposition_by_starting_index_lengths.append(float('inf'))
        else:
            piece_decomposition_by_starting_index_lengths.append(len(D))
    min_piece_length=min(piece_decomposition_by_starting_index_lengths)
    min_length_index=piece_decomposition_by_starting_index_lengths.index(min_piece_length)
    return piece_decomposition_by_starting_index[min_length_index]
            
def corner_remainder_segment(rels,corner):
    """
    Given a corner, return the segment that is the remainder of the relator containing the corner after removing the two legs of the corner.
    """
    relator_index=corner[0][0]
    in_leg=reverse_segment(rels,corner[0])
    out_leg=corner[1]
    direction=out_leg[2]
    startvertex=(out_leg[1]+direction*out_leg[3])%len(rels[relator_index])
    length=len(rels[relator_index])-in_leg[3]-out_leg[3]
    return (relator_index,startvertex,direction,length)


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
    """
    Returns the set of all piece-segments that are not a subsegment of a longer piece-segment.
    """
    mpc=dict()
    for segment in piecesegments:
        mpc.setdefault(segment[:3],set([])).add(segment[3])
    return {k+(max(mpc[k]),) for k in mpc}

def unit_piecesegments(piecesegments):
    """
    Return the set of piece-segments of length 1.
    """
    return {segment[:3]+(1,) for segment in piecesegments}
    
def successor_pieces(rels,thepiecesegments,thesegment):
    """
    Given a segment, yield all piece-segments that begin where the given segment ends. 
    """
    (r,v,e,l)=thesegment
    endvertex=(v+e*l)%len(rels[r])
    for nextl in range(1,len(rels[r])-l+1):
        nextsegment=(r,endvertex,e,nextl)
        if nextsegment in thepiecesegments:
            yield nextsegment

def interior_corners(rels,thepiecesegments):
    """
    Yield interior corners. They are pairs of piece-segments based at the same vertex of a relator and pointing away from the common vertex, such that the sum of their lengths is at most the length of the relator. 
    """
    for firstsegment in thepiecesegments:
        relatorlength=len(rels[firstsegment[0]])
        for secondsegment in [secondsegment for secondsegment in successor_pieces(rels,thepiecesegments,firstsegment) if firstsegment[3]+secondsegment[3]<=relatorlength]:
            yield (reverse_segment(rels,firstsegment),secondsegment)

def corner_angle_from_generator_weights(rels,corner,generator_weights,precomputed_piecedict=None):
    """
    Return the interior angle of the corner expressed in turns.
    weights is tuple of integral weights of the group generators
    """
    leg0=weighted_word_length(subword(rels,corner[0]),generator_weights)
    leg1=weighted_word_length(subword(rels,corner[1]),generator_weights)
    remainder_segment=corner_remainder_segment(rels,corner)
    restlength=weighted_word_length(subword(rels,remainder_segment),generator_weights)
    return Fraction(1,2)*(1-Fraction(leg0+leg1,leg0+leg1+restlength))


def corner_angle_from_piece_weights_with_known_boundary(rels,corner,boundary_weight,piece_sgement_weights):
    return Fraction(piece_segment_weights(corner[0])+piece_segment_weights(corner[1]),2*boundary_weight)
    

def weighted_word_length(theword,theweights):
    """
    Return the weighted word length the in input work, where generator i has length theweights[i-1].
    """
    wwl=0
    for x in theword:
        wwl+=theweights[abs(x)-1]
    return wwl

def piece_segment_weights_from_piece_weights(rels,thepiecedict,thepieceweights):
    return {ps:thepieceweights[subword(rels,ps)] for ps in piecesegments(rels,thepiecedict)}

def unweighted_corner_graph(rels,somepiecesegments,piece_up_to_automorphism=True):
    """
    Returns a digraph whose vertices are corners such that corner1 and corner2 are connected by a directed edge if the second leg of corner1 and the first leg of corner2 are distinct segments that define the same piece. A directed cycle in this graph corresponds to the arc neighborhood of an interior vertex in a reduced van Kampen diagram. 
    """
    if piece_up_to_automorphism:
        roots,powers=zip(*[maxroot(relator) for relator in rels])
    G=nx.DiGraph()
    for corner in interior_corners(rels,somepiecesegments):
        G.add_node(corner)
    if piece_up_to_automorphism and any(x>1 for x in powers):
            for u in G:
                for v in [v for v in G if v[0]!=u[1] and subword(rels,v[0])==subword(rels,u[1])]:
                    if piece_up_to_automorphism and v[0][0]==u[1][0] and powers[v[0][0]]>1 and (v[0][1]%len(roots[v[0][0]]))==(u[1][1]%len(roots[u[1][0]])): # second leg of corner u differs from  first leg of corner v by rotation of the relator, no edge here, would make reduced diagram
                        continue
                    else:
                        G.add_edge(u,v)
    else:
        for u in G:
            for v in [v for v in G if v[0]!=u[1] and subword(rels,v[0])==subword(rels,u[1])]:
                G.add_edge(u,v)
    return G

def weight_corner_graph(rels,precomputed_corner_graph,corner_angle_function,precomputed_piecedict=None):
    """
    Given a corner graph and an angle_function that takes input a corner and outputs an angle, compute the interior angle of each corner and write it in node attribute 'weight'. Also write to each edge attribute 'weight' that is the average of its two node weights. 
    """
    G=precomputed_corner_graph
    for corner in G:
        G.nodes[corner]['weight']=corner_angle_function(corner)
    for (u,v) in G.edges:
        G.edges[u,v]['weight'] = Fraction(G.nodes[u].get('weight')+G.nodes[v].get('weight'),2)
    return 



def light_loops_with_initial_segment(G,initialsegment,maxweight=1,simple=False,upperdensitybound=float('inf'),upperlengthbound=float('inf')):
    """
    Given a digraph with weighted vertices and a nonempty list initialsegment of vertices corresponding to a directed segment in the graph, generate all possible loops begining with that initial segment such that the sum of vertex weights of the loop does not exceed maxweight.
    If 'simple' is True only yield simple loops. 
    """
    if not initialsegment:
        raise ValueError("Must give nonempty initial segment.")
    if maxweight==float('inf') and upperlengthbound==float('inf') and not simple:
        raise ValueError('Must give weight or length bound or restrict to simple loops.')
    initialweight=sum(G.nodes[v]['weight'] for v in initialsegment)
    if initialweight>maxweight or len(initialsegment)>upperlengthbound or Fraction(initialweight,len(initialsegment))>upperdensitybound: # given initial segment already violates one of the restrictions
        return
    if initialweight==maxweight or len(initialsegment)==upperlengthbound:
        if initialsegment[0] in G[initialsegment[-1]]:
            if not simple or  len(set(initialsegment))==len(initialsegment):
                return initialsegment
            else:
                return
        else:
            return
    if simple and len(set(initialsegment))<len(initialsegment):# can't make a simple path if given initialsegment is nonsimple
        return
    lastvertex=initialsegment[-1]
    potentialnext=sorted([v for v in G[lastvertex] if v==initialsegment[0] or initialweight+G.nodes[v]['weight']<maxweight or (initialweight+G.nodes[v]['weight']==maxweight and initialsegment[0] in G[v])],key=lambda v:-G.nodes[v]['weight'])
    currentpath=[[lastvertex,potentialnext]]
    maxsuffixlength=upperlengthbound-len(initialsegment)+1
    maxsuffixweight=maxweight-initialweight
    def advance(path):
        currentpathweight=sum(G.nodes[p[0]]['weight'] for p in path[1:])
        if (len(path)==maxsuffixlength or currentpathweight==maxsuffixweight) and initialsegment[0] in path[-1][1]:
            return path[:-1]+[[path[-1][0],[]]]+[[initialsegment[0],[]]]
        elif len(path)<maxsuffixlength and currentpathweight<maxsuffixweight and path[-1][1] and (not simple or path[-1][0]!=initialsegment[0]): # try to add another vertex to the path
            potentialnext=path[-1][1].pop()
            if potentialnext==initialsegment[0]:
                if simple or len(path)+1==maxsuffixlength or G.nodes[potentialnext]['weight']+currentpathweight>=maxsuffixweight: # can close to a loop here, but cannot extend any further because it would go out of bounds.
                    return path+[[potentialnext,[]]]
                else:
                    return path+[[potentialnext,sorted([v for v in G[potentialnext] if currentpathweight+G.nodes[potentialnext]['weight']+G.nodes[v]['weight']<maxsuffixweight or (currentpathweight+G.nodes[potentialnext]['weight']+G.nodes[v]['weight']==maxsuffixweight and initialsegment[0] in G[v])],key=lambda v:-G.nodes[v]['weight'])]]
            else:
                if simple and potentialnext in initialsegment+[p[0] for p in path]:
                    return path
                elif len(path)+1==maxsuffixlength or currentpathweight+G.nodes[potentialnext]['weight']==maxsuffixweight:
                    if initialsegment[0] in G[potentialnext]:
                        return path+[[potentialnext,[initialsegment[0]]]]
                    else:
                        return path+[[potentialnext,[]]]
                else:
                    return path+[[potentialnext,sorted([v for v in G[potentialnext] if v==initialsegment[0] or currentpathweight+G.nodes[potentialnext]['weight']+G.nodes[v]['weight']<=maxsuffixweight],key=lambda v:G.nodes[v]['weight'])]]
        else:
            return path[:-1]
             
    while currentpath:
        if len(currentpath)==1:
            currentpath=advance(currentpath)
        else:
            if initialsegment[0]==currentpath[-1][0] and currentpath[-1][1]==[]: # this is a closed loop, and we have already considered any allowed extensions
                yield initialsegment+[p[0] for p in currentpath[1:-1]]
            currentpath=advance(currentpath)
            
def light_loops(G,maxweight=1,simple=False):
    """
    Generator that yields unbased loops of length at least 3 in vertex weighted digraph G whose weight does not exceed maxweight.
    Loop is given as a list of vertices such that that successive vertices have a directed edge between them in the graph, and such that there is an edge in the graph from last vertex to first vertex.
    If 'simple' is True only yield loops in which no vertex is repeated.
    """
    if maxweight==float('inf') and not simple:
        raise ValueError('Must give finite maxweight or restrict to simple loops.')
    workingG=G.copy()
    # The corner graph we are actually interested in has a loop symmetry reversing orientation. Normalize with chosen orientation such that the first vertex of each loop is a corner whose directions are -1,1. To achieve this sort the nodes such that vertices with v[0][2]==-1 come first, and only yield loops that start with one of those vertices.  
    nodes=sorted(G, key=lambda v: (v[0][2],G.nodes[v]['weight'],v))
    for i in range(len(nodes)):
        if nodes[i][0][2]==-1:
            for path in light_loops_with_initial_segment(workingG,[nodes[i],],maxweight,simple=simple):
                if len(path)>2: # because we only are interested in essential vertices
                    yield path
        workingG.remove_node(nodes[i]) # By rebasing, we may assume loop starts at node with least index. Thus, after generating loops starting at node i we only need to continue generating loops that do not pass through node i. 


        
def restricted_light_loops(G,possiblefirstcornerleg1=None,possiblefirstcornerleg2=None,possiblelastcornerleg1=None,possiblelastcornerleg2=None,maxweight=1,simple=False):
    """
    Generator that yields loops in weighted corner graph G whose weight does not exceed maxweight, subject to the constraints that the first vertex is of the form (a,b) where a is in possiblefirstcornerleg1 and b is in possiblefirstcornerleg2 and the last vertex is of the form (c,d) where c is in possiblelastcornerleg2 and d is in possiblelastcornerleg2. If any of the possibility sets is None that means no restriction on that coordinate. 
    """
    if None in [possiblefirstcornerleg1,possiblefirstcornerleg2,possiblelastcornerleg1,possiblelastcornerleg2]:
        allpiecesset={x for (x,y) in G}
        allpieces=sorted(allpiecesset,key=lambda p: -p[3])
    if possiblefirstcornerleg1 is None:
        A=allpieces
    else:
        A=possiblefirstcornerleg1
    if possiblefirstcornerleg2 is None:
        B=allpieces
    else:
        B=possiblefirstcornerleg2
    if possiblelastcornerleg1 is None:
        C=allpieces
    else:
        C=possiblelastcornerleg1
    if possiblelastcornerleg2 is None:
        D=allpieces
    else:
        D=possiblelastcornerleg2
    for possible_start_vertex in [(a,b) for (a,b) in G if a in A and b in B]:
        for possible_end_vertex in [(c,d) for (c,d) in G if c in C and d in D]:
            for path in light_loops_with_initial_segment(G,[possible_end_vertex,possible_start_vertex],maxweight,simple=simple):
                if len(path)>2: # because we only are interested in essential vertices
                    yield path[1:]+path[0:1] # move what was supposed to be possible_end_vertex to the end

def negative_curvature_check(relator_list,angle_function,precomputed_piecedict=None,precomputed_corner_graph=None,noparse=False,verbose=False):
    """
    Check if every interior vertex of every reduced diagram has negative curvature.
    If verbose and False, also return the curvature of the first non-negatively curved vertex found, and the loop in the corner graph describing that vertex.
    angle_function should be a function that takes a corner as input and outputs the angle assigned to that corner. 
    """
    if noparse:
        rels=relator_list
    else:
        if verbose:
            print("Parsing input words.")
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        if verbose:
            print("Computing pieces.")
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    thepiecesegments=piecesegments(rels,thepiecedict)
    if precomputed_corner_graph is None:
        if verbose:
            print("Constructing  corner graph.")
        G=unweighted_corner_graph(rels,thepiecesegments)
        if verbose:
            print("Corner graph has "+str(len(G))+" vertices.")
    else:
        G=precomputed_corner_graph
    if verbose:
        print("Computing angle structure") 
    weight_corner_graph(rels,G,angle_function,precomputed_piecedict=thepiecedict)
    ll=light_loops(G,maxweight=1)
    for loop in ll:
        if verbose:
            return False,1-sum(G.nodes[v].get('weight') for v in loop),loop
        else:
            return False
    else:
        return True

def piece_segment_weight_guesser(rels,precomputed_piecedict=None,initial_ps_weights=None,iterations=10,corner_adjustment_factor=Fraction(1,10),complement_adjustment_factor=Fraction(1,10),number_target_loops=2,precision=3,verbose=False):
    if precomputed_piecedict is None:
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    allpiecesegments=piecesegments(rels,thepiecedict)
    G=unweighted_corner_graph(rels,allpiecesegments)
    if initial_ps_weights is None:
        psweightdict={x:1 for x in allpiecesegments}
    elif callable(initial_ps_weights):
        psweightdict={x:initial_ps_weights(x) for x in allpiecesegments}
    elif type(initial_ps_weights)==dict:
        psweightdict=initial_ps_weights
    def blind_corner_angle(corner):
        complement=corner_remainder_segment(rels,corner)
        assert(complement[3]>0) # presentation is C(3), no bigons 
        return Fraction(1,2)*(1-Fraction(psweightdict[corner[0]]+psweightdict[corner[1]],psweightdict[corner[0]]+psweightdict[corner[1]]+segment_piece_length(rels,complement,thepiecedict,lambda x:psweightdict[x])))
    CAF=blind_corner_angle
    iteration=0
    while iteration<=iterations-1:
        iteration+=1
        weight_corner_graph(rels,G,CAF,thepiecedict)
        ll=light_loops(G,maxweight=1)
        lightest=sorted(ll,key=lambda x:sum(G.nodes[v]['weight'] for v in x))
        corner_segment_frequency=dict()
        complement_segment_frequency=dict()
        num_loops=min(len(lightest),number_target_loops)
        if not num_loops:
            if verbose:
                print("No light loops.")
            return psweightdict
        else:
            if verbose:
                print('Iteration '+str(iteration)+' of '+str(iterations)+'. Lightest loops '+" ".join(f"{float(sum(G.nodes[v]['weight'] for v in lightest[i])):.3g}" for i in range(num_loops))+'.')
        for i in range(num_loops):
            for c in lightest[i]:
                corner_segment_frequency[c[0]]=corner_segment_frequency.setdefault(c[0],0)+Fraction(1,num_loops)
                corner_segment_frequency[c[1]]=corner_segment_frequency.setdefault(c[1],0)+Fraction(1,num_loops)
                complement=corner_remainder_segment(rels,c)
                startvertex=complement[1]
                endvertex=complement[1]+complement[2]*complement[3]
                spg=segment_piece_graph(rels,complement,thepiecedict)
                paths=nx.all_simple_paths(spg,startvertex,endvertex)
                pathcount=0
                segcount=dict()
                for path in paths:
                    pathcount+=1
                    for e in range(len(path)-1):
                        pslabel=spg[path[e]][path[e+1]]['label']
                        segcount[pslabel]=segcount.setdefault(pslabel,0)+1
                for ps in segcount:
                    complement_segment_frequency[ps]=complement_segment_frequency.setdefault(ps,0)+Fraction(segcount[ps],pathcount)*Fraction(1,num_loops)
        newpsweightdict=dict()
        for ps in allpiecesegments:
            if ps not in complement_segment_frequency and ps not in corner_segment_frequency:
                newpsweightdict[ps]=psweightdict[ps]
            elif ps not in complement_segment_frequency:
                newpsweightdict[ps]=psweightdict[ps]*(1-corner_adjustment_factor*Fraction(1,iteration)*min(1,corner_segment_frequency[ps]))
            elif ps not in corner_segment_frequency:
                newpsweightdict[ps]=psweightdict[ps]*(1+complement_adjustment_factor*Fraction(1,iteration)*min(1,complement_segment_frequency[ps]))
            else:
                newpsweightdict[ps]=psweightdict[ps]*(1-corner_adjustment_factor*Fraction(1,iteration)*min(1,corner_segment_frequency[ps]))*(1+complement_adjustment_factor*Fraction(1,iteration)*min(1,complement_segment_frequency[ps]))
        for ps in allpiecesegments:
            psweightdict[ps]=Fraction(round(newpsweightdict[ps]*10**precision),10**precision)
        def blind_corner_angle(corner):
            complement=corner_remainder_segment(rels,corner)
            assert(complement[3]>0) # presentation is C(3), no bigons 
            return Fraction(1,2)*(1-Fraction(psweightdict[corner[0]]+psweightdict[corner[1]],psweightdict[corner[0]]+psweightdict[corner[1]]+segment_piece_length(rels,complement,thepiecedict,lambda x:psweightdict[x])))
        CAF=blind_corner_angle
    return psweightdict
        
        

def generator_weight_guesser(rels,precomputed_piecedict=None,initial_generator_weights=None,iterations=10,corner_adjustment_factor=Fraction(1,10),complement_adjustment_factor=Fraction(1,10),number_target_loops=2,precision=3,verbose=False):
    if precomputed_piecedict is None:
        thepiecedict = piecedict(rels)
    else:
        thepiecedict = precomputed_piecedict
    allpiecesegments = piecesegments(rels, thepiecedict)
    G = unweighted_corner_graph(rels, allpiecesegments)
    maxgen = max(abs(x) for w in rels for x in w)
    if initial_generator_weights is None:
        generator_weights = [1 for x in range(maxgen)]
    elif isinstance(initial_generator_weights, (list, tuple)):
        generator_weights =  initial_generator_weights
    else:
        raise TypeError("initial_generator_weights must be None, list/tuple")
    def CAF(corner):
        return corner_angle_from_generator_weights(rels,corner,generator_weights,precomputed_piecedict=thepiecedict)
    iteration = 0
    while iteration <= iterations-1:
        iteration += 1
        weight_corner_graph(rels, G, CAF, thepiecedict)
        ll = light_loops(G, maxweight=1)
        lightest = sorted(ll, key=lambda loop: sum(G.nodes[v]['weight'] for v in loop))
        num_loops = min(len(lightest), number_target_loops)
        if not num_loops:
            if verbose:
                print("No light loops.")
            return generator_weights
        if verbose:
            print(
                f"Iteration {iteration} of {iterations}. Lightest loops "
                + " ".join(
                    f"{float(sum(G.nodes[v]['weight'] for v in lightest[i])):.3g}"
                    for i in range(num_loops)
                )
                + "."
            )
        leg_freq = {i: Fraction(0,1) for i in range(1, maxgen+1)}
        comp_freq = {i: Fraction(0,1) for i in range(1, maxgen+1)}
        for i in range(num_loops):
            loop = lightest[i]
            for c in loop:
                for x in subword(rels, c[0]):
                    leg_freq[abs(x)] += Fraction(1, num_loops)
                for x in subword(rels, c[1]):
                    leg_freq[abs(x)] += Fraction(1, num_loops)
                complement = corner_remainder_segment(rels, c)
                assert(complement[3]>0)
                frac_unit = Fraction(1, complement[3]*num_loops)
                for x in subword(rels, complement):
                    comp_freq[abs(x)] += frac_unit
        new_weights = list()
        for i in range(1, maxgen+1):
            old = generator_weights[i-1]
            lf_clamped = min(Fraction(1,1), leg_freq[i])
            cf_clamped = min(Fraction(1,1), comp_freq[i])
            if lf_clamped == 0 and cf_clamped == 0:
                factor = 1
            elif lf_clamped == 0:
                factor = (1 + complement_adjustment_factor * Fraction(1, iteration) * cf_clamped)
            elif cf_clamped == 0:
                factor = (1 - corner_adjustment_factor * Fraction(1, iteration) * lf_clamped)
            else:
                factor = (1 - corner_adjustment_factor * Fraction(1, iteration) * lf_clamped) * (1 + complement_adjustment_factor * Fraction(1, iteration) * cf_clamped)
            new_weights.append(Fraction(round(old*factor*10**precision),10**precision))
        generator_weights = new_weights
        def CAF(corner):
            return corner_angle_from_generator_weights(rels,corner,generator_weights,precomputed_piecedict=thepiecedict)
    return generator_weights


def second_small_cancellation_check(relator_list,corner_angle_function=None,additive_piece_weight=False, generator_weights=None,piece_weights=None,piece_segment_weights=None,precomputed_piecedict=None,precomputed_corner_graph=None, noparse=False,verbose=False,minimum_heavy_link_weight=Fraction(1,1),maximum_light_link_weight=Fraction(2,1),guess_weights=False):
    """
    Decide if van Kampen diagrams over C(3) presentation with given relator_list have the property that with the given weights on the lengths of edges corresponding to each generator, all deep vertices in the diagram have average negative curvature. If yes, return True, and the group is hyperbolic. If algorithm fails then False is returned and hyperbolicity of the group is not determined. 
    
    Algorithm works by enumerating possible 2-neighborhoods of vertices in reduced van Kampen diagrams, computing curvatures of central vertex and its neighbors, and, if central vertex is nonnegatively curved, taking donation from each of its negatively curved neighbors in the amount of curvature/degree. If the result is always negative then return True. Return false if a vertex is found that has nonnegative curvature and for which it is likely that a neighborhood can be found for which the neighboring vertices are not able to donate sufficient curvature to make the central vertex negatively curved. This can mean that either an explicit failing neighborhood has been found, or that ruling out such a neighborhood would take an excessively long time. This decision is controlled by parameters minimum_heavy_link_weight and maximum_light_link_weight. Must have 1<=minimum_heavy_link_weight<=maximum_light_link_weight. Smaller values lead to faster execution but potentially more false negative results. 

    By default, relator_list is parsed and put into standard form. Use noparse=True to skip this step if the relator_list is already known to be in standard form.

    If corner_angle_function is supplied, use it to measure corner angles. Otherwise:
        If additive_piece_weight is True:
            If generator_weights are supplied, define corner angles using weighted word length. Otherwise, define corner angles using word length.
        Otherwise, if piece_weights are supplied, use them to define corner angles. If not, use constant piece weight 1.
    
    """
    # in the additive case the corner angles can be computed initially and do not depend on the diagram. They are stored in the weighted corner graph G and only have to be looked up. In the nonadditive case, computing the angle of a corner in a diagram requires a decomposition of a face into pieces, not just  the two pieces at the corner.  In this case, the angle recorded in the corner graph is the lower bound obtained by assuming that the complement of the corner is decomposed into pieces as coarsely as possible. However, when we start enumerating links of neighboring vertices, this puts further constraints on the faces incident to the central vertex that can decrease its curvature. So the central vertex curvature must be recomputed several times in the computation.   
    if noparse:
        rels=relator_list
    else:
        if verbose:
            print("Parsing input words.")
        rels=parseinputwords(relator_list)
    if precomputed_piecedict is None:
        if verbose:
            print("Computing pieces.")
        thepiecedict=piecedict(rels)
    else:
        thepiecedict=precomputed_piecedict
    p,q=CT(rels,precomputed_piecedict=thepiecedict,noparse=True)
    if p<3: # when p=3 some of the partial diagram constructed in this algorithm may not actually be realizable, because we have not checked for triangles in the second shell. 
        if verbose:
            print("This presentation is a C"+str(p)+" presentation. Algorithm requires at least C3.")
        return False
    if Fraction(1,p)+Fraction(1,q)<Fraction(1,2):
        if verbose:
            print("Presentation is hyperbolic C"+str(p)+"-T"+str(q)+".")
        # here we should shortcircuit and return True, but for debug checking we continue
    elif Fraction(1,p)+Fraction(1,q)==Fraction(1,2):
        if verbose:
            print("Presentation is non-positively curved C"+str(p)+"-T"+str(q)+".")
    else:
         if verbose:
            print("Presentation CT does not rule out positively curved vertices: C"+str(p)+"-T"+str(q)+".")
    thepiecesegments=piecesegments(rels,thepiecedict)
    if precomputed_corner_graph is None:
        if verbose:
            print("Constructing corner graph.")
        G=unweighted_corner_graph(rels,thepiecesegments)
        if verbose:
            print("Corner graph has "+str(len(G))+" vertices.")
    else:
        G=precomputed_corner_graph
    if corner_angle_function is not None:
        CAF=corner_angle_function
    elif not additive_piece_weight:
        if piece_segment_weights is not None:
            if callable(piece_segment_weights):
                thepsweights=piece_segment_weights
            elif type(piece_segment_weights)==dict:
                def thepsweights(ps):
                    return piece_segment_weights[ps]
            else:
                raise TypeError("piece_segments_weights should be function or dict.")
        else:
            if piece_weights is None:
                if guess_weights:
                    if verbose:
                        print('Guessing piece-segment weights.')
                    psweightdict=piece_segment_weight_guesser(rels,precomputed_piecedict=thepiecedict,initial_ps_weights=None,iterations=10,corner_adjustment_factor=Fraction(1,30),complement_adjustment_factor=Fraction(1,30),number_target_loops=2,precision=2,verbose=False)
                    thepsweights=lambda x:psweightdict[x]
                else:
                    thepsweights=lambda x:1
            elif callable(piece_weights) or type(piece_weights)==dict:
                if callable(piece_weights):
                    _psweights={ps:piece_weights(subword(rels,ps)) for ps in piecesegments(rels,thepiecedict)}
                else:
                    _psweights={ps:piece_weights[subword(rels,ps)] for ps in piecesegments(rels,thepiecedict)}
                def thepsweights(piece_segment):
                    return _psweights[piece_segment]
            else:
                raise TypeError("piece_weights should be function or dict.")
        _segment_weight_cache=dict()
        def cached_segment_weight(segment):
            if segment in _segment_weight_cache:
                return _segment_weight_cache[segment]
            theweight = segment_piece_length(rels,segment,precomputed_piecedict=thepiecedict, piece_segment_weights=thepsweights)
            _segment_weight_cache[segment] = theweight
            return theweight
        def blind_corner_angle(corner):
            complement=corner_remainder_segment(rels,corner)
            assert(complement[3]>0) # presentation is C(3), no bigons 
            return Fraction(1,2)*(1-Fraction(thepsweights(corner[0])+thepsweights(corner[1]),thepsweights(corner[0])+thepsweights(corner[1])+cached_segment_weight(complement)))
        CAF=blind_corner_angle
    else: # additive piece weights, compute corner angles once and for all as weights in corner graph
        if generator_weights is not None:
            CAF=functools.partial(corner_angle_from_generator_weights,rels,generator_weights=weights,precomputed_piecedict=thepiecedict)
        elif guess_weights:
            genweights=generator_weight_guesser(rels,precomputed_piecedict=thepiecedict,iterations=10,corner_adjustment_factor=Fraction(1,30),complement_adjustment_factor=Fraction(1,30),number_target_loops=2,precision=3,verbose=False)
            if verbose:
                print('Guessed generator weights: '+str(genweights))
            CAF=functools.partial(corner_angle_from_generator_weights,rels,generator_weights=genweights,precomputed_piecedict=thepiecedict)
        else:
            CAF=functools.partial(corner_angle_from_generator_weights,rels,generator_weights=[1 for i in range(max([max([abs(x) for x in rel]) for rel in rels]))], precomputed_piecedict=thepiecedict)
    if verbose:
        print("Computing angle structure.")
    weight_corner_graph(rels,G,CAF,precomputed_piecedict=thepiecedict)
    the_light_links=light_loops(G,maxweight=1)
    try:
        first_light_link=next(the_light_links)
        LL=itertools.chain([first_light_link], the_light_links)
    except StopIteration:
        if verbose:
            print("All vertices have negative curvature.")
            return (True,)
        else:
            return True
    if verbose:
        print("Analyzing angle structure.")
    minangle=min(G.nodes[v]['weight'] for v in G)
    mindensity=min(Fraction(G.edges[e]['weight']+G.edges[f]['weight']+G.edges[g]['weight']+G.edges[h]['weight'],4) for e in G.edges() for f in G.out_edges(e[1]) for g in G.out_edges(f[1]) for h in G.out_edges(g[1])) # calculate minimum of (1/4)(path weight) over paths of length 4. This is lower bound for denisty=total weght/length of closed loops. 
    if verbose:
        print("Smallest angle is "+str(minangle)+".")
        print("Lower bound on loop density is "+str(mindensity)+".")
    if mindensity<=0:
        raise ValueError('Corner graph has loops with nonpositive weight.') # Algorithm doesn't work with nonpositive densities. This can happen for some choices of nonpositive weights. Should not happen otherwise. 
    if verbose:
        print("Searching for nonnegative curvature.")
    for centerlink in LL: # each possible centerlink describes link of vertex with nonnegative curvature
        faceangles=[G.nodes[v].get('weight') for v in centerlink]
        centerangle=sum(faceangles)
        centercurvature=1-centerangle
        if centercurvature>=len(centerlink)*mindensity: # can only expect asymptotically that neighbors donate -mindensity, so if this test is true then the center is too positively curved for this algorithm to cancel it out via neighbor donations.
            if verbose:
                return False,"link"+str(centerlink)+"has curvature "+str(centerurvature)+", too small to be balanced by neighbors."
            else:
                return False
        if verbose:
            print("Found link of curvature +"+str(centercurvature)+", based on face angles "+str(faceangles)+". Soliciting donations.")
        neighbors=[]
        # neighbors is a list whose entry at index i describes the neighborhood of the vertex opposite the central vertex along the arc common to corners i and i+1 of centerlink. 
        # entry is a dict:
        #{
        #'firstcorner': first corner of the link clockwise starting from the central arc,
        #'firstcornerangle': current lower bound on angle of the first corner,
        #'lastcorner': last corner of the link clockwise starting from the central arc,
        #'lastcornerangle': current lower bound on angle of the last corner,
        #'firstlastgen': generator of other possible first/last corner choices, in case we backtrack and replace the current first/last corner pair
        #'thislightlink':vertex link that starts with firstcorner, ends with lastcorner and whose weight is at most lightlinkcompletionweightlimit, or None if no such light links exist.
        #'lightlinkcompletionweightlimit': used to define generator of light links. 
        #'lightlinkgen':generator of other possible light links for this first/last pair or None if there are no such,
        #'donation':curvature donation of this vertex, which is either min(0,(1-linkweight)/valence) if there is a light link, or a commonly computed negative upper bound that applies to all heavy links with this first/last pair].

        if not additive_piece_weight:
            def corner_angles_from_face_corner_sequence(thisfacecornersequence):
                faceinnerwordlength=sum(c[0][3] for c in thisfacecornersequence)+thisfacecornersequence[-1][1][3]
                relatorindex=thisfacecornersequence[0][0][0]
                assert(all(relatorindex==c[0][0] and relatorindex==c[1][0] for c in thisfacecornersequence))
                facetotalwordlength=len(rels[relatorindex])
                if facetotalwordlength>=faceinnerwordlength:
                    thisfaceistriangle=False
                    if thisfacecornersequence[0][0][2]==1:
                        complementarysegment=(relatorindex,(thisfacecornersequence[0][0][1]+thisfacecornersequence[0][0][3])%len(rels[relatorindex]),1,facetotalwordlength-faceinnerwordlength)
                    else:
                        complementarysegment=(relatorindex,(thisfacecornersequence[0][0][1]-thisfacecornersequence[0][0][3])%len(rels[relatorindex]),-1,facetotalwordlength-faceinnerwordlength)
                elif len(thisfacecornersequence)==3 and thisfacecornersequence[-1][1][3]==thisfacecornersequence[0][0][3] and facetotalwordlength==sum(c[0][3] for c in thisfacecornersequence):
                        thisfaceistriangle=True
                else:
                    raise RuntimeError # should not have allowed such a choice of corners
                if thisfaceistriangle:
                    faceboundaryweight=sum(thepsweights(c[0]) for c in thisfacecornersequence)
                else:
                    faceboundaryweight=sum(thepsweights(c[0]) for c in thisfacecornersequence)+thepsweights(thisfacecornersequence[-1][1])+cached_segment_weight(complementarysegment)
                return [Fraction(1,2)*(1-Fraction(thepsweights(c[0])+thepsweights(c[1]),faceboundaryweight)) for c in thisfacecornersequence]
                
            def recomputecurvatures(verbose=False):
                # In the non-additive case we need to know the piece decomposition of a face to compute the angles of its corners. We compute lower bounds on angles by taking the known pieces in the boundary of the face and taking the coarsest piece decomposition of the remainder segment. As the diagram is extended, more of the actual piece decomposition of the faces becomes determined, so the possibilites for piece decompositions are further constrained. This can lead to better lower bounds for corner angles, hence, lower curvatures. This function recomputes all the curvature bounds given the current state of the diagram. 
                faceangles=[]#recomputed angles of the central corners
                if len(neighbors)==0: # no neighbors yet, use angles from corner graph
                    return 1-sum(G.nodes[corner]['weight'] for corner in centerlink)
                elif len(neighbors)==len(centerlink): # links for all neighbors are chosen
                    for i in range(len(centerlink)):
                        thisfacecornersequence=[neighbors[i-1]['lastcorner'],centerlink[i],neighbors[i]['firstcorner']]
                        thisfaceangles=corner_angles_from_face_corner_sequence(thisfacecornersequence)
                        
                        faceangles.append(thisfaceangles[1])
                        neighbors[i-1]['lastcornerangle']=thisfaceangles[0]
                        neighbors[i]['firstcornerangle']=thisfaceangles[2]
                    for i in range(len(centerlink)): # recompute curvature donations
                        if 'thislightlink' in neighbors[i] and neighbors[i]['thislightlink'] is not None:
                            neighbors[i]['donation']=min(0,Fraction(1-neighbors[i]['firstcornerangle']-neighbors[i]['lastcornerangle']-sum(G.nodes[c]['weight'] for c in neighbors[i]['thislightlink'][1:-1]),len(neighbors[i]['thislightlink'])))
                        elif i<len(neighbors)-1: # we have not reached the end of neighbors yet, so for this neighbor light links should have already been defined, and we have exhausted them. This is a heavy link.                                 
                            neighbors[i]['donation']=-mindensity*(1-Fraction(1,neighbors[i]['lightlinkcompletionweightlimit']+neighbors[i]['firstcornerangle']+neighbors[i]['lastcornerangle'])) # heavy link donation
                        else: # we have reached last neighbor, and either 'thislightlink' is undefined or None
                            assert(i==len(neighbors)-1)
                            if 'thislightlink' not in neighbors[-1]:
                                 neighbors[-1]['donation']=0
                            if 'lightlinkcompletionweightlimit' not in neighbors[-1]:
                                requireddonation=min(0,-(centercurvature+sum(neighbors[i]['donation'] for i in range(len(neighbors)-1))))
                                mw=Fraction(mindensity,mindensity+requireddonation)
                                neighbors[-1]['lightlinkcompletionweightlimit']=mw-faceangles[0]-faceangles[-1]
                                neighbors[-1]['donation']=-mindensity*(1-Fraction(1,neighbors[-1]['lightlinkcompletionweightlimit']+neighbors[-1]['firstcornerangle']+neighbors[-1]['lastcornerangle'])) # heavy link donation
                else: # links for some, but not all, neighbors have been chosen
                    # set new face angles around center vertex
                    # corner 0, last neighbor not chosen yet, only use neighbor[0]
                    thisfacecornersequence=[centerlink[0],neighbors[0]['firstcorner']]
                    thisfaceangles=corner_angles_from_face_corner_sequence(thisfacecornersequence)
                    faceangles.append(thisfaceangles[0])
                    neighbors[0]['firstcornerangle']=thisfaceangles[1]
                    # remaining corners for which at least firstlast of link are chosen
                    for i in range(1,len(neighbors)):
                        thisfacecornersequence=[neighbors[i-1]['lastcorner'],centerlink[i],neighbors[i]['firstcorner']]
                        thisfaceangles=corner_angles_from_face_corner_sequence(thisfacecornersequence)
                        faceangles.append(thisfaceangles[1])
                        neighbors[i-1]['lastcornerangle']=thisfaceangles[0]
                        neighbors[i]['firstcornerangle']=thisfaceangles[2]
                    # first corner for which link has not been chosen, still gets some info from previous
                    thisfacecornersequence=[neighbors[-1]['lastcorner'],centerlink[len(neighbors)]]
                    thisfaceangles=corner_angles_from_face_corner_sequence(thisfacecornersequence)
                    faceangles.append(thisfaceangles[1])
                    neighbors[-1]['lastcornerangle']=thisfaceangles[0]
                    # rest of the corners
                    for i in range(1+len(neighbors),len(centerlink)):
                        faceangles.append(G.nodes[centerlink[i]]['weight']) # this weight already stored in corner graph
                    
                    # set donation for neighbors whose first and last corners are set
                    for i in range(len(neighbors)):
                        if 'thislightlink' in neighbors[i] and neighbors[i]['thislightlink'] is not None:
                            neighbors[i]['donation']=min(0,Fraction(1-neighbors[i]['firstcornerangle']-neighbors[i]['lastcornerangle']-sum(G.nodes[c]['weight'] for c in neighbors[i]['thislightlink'][1:-1]),len(neighbors[i]['thislightlink'])))# donation is angle deficiency/valence
                        elif i<len(neighbors)-1: # heavy link
                            neighbors[i]['donation']=-mindensity*(1-Fraction(1,neighbors[i]['lightlinkcompletionweightlimit']+neighbors[i]['firstcornerangle']+neighbors[i]['lastcornerangle'])) # heavy link donation
                        elif i==len(neighbors)-1 and 'thislightlink' not in neighbors[-1]: # haven't set thislightlink yet for this neighbor
                            neighbors[-1]['donation']=0
                        else: # heavy link
                            assert(i==len(neighbors)-1)
                            if 'lightlinkcompletionweightlimit' not in neighbors[-1]:
                                howmanyafterthis=len(centerlink)-len(neighbors)
                                requireddonation=min(0,-(centercurvature+sum(neighbors[i]['donation'] for i in range(len(neighbors)-1))-mindensity*howmanyafterthis)) # this is donation needed from this vertex to get to 0 if all of the remaining only donate -mindensity
                                mw=Fraction(mindensity,mindensity+requireddonation)
                                maxweighttocheck=max(heavyparameter,mw) # every link with weight > maxweighttocheck is guaranteed to dontate requireddontation. Proof requires separate estimates for case that length of the link is at least maxweighttocheck/mindensity or less than that.
                                neighbors[-1]['lightlinkcompletionweightlimit']=maxweighttocheck-neighbors[-1]['firstcornerangle']-neighbors[-1]['lastcornerangle']
                                neighbors[-1]['donation']=-mindensity*(1-Fraction(1,neighbors[-1]['lightlinkcompletionweightlimit']+neighbors[-1]['firstcornerangle']+neighbors[-1]['lastcornerangle'])) # heavy link donation
                centerangle=sum(faceangles)
                newcentercurvature=1-centerangle
                if verbose:
                    print('recomputed center curvature',newcentercurvature,'new face angles',faceangles )
                return newcentercurvature
            ############### end of recomputecurvature
            
        # Iterate through the neighbors of the centeral vertex. For each of them we make a generator that yields links of that vertex that agrees with link of central vertex and with preceding neighbor. The activeneighbor is the one we are currently working on.
        activeneighbor=0
        while activeneighbor>=0:
            if len(neighbors)==activeneighbor or 'lightlinkgen' not in neighbors[activeneighbor]: # either there is nothing in this entry yet, or there is a choice of first/last outgoing arc, but no link completions
                if len(neighbors)==activeneighbor: # no entry yet, create and populate the dict
                    firstcornerfirstleg=reverse_segment(rels,centerlink[activeneighbor][1])
                    lastcornersecondleg=reverse_segment(rels,centerlink[(activeneighbor+1)%len(centerlink)][0])
                    if activeneighbor==0: # this is the first neighbor of the central vertex whose link we are trying to complete. Its first corner belongs to face with two sides already specified by centerlink[0][0] and centerlink[0][1]. firstcornerfirstleg is already computed. second leg can be any successor that is short enough with respect to two existing sides. 
                        firstboundarycomplementlength=len(rels[centerlink[0][0][0]])-centerlink[0][0][3]-centerlink[0][1][3]
                        assert(firstboundarycomplementlength>0) # because p>2
                        possiblefirstcorners=((a,b) for (a,b) in G if a==firstcornerfirstleg and b[3]<=firstboundarycomplementlength)
                    else: # not the first neighbor, so cell containing the first corner of this vertex link already has 3 of its faces decided by the centerlink and the last corner of the previous neighbor
                        firstboundarycomplementlength=len(rels[centerlink[activeneighbor][0][0]])-neighbors[activeneighbor-1]['lastcorner'][0][3]-centerlink[activeneighbor][0][3]-centerlink[activeneighbor][1][3]
                        assert(firstboundarycomplementlength>0 or (firstboundarycomplementlength==0 and p==3))
                        if firstboundarycomplementlength==0: # triangular face potentially happens in a C(3) presentation
                            possiblefirstcorners=[(firstcornerfirstleg,reverse_segment(rels,neighbors[activeneighbor-1]['lastcorner'][0])),]
                        else:
                            possiblefirstcorners=((a,b) for (a,b) in G if a==firstcornerfirstleg and b[3]<=firstboundarycomplementlength)
                    if activeneighbor==len(centerlink)-1: # this is the last neighbor of the central vertex. extra care here because the last corner of this link is influenced by the first corner of the link of neighbor 0.
                        lastboundarycomplementlength=len(rels[centerlink[0][0][0]])-centerlink[0][0][3]-centerlink[0][1][3]-neighbors[0]['firstcorner'][1][3]
                        assert(lastboundarycomplementlength>0 or (lastboundarycomplementlength==0 and p==3))
                        if lastboundarycomplementlength==0:  # triangular face potentially happens in a C(3) presentation
                            possiblelastcorners=[(reverse_segment(rels,neighbors[0]['firstcorner'][1]),lastcornersecondleg),]
                        else:
                            possiblelastcorners=((c,d) for (c,d) in G if d==lastcornersecondleg and c[3]<=lastboundarycomplementlength)
                    else: # some neighbor in the middle
                        lastboundarycomplementlength=len(rels[centerlink[(activeneighbor+1)%len(centerlink)][0][0]])-centerlink[(activeneighbor+1)%len(centerlink)][0][3]-centerlink[(activeneighbor+1)%len(centerlink)][1][3]
                        possiblelastcorners=((c,d) for (c,d) in G if d==lastcornersecondleg and c[3]<=lastboundarycomplementlength)
                    firstlastgen=itertools.product(possiblefirstcorners,possiblelastcorners) 
                    thisfirstlast=next(firstlastgen) # If the presentation is C(2) there should be at least one such possibility, so should not get a StopIteration here.
                    neighbors.append({'firstcorner':thisfirstlast[0],'lastcorner':thisfirstlast[1],'firstlastgen':firstlastgen})
                    if not additive_piece_weight: # we added new neighbor corners, so recompute curvatures, and return new curvature of center vertex
                        centercurvature=recomputecurvatures()
                        if verbose>=2:
                            print("active neighbor "+str(activeneighbor)+". Recompute curvatures")
                        backtrack=False
                        while centercurvature<0 and not backtrack:# center curvature came out negative, so we're happy with any diagram that agrees with this partial description. Move on; try different choices of first and last corner at activeneighbor. 
                            try:
                                thisfirstlast=next(neighbors[activeneighbor]['firstlastgen'])
                                neighbors[activeneighbor]['firstcorner']=thisfirstlast[0]
                                neighbors[activeneighbor]['lastcorner']=thisfirstlast[1]
                                centercurvature=recomputecurvatures()
                                if verbose>=2:
                                    print("active neighbor "+str(activeneighbor)+". Recompute curvatures")
                                # this also sets neighbors[activeneighbor]['firstcornerangle'] and neighbors[activeneighbor]['lastcornerangle']
                            except StopIteration:
                                backtrack=True
                        if backtrack:
                            neighbors=neighbors[:activeneighbor]
                            activeneighbor-=1
                            if verbose>=2:
                                print('backtrack to activeneighbor '+str(activeneighbor))
                            continue
                    else: # in the additive case adding neighbor corners does not change computation of centercurvature
                        neighbors[activeneighbor]['firstcornerangle']=G.nodes[neighbors[activeneighbor]['firstcorner']]['weight']
                        neighbors[activeneighbor]['lastcornerangle']=G.nodes[neighbors[activeneighbor]['lastcorner']]['weight']

                # find link completions of this neighbor
                howmanyafterthis=len(centerlink)-activeneighbor-1
                requireddonation=min(0,-(centercurvature+sum(neighbors[i]['donation'] for i in range(activeneighbor))-mindensity*howmanyafterthis)) # this is donation needed from this vertex to get to 0 if all of the remaining only donate -mindensity
                if mindensity<-requireddonation:
                    assert(False) # this is a fail state, but should have been caught when the previous neighbor link was evaluated. Something is wrong. 
                elif mindensity==-requireddonation: # this is a fail state. The donation of a minimial density link is (1-weight)/length= -mindensity + 1/length, so dontation are always slightly less than mindensity no matter how heavy the link. 
                    if verbose:
                        return False,centercurvature,newcentercurvature, centerlink, neighbors  ########## maybe some text would be in order to explain this output
                    else:
                        return False
                mw=Fraction(mindensity,mindensity+requireddonation)
                maxweighttocheck=max(minimum_heavy_link_weight,mw) # every link with weight > maxweighttocheck is guaranteed to dontate requireddontation. Proof requires separate estimates for case that length of the link is at least maxweighttocheck/mindensity or less than that.
                if maxweighttocheck>maximum_light_link_weight: # search takes too long if demand to check all possible links with weights that go past maximum_light_link_weight. This can easily happen for nonhyperbolic examples. 
                    if verbose:
                        return False,centercurvature,newcentercurvature, centerlink, neighbors  ########## maybe some text would be in order to explain this output
                    else:
                        return False
                        
                neighbors[activeneighbor]['lightlinkcompletionweightlimit']=maxweighttocheck-neighbors[activeneighbor]['firstcornerangle']-neighbors[activeneighbor]['lastcornerangle']
                possiblelightlinks=restricted_light_loops(G,possiblefirstcornerleg1=set([neighbors[activeneighbor]['firstcorner'][0],]),possiblefirstcornerleg2=set([neighbors[activeneighbor]['firstcorner'][1],]),possiblelastcornerleg1=set([neighbors[activeneighbor]['lastcorner'][0],]),possiblelastcornerleg2=set([neighbors[activeneighbor]['lastcorner'][1],]),maxweight=neighbors[activeneighbor]['lightlinkcompletionweightlimit']+G.nodes[neighbors[activeneighbor]['firstcorner']]['weight']+G.nodes[neighbors[activeneighbor]['lastcorner']]['weight']) # this generates links starting with first and ending with last whose total weight is at most maxweighttocheck. All other potential link completions are 'heavy', and we have already given upper bound for their donation. 
                try:
                    thislightlink=next(possiblelightlinks)
                    thislightlinkcurvature=1-sum(G.nodes[v]['weight'] for v in thislightlink[1:-1])-neighbors[activeneighbor]['firstcornerangle']-neighbors[activeneighbor]['lastcornerangle']
                    thislightlinkdonation=min(0,Fraction(thislightlinkcurvature,len(thislightlink)))
                    neighbors[activeneighbor]['thislightlink']=thislightlink
                    neighbors[activeneighbor]['lightlinkgen']=possiblelightlinks
                    neighbors[activeneighbor]['donation']=thislightlinkdonation
                    newcentercurvature=centercurvature+sum(neighbors[i]['donation'] for i in range(activeneighbor+1))
                    usedaheavylink=any(neighbors[i]['thislightlink'] is None for i in range(activeneighbor+1))
                    if newcentercurvature<0 or (newcentercurvature==0 and usedaheavylink): # we have enough donated curvature already to make the center vertex negatively curved, no matter what happens with the rest of the neighborhood. No need to continue; must look elsewhere for bad diagrams.
                        pass # do not change activeneighbor, next time through the loop will change the choice of link completing this first/last choice
                    elif newcentercurvature>(len(centerlink)-activeneighbor-1)*mindensity or ( newcentercurvature==(len(centerlink)-activeneighbor-1)*mindensity and not usedaheavylink): # we have fallen behind to the point where it will not be possible to get the center vertex to negative curvature by completing this diagram neighborhood with heavy links.
                        if verbose:
                            return False,centercurvature,newcentercurvature, centerlink, neighbors  ########## maybe some text would be in order to explain this output
                        else:
                            return False
                    else:
                        activeneighbor+=1
                except StopIteration: # no way to complete to light link
                    neighbors[activeneighbor]['thislightlink']=None
                    neighbors[activeneighbor]['lightlinkgen']=None
                    neighbors[activeneighbor]['donation']=-mindensity*(1-Fraction(1,neighbors[activeneighbor]['lightlinkcompletionweightlimit']+neighbors[activeneighbor]['firstcornerangle']+neighbors[activeneighbor]['lastcornerangle'])) # heavy link donation
                    if not additive_piece_weight:
                        centercurvature=recomputecurvatures()
                        if verbose>=2:
                            print("active neighbor "+str(activeneighbor)+". Recompute curvatures")
                    newcentercurvature=centercurvature+sum(neighbors[i]['donation'] for i in range(activeneighbor+1))
                    usedaheavylink=True
                    if activeneighbor==len(centerlink)-1: # found links for all neighbors
                        assert(newcentercurvature<=0) # this is a success case, backtrack to continue searching for fail cases
                        pass # do not change activeneighbor. The next iteration of the main loop will increment the first/last choice
                    else:
                        activeneighbor+=1
            else: #  Increment the choice of link completion, or, if current link completion is already heavy, increment the choice of fistlast.
                neighbors=neighbors[:1+activeneighbor]# if we got here by backtracking, throw away all forward indices
                if neighbors[activeneighbor]['thislightlink'] is not None: # current completion with a light link, try for another
                    try:
                        thislightlink=next(neighbors[activeneighbor]['lightlinkgen'])
                        thislightlinkcurvature=1-sum(G.nodes[v]['weight'] for v in thislightlink[1:-1])-neighbors[activeneighbor]['firstcornerangle']-neighbors[activeneighbor]['lastcornerangle']
                        thislightlinkdonation=min(0,Fraction(thislightlinkcurvature,len(thislightlink)))
                        neighbors[activeneighbor]['thislightlink']=thislightlink
                        neighbors[activeneighbor]['donation']=thislightlinkdonation
                        newcentercurvature=centercurvature+sum(neighbors[i]['donation'] for i in range(activeneighbor+1))
                        usedaheavylink=any(neighbors[i]['thislightlink'] is None for i in range(activeneighbor+1))
                        if newcentercurvature<0 or (newcentercurvature==0 and usedaheavylink): # we have enough donated curvature already to make the center vertex negatively curved, no matter what happens with the rest of the neighborhood. No need to continue; must look elsewhere for bad diagrams.
                            pass # do not change activeneighbor, next time through the loop will change the choice of link completing this first/last choice
                        elif newcentercurvature>=(len(centerlink)-activeneighbor-1)*mindensity:# or ( newcentercurvature==(len(centerlink)-activeneighbor-1)*mindensity and not usedaheavylink): # we have fallen behind to the point where it will not be possible to get the center vertex to negative curvature by completing this diagram neighborhood with heavy links.
                            if verbose:
                                return False,centercurvature,newcentercurvature, centerlink, neighbors  ########## maybe some text would be in order to explain this output
                            else:
                                return False
                        else:
                            activeneighbor+=1
                    except StopIteration: # no way to complete to light link, use a heavy one
                        neighbors[activeneighbor]['thislightlink']=None
                        neighbors[activeneighbor]['lightlinkgen']=None
                        neighbors[activeneighbor]['donation']=-mindensity*(1-Fraction(1,neighbors[activeneighbor]['lightlinkcompletionweightlimit']+neighbors[activeneighbor]['firstcornerangle']++neighbors[activeneighbor]['lastcornerangle'])) # heavy link donation
                        newcentercurvature=centercurvature+sum(neighbors[i]['donation'] for i in range(activeneighbor+1))
                        usedaheavylink=True
                        if activeneighbor==len(centerlink)-1: # found links for all neighbors
                            assert(newcentercurvature<=0) # this is a success case, backtrack to continue searching for fail cases
                            pass # do not change activeneighbor. The next iteration of the main loop will increment the first/last choice
                        else:
                            activeneighbor+=1
                else: # current completion is with a heavy link. Inrement first/last choice.
                    foundnext=False
                    try:
                        thisfirstlast=next(neighbors[activeneighbor]['firstlastgen'])
                        foundnext=True
                    except StopIteration:
                        activeneighbor-=1 # no more choices of firstlast at this neighbor. Since we didn't find any fail states, backtrack and loop. Since len(neighbors[activeneighbor])==5 we will take second branch and inrement the link completion, or, if activeneighbor==-1 we will be done
                        if verbose>=2:
                            print("backtrack to activeneighbor "+str(activeneighbor))
                    if foundnext:
                        neighbors[activeneighbor]['firstcorner']=thisfirstlast[0]
                        neighbors[activeneighbor]['lastcorner']=thisfirstlast[1]
                        neighbors[activeneighbor]={k:neighbors[activeneighbor][k] for k in {'firstlastgen', 'firstcorner','lastcorner'}}
                        if not additive_piece_weight: # in this case the new choice of firstlast can change the angles defined for center and previous neighbor. recompute
                            centercurvature=recomputecurvatures()
                            if verbose>=2:
                                print("active neighbor "+str(activeneighbor)+". Recompute curvatures")
                            backtrack=False
                            while centercurvature<0 and not backtrack:
                                try:
                                    thisfirstlast=next(neighbors[activeneighbor]['firstlastgen'])
                                    neighbors[activeneighbor]['firstcorner']=thisfirstlast[0]
                                    neighbors[activeneighbor]['lastcorner']=thisfirstlast[1]
                                    centercurvature=recomputecurvatures()
                                    if verbose>=2:
                                        print("active neighbor "+str(activeneighbor)+". Recompute curvatures")
                                except StopIteration:
                                    backtrack=True
                            if backtrack:
                                neighbors=neighbors[:activeneighbor]
                                activeneighbor-=1
                                if verbose>=2:
                                    print("backtrack to activeneighbor "+str(activeneighbor))
                                continue
                        else:
                            neighbors[activeneighbor]['firstcornerangle']=G.nodes[neighbors[activeneighbor]['firstcorner']]['weight']
                            neighbors[activeneighbor]['lastcornerangle']=G.nodes[neighbors[activeneighbor]['lastcorner']]['weight']                    
        assert(activeneighbor==-1)
        # We have now exited the main loop.
        # This is because we have explored all diagram neighborhoods without finding a fail state in which the algorithm does not succeed in giving the central vertex negative curvature.
        if verbose:
            print("Received enough negative curvature from neighbors.")
    if verbose:
        return (True,)
    else:
        return True
        
                    
                    
   
        
                              
                              
    

        

            

# basics of pieces and piecelength
        
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

def piece_length(theword,thepieces,quit_at=float('inf')):
    """
    Calculate minimal length of the given input string theword as a concatenation of thepieces. 

    Returns float('inf') if the input string cannot be written as a concatentation of thepieces. 
    """
    shortest_expression=shortest_piece_expression(theword,thepieces,quit_at)
    return len(shortest_expression) if shortest_expression is not None else float('inf')


def shortest_piece_expression(theword,thepieces,quit_at=float('inf'),as_cyclic_word=False):
    """
    Recursive determination of shortest expression of string theword as concatentation of thepieces.
  
    Return None if no such expression exists.

    Returned expression is guaranteed to be shortest possible, but is not necessarily unique piece expression of this length. 

    If as_cyclic_word=True then theword is treated as cyclic word, so allow the possibility that a minimal expression is actually an expression of a cyclic permutation of theword. 
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
               
def all_piece_expressions(theword,thepieces, as_cyclic_word=False, quit_at=float('inf')):
    """
    Recursively yield lists of words of length at most quit_at whose elements are in thepieces and whose concatenation is either theword, or, when as_cyclic_word=True, a cyclic permutation of theword. 
    When as_cyclic_word=True the yielded lists are normalized so that the first element of the list has a (possibly non-proper) suffix that agrees with a prefix of relator. 

    >>> [ expr for expr in all_piece_expressions('aba',['aa','b'],as_cyclic_word=True)]
    [['aa', 'b']]
    >>> [ expr for expr in all_piece_expressions('aba',['aa','b'],as_cyclic_word=False)]
    []
    >>> sorted([ expr for expr in all_piece_expressions('aaba',['a','aa','b'],as_cyclic_word=True)])
    [['a', 'a', 'b', 'a'], ['aa', 'a', 'b'], ['aa', 'b', 'a']]
    >>> sorted([ expr for expr in all_piece_expressions('aaba',['a','aa','b'],as_cyclic_word=True,quit_at=3)]) # only expressions of length <=3
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
        



 
########## Utility functions that are not particular to small cancellation
    
def shortest_cycle_length(inputgraph,at_vertex=None, first_edge=None,weight=None,nobigon=False):
    """
    Return the length of the shortest directed cycle in the graph, or float('inf') if none exists. 
    If weight=string (typically string='weight')  is given then use edge attribte string as edge lengths. If weight=None then all edges have length 1.

    If one of at_vertex = v or first_edge=e then give the length of the shortest cycle that either starts at v, or has first edge e, respectively. 

    If nobigon=True do not count cycles of only 2 edges. 
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
        if nobigon:
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
        if nobigon and reverse_edge_is_present:
            shortest_allowed_cycle_using_e=min(shortest_simple_cycle_using_e,4) # bigon is not allowed, but double cover of bigon is allowed.
            if weight is None:
                G.add_edge(e[1],e[0])
            else:
                G.add_edge(e[1],e[0],weight=reverse_edge_weight)
        else:
            shortest_allowed_cycle_using_e=shortest_simple_cycle_using_e       
        shortestcyclelength=min(shortestcyclelength,shortest_allowed_cycle_using_e)
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
    
            
    
