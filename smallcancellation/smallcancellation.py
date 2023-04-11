from collections import deque
from fractions import Fraction
import itertools
import networkx as nx
import math

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
    theT=T(rels,noparse=True)
    Cest=int(math.ceil(Fraction(theCprimebound.denominator,theCprimebound.numerator))) #C'(1/L) => C(L+1), quick check without computing C value
    if (Cest>=5 and theT>=4) or (Cest>=4 and theT>=5) or (Cest>=3 and theT>=7):
        return True
    theC=C(rels,7) # sometimes the C value is better than the estimate given by the C' value, compute it for real
    if (theC>=7) or (theC>=5 and theT>=4) or (theC>=4 and theT>=5) or (theC>=3 and theT>=7):
        return True
    else:
        return False

def Cprimebound(relatorlist,Lambda=1, noparse=False):
    """
    The largest ratio of piece length to length of relator containing it.

    Stop and return 1 if we find any such ratio >= 1/Lambda.

    >>> Cprimebound(['abABcdCD'])
    Fraction(1, 8)
    >>> Cprimebound(['aCbcABCac'])
    Fraction(2, 9)
    """
    if noparse:
        rels=relatorlist
    else:
        rels=parseinputwords(relatorlist,asrelatorlist=True)
    biggestratio=Fraction(1,min(len(r) for r in rels))
    if biggestratio>=Fraction(1,Lambda):
        return 1
    rels.sort(key=len) # sort list of relators with shortest first
    irels=[intlisttostring(rel) for rel in itertools.chain.from_iterable(zip(rels,[inverse(w) for w in rels]))] # arrange relators and inverses in a list of the form relator1, inverse of relator1, relator2, inverse of relator2,...
    # irels are in string representaiton
    drels=[x+x for x in irels] # double the relators to look for pieces that would have wrapped
    for relatorindex in range(len(rels)):
        relator=irels[2*relatorindex]
        foundbiggest=False
        for L in range(len(relator),int(biggestratio*len(relator)),-1):# only check subwords of length L that would give biggest ratio if they were a piece
            for startingindex in range(len(relator)):
                p=(relator+relator)[startingindex:startingindex+L] # the subword of length L starting at index i in relator as a cyclic word
                # now we need to check if p is a piece
                # we do not need to check lower relatorindices, because we already scanned those relators for pieces
                if any(p in x for x in [(relator+relator)[startingindex+1:len(relator)+startingindex+L-1]]+[drels[i] for i in range(2*relatorindex+1,len(drels))]):# look in this (doubled) relator at higher starting indices, and in all later relators, for other copies of p. If found a matching subword, p is a piece.
                    biggestratio=Fraction(len(p),len(relator))
                    foundbiggest=True # we search pieces by decreasing length, so first one we find is longest
                    if biggestratio>=Fraction(1,Lambda):
                        return 1
                    break
            if foundbiggest: # if we found the biggest piece in this relator we can move on to the next relator. 
                break
    return biggestratio

def T(relatorlist,noparse=False):
    """
    Find the minimum degree of an essential interior vertex in a reduced van Kampen diagram.

    >>> T([[1,2,1,2,2,1,2,2,2]]) 
    inf
    >>> T(['cacbcbcbcabacbcaba'])
    6
    >>> T(['aB','bC','cA'])
    3
    """
    # equal to the shortest embedded cycle in the reduced Whitehead graph
    if noparse:
        rels=relatorlist
    else:
        rels=parseinputwords(relatorlist,asrelatorlist=True)
    G=simple_Whitehead_graph(rels)
    theedges=[e for e in G.edges()]
    shortestcycle=float('inf')
    for e in theedges:
        G.remove_edge(*e)
        try:
            shortestcycleusing_e=1+nx.shortest_path_length(G,*e) # compute distance between endpoints of e in G-e then add 1
        except nx.NetworkXNoPath:
            shortestcycleusing_e=float('inf')
        G.add_edge(*e)
        shortestcycle=min(shortestcycle,shortestcycleusing_e)
    return shortestcycle


    
def C(relatorlist,quit_at=float('inf'),piece_up_to_automorphism=False,precomputed_pieces=None,noparse=False):
    """
    FInd the minimum number p such that there exists some cyclic permutation of some relator that can be expressed as a freely reduced product of p pieces.

    If quit_at=q is specified the algorithm will stop and return q once it is determined that p>=q.

    If piece_up_to_automorphism=True then a word only counts as a pieces if it occurs in places in the relator list that are distinct up to automorphism. This means that copies of the root word in a relation that is a proper power do not yield pieces. 
    """
    if noparse:
        rels=relatorlist
    else:
        rels=parseinputwords(relatorlist,asrelatorlist=True)
    if precomputed_pieces is None:
        thepieces=pieces(rels,piece_up_to_automorphism,noparse=False)
    else:
        thepieces=precomputed_pieces
    relators_as_strings=[intlisttostring(rel) for rel in rels]
    pieces_as_strings=[intlisttostring(piece) for piece in thepieces]
    return min((piece_length(therelator,pieces_as_strings,quit_at) for therelator in relators_as_strings))

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
        
        
def pieces(relatorlist,piece_up_to_automorphism=False,noparse=False,asstring=False):
    """
    Given input container of relators, return set of pieces, which are subwords occuring more than once in relators or their inverses, as cyclic words.

    If piece_up_to_automorphism=True then do not count as a piece a subword of a relator that occurs only periodically.

    >>> pieces(['aa'],asstring=True,piece_up_to_automorphism=True)
    set()
    >>> sorted(pieces(['aa'],asstring=True))
    ['A', 'AA', 'a', 'aa']
    >>> sorted(pieces(['abABcdCD'],asstring=True))
    ['A', 'B', 'C', 'D', 'a', 'b', 'c', 'd']
    >>> sorted(pieces(['abbccbABCCB'],asstring=True))
    ['A', 'B', 'BC', 'BCC', 'BCCB', 'C', 'CB', 'CC', 'CCB', 'a', 'b', 'bc', 'bcc', 'bccb', 'c', 'cb', 'cc', 'ccb']
    >>> sorted(pieces([[1,2,3,1,2],[2,2,1,1,2,3]]))
    [[-3], [-3, -2], [-3, -2, -1], [-2], [-2, -1], [-1], [-1, -2], [1], [1, 2], [1, 2, 3], [2], [2, 1], [2, 3], [3]]
    """
    if noparse:
        rels=relatorlist
    else:
        rels=parseinputwords(relatorlist,asrelatorlist=True)
    pieces=set()
    irels=[intlisttostring(rel) for rel in itertools.chain.from_iterable(zip(rels,[inverse(w) for w in rels]))] # arrange relators and inverses in a list of the form relator1, inverse of relator1, relator2, inverse of relator2,...
    drels=[x+x for x in irels]
    for relatorindex in range(len(rels)): # only need to search relators for candidate pieces, since a piece contained in inverse will be inverse of piece contained in relator
        relator=irels[2*relatorindex]
        for L in range(1,1+len(relator)): # L is length of prospective piece
            for startingindex in range(len(relator)):
                p=(relator+relator)[startingindex:startingindex+L] # the subword of length L starting at index i in reltaor as a cyclic word
                # now we need to check if p is a piece
                # we do not need to check lower relatorindices, because we already scanned those relators for pieces
                # first check higher indices
                if any(p in x for x in [drels[i] for i in range(2*relatorindex+1,len(drels))]):# found a matching subword, p is a piece
                    pieces.add(p)
                    pieces.add(''.join(reversed(p.swapcase())))
                    continue
                # didn't find any occurences in p in higher index relators. check in different places in relator
                if not piece_up_to_automorphism:
                    if any(p in x for x in [(relator+relator)[startingindex+1:len(relator)+startingindex+L-1]]):
                        pieces.add(p)
                        pieces.add(''.join(reversed(p.swapcase())))
                        continue
                else:
                    the_root,the_power=maxroot(relator)
                    if any(p in x for x in [(relator+relator)[startingindex+1:len(the_root)+startingindex+L-1]]):
                        pieces.add(p)
                        pieces.add(''.join(reversed(p.swapcase())))
                        continue
    if asstring:
        return pieces
    else:
        return [stringtointlist(pieces) for pieces in pieces]



def parseinputwords(inputwords,asrelatorlist=True):
    """
    Take as input two different possible representations of lists of words in a free group and return the same represented as lists of lists of nonzero integers where positive i indicates the i-th generator and -i indicates its inverse. 

    If input is already in the form of list of lists of integers then it is returned. If input is a list of alphabetic strings then it is convertex to list of lists of integers with a -> 1, A -> -1, b -> 2, etc.

    Raise an exception if input is not one of two types described above.

    If asrelatorlist=True, raise an exception if some word in the list is not freely or cyclically reduced, or if some word in the list has length less than 2, or if some pair of words are conjugates or inverse conjugates.  
    """
    if all(type(w)==str and (w=='' or w.isalpha()) for w in inputwords):
        rels=[stringtointlist(w) for w in inputwords]
    elif all(all(type(x)==int and x!=0 for x in w) for w in inputwords):
        rels=inputwords
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
    []
    >>> freelyreduce([1,2,3])
    [1, 2, 3]
    >>> freelyreduce([1,2,1,-1,2])
    [1, 2, 2]
    >>> freelyreduce([1,2,-2,3,1,2,-2,-1,-3,-1])
    []
    """
    reduced=intlist.copy()
    if len(reduced)<2:
        return reduced
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
    return reduced
        

def cyclicreduce(intlist):
    """
    >>> cyclicreduce([])
    []
    >>> cyclicreduce([1,-1])
    []
    >>> cyclicreduce([2,1,-1,2,3,-2,-2])
    [3]
    >>> cyclicreduce([1,2,1,2,-1,-2,-2,-1])
    [1, 2, -1, -2]
    """
    theword=freelyreduce(intlist)
    if len(theword)<2:
        return theword
    conjugatorindex=0
    while conjugatorindex<=len(theword)//2 and theword[conjugatorindex]==-theword[-conjugatorindex-1]:
        conjugatorindex+=1
    return theword[conjugatorindex:len(theword)-conjugatorindex]

def inverse(input_word):
    """
    Return the inverse of the given input_word in the same form as given.

    >>> inverse([])
    []
    >>> inverse([1])
    [-1]
    >>> inverse([-2])
    [2]
    >>> inverse([1,2,-3])
    [3, -2, -1]
    >>> inverse('')
    ''
    >>> inverse('a')
    'A'
    >>> inverse('B')
    'b'
    >>> inverse('abC')
    'cBA'
    """
    if type(input_word)==list and all(type(x)==int for x in input_word):
        return list(map(lambda x:-1*x, input_word[::-1]))
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
    
    
                
