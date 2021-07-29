from math import gcd
import numpy as np
import itertools
import grouptheory.group

# class permutation has main representation as a dictionary of input:ouput pairs. 

class permutation(object):
    """
    Defines a permutation. Initialization by list of tuples representing disjoint cycles or by dictionary containing input:output pairs. Argument 'chirality' can be 'L' or 'R' according to whether permutation action is left action or right action. Left action means permutation acts on the left of set that is being permuted, so in a product of permutations the right most factor acts first. Default chirality is left.

    >>> p=permutation([(1,2,3),(4,5)])
    >>> q=permutation([(3,4)])
    >>> p(1)
    2
    >>> p(6)
    6
    >>> p**(-1)
    [1, 3, 2][4, 5]
    >>> (p*q)(3)
    5
    >>> p=permutation([(1,2,3),(4,5)],chirality='R')
    >>> q=permutation([(3,4)],chirality='R')
    >>> (p*q)(3)
    1
    """
    def __init__(self,thecycles=None,thedict=None,chirality='L'):
        if chirality=='L' or chirality=='l':
            self.chirality='L'
        elif chirality=='R' or chirality=='r':
            self.chirality='R'
        else:
            raise InputError
        if thedict:
            self.thedict=dict({k:v for k, v in thedict.items() if k!=v})
        elif thecycles:
            thisdict=dict()
            for cycle in thecycles:
                if len(cycle)>1:
                    for i in range(len(cycle)):
                        thisdict[cycle[i]]=cycle[(i+1)%len(cycle)]
            self.thedict=dict(thisdict)
        else:
            self.thedict=dict()
        self.thefrozenset=frozenset({(k,self.thedict[k]) for k in self.thedict})
        thecycletypes=dict()
        for cycle in self.cycles():
            try:
                thecycletypes[len(cycle)]+=1
            except KeyError:
                thecycletypes[len(cycle)]=1
        self.thecycletypes=dict(thecycletypes)
        self.theorder=lcm(*[l for l in self.thecycletypes.keys()])
        self.thehash=hash(self.thefrozenset)
        
    def __repr__(self):
        thestring=''
        for cycle in self.cycles():
                thestring+=str(cycle)
        if thestring=='':
            return '[]'
        else:
            return thestring

    def __call__(self,arg):
        try:
            return self.thedict[arg]
        except KeyError:
            return arg
    
    def __mul__(self,other):
        assert(self.chirality==other.chirality)
        if self.chirality=='L': # 'other' acts first
            secondkeys=[k for k in self.thedict.keys() if k not in other.thedict.keys()]
            firstkeys=[k for k in other.thedict.keys()]
            newdict=dict()
            for k in firstkeys:
                newdict[k]=self(other(k))
            for k in secondkeys:
                newdict[k]=self(k)
        else: # self.chirality=='R' means 'self' acts first
            secondkeys=[k for k in other.thedict.keys() if k not in self.thedict.keys()]
            firstkeys=[k for k in self.thedict.keys()]
            newdict=dict()
            for k in firstkeys:
                newdict[k]=other(self(k))
            for k in secondkeys:
                newdict[k]=other(k)
        return permutation(None,newdict,chirality=self.chirality)
    
    def inverse(self):
        return permutation(None,{v: k for k, v in self.thedict.items()},chirality=self.chirality)

    def order(self):
        return self.theorder

    def __pow__(self,exponent):
        thepow=permutation(chirality=self.chirality) # the trivial permutation
        for i in range(exponent%self.theorder):
            thepow=thepow*self
        return thepow

    def is_trivial(self):
        return self.theorder==1

    def is_conjugate(self,other):
        if self.theorder!=other.theorder:
            return False
        else:
            return self.thecycletypes==other.thecycletypes

    def __eq__(self,other):
        if self.is_conjugate(other):
            return self.thedict==other.thedict
        else:
            return False

    def __hash__(self):
        return self.thehash

    def __ne__(self,other):
        return not self == other

    def cycles(self):
        elements=self.support()
        cycles=[]
        while elements:
            firstelement=elements.pop()
            thiscycle=[firstelement]
            currentelement=firstelement
            nextelement=self.thedict[currentelement]
            while nextelement!=firstelement:
                thiscycle.append(nextelement)
                currentelement=nextelement
                elements.remove(currentelement)
                nextelement=self.thedict[currentelement]
            cycles.append(thiscycle)
        return cycles
    
    def support(self):
        return set(self.thedict.keys())

    def sign(self):
        return sum(self.thecycletypes[k]*(k-1) for k in self.thecycletypes)%2
        

def symmetric_group_gen(degree,chirality='L'):
    """
    Generator that yields elements of the symmetric group of specified degree.
    """
    P=itertools.permutations(range(1,1+degree))
    for p in P:
        d={i+1:p[i] for i in range(degree)}
        yield permutation(thedict=d,chirality=chirality)

def alternating_group_gen(degree,chirality='L'):
    """
    Generator that yields elements of the symmetric group of specified degree.
    """
    P=itertools.permutations(range(1,1+degree))
    for p in P:
        d={i+1:p[i] for i in range(degree)}
        p=permutation(thedict=d,chirality=chirality)
        if p.sign()==0:
            yield p



def generate_permutation_representations(G,n):
    """
    Yields lists of permutations of n elements such that the map sending the i-th generator of the finitely presented group G to the i-th element of the list defines a homomoprhism into Sym(n).
    """
    for T in itertools.product(symmetric_group_gen(n),repeat=len(G.gens)):
        if isrepresentation(G,T):
            yield T

def generate_alternating_representation(G,n):
    """
    Yields lists of permutations of n elements such that the map sending the i-th generator of the finitely presented group G to the i-th element of the list defines a homomoprhism into Alt(n).
    """
    for T in itertools.product(symmetric_group_gen(n),repeat=len(G.gens)):
        if sum(x.sign() for x in T)==0 and isrepresentation(G,T):
            yield T

def is_permutation_representation(G,T):
    """
    G is a finitely presented group. T is a list of permutations of some set.
    Return True if the mapping sending the i-th generator of the finitely presented group G to the i-th entry in the list defines a homomorphism. 
    """
    Tdict=dict()
    for i in range(len(T)):
        Tdict[i+1]=T[i]
        Tdict[-i-1]=T[i].inverse()
    return all((reduce(lambda x, y: x*y, [Tdict[z] for z in r.letters])).is_trivial() for r in G.rels)
    
def group_generated_by(T):
    """
    Return as a set the subgroup of the symmetric group generated by the permutations in T.
    """
    if not T:
        return set([])
    thegroup=set([T[0]*(T[0]).inverse()]) # the set consisting of the trivial element in the group
    gens=set(T)|set([t.inverse() for t in T])
    gens-=thegroup #discard trivial element from gen set
    new=set(T)
    thegroup|=new
    while new:
        x=new.pop()
        for y in gens:
            z=x*y
            if z not in thegroup:
                new.add(z)
                thegroup.add(z)
    return thegroup

    

def action(input_tuple,thepermutation):
    return tuple(thepermutation(i) for i in input_tuple)

def lcm(*numbers):
    """Return least common multiple."""    
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    if not numbers:
        return 1
    else:
        thelcm=1
        for n in numbers:
            thelcm=lcm(thelcm,n)
        return thelcm

if __name__ == "__main__":
    import doctest
    doctest.testmod()
