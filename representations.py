import itertools
import permutationgroups
import group

# G should be a group.FPGroup and T should be a list of permutationgroups.permutation of some set.


def is_representation(G,T):
    """
    Return True if the mapping sending the i-th generator of the finitely presented group G to the i-th entry in the list of permutations T defines a homomorphism to the group of permutations of the set that T acts on. 
    """
    Tdict=dict()
    for i in range(len(T)):
        Tdict[i+1]=T[i]
        Tdict[-i-1]=T[i].inverse()
    return all((reduce(lambda x, y: x*y, [Tdict[z] for z in r.letters])).is_trivial() for r in G.rels)
    
def generate_representation(G,permutation_generator):
    """
    Yields a list of permutations of n elements such that the map sending the i-th generator of the finitely presented group G to the i-th element of the list defines a homomoprhism into Sym(n).
    """
    for T in itertools.product(permutation_generator,repeat=len(G.gens)):
        if is_representation(G,T):
            yield T

def group_generated_by(T):
    """
    Return as a set the subgroup of the symmetric groups generated by the permutations in T.
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