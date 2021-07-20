import freegroups.whiteheadgraph as wg
import copy
import grouptheory.group as group
#import freegroups.whiteheadgraph.orderedmultigraph as omg
import freegroups.AutF as aut
import freegroups.freegroup as freegroup
from knownexamples import *





def cutpairtest(maxlength,verbose,debug,randomautomorphismlength,examplename,freegroup,wordlist,splitsfreely,iscircle,isrigid,cutpoints,uncrossed):
    nonefailed=True
    # take a known example and mix it up with an automorphism alpha
    F=freegroup
    rank=F.rank
    alpha,alphainv=aut.random_automorphism_pair(F,randomautomorphismlength)

    if verbose:
        print("Trying example ", examplename, " changed by automorphism:\n", alpha)
    
    wm=wg.whitehead_minimal(F,[alpha(w) for w in wordlist], verbose=verbose)
    minimizingautomorphism=wm['minimizingautomorphism']
    newwordlist=wm['wordlist']
    W=wg.WGraph(newwordlist, simplified=True, verbose=verbose)
    wholeautomorphism=minimizingautomorphism*alpha
    
    newcutpoints=[wholeautomorphism(cutpoint) for cutpoint in cutpoints] 
    newuncrossed=[wholeautomorphism(uncross) for uncross in uncrossed]
    
    if not F.are_equivalent_wordlists(W.get_wordlist(),newwordlist):
        if verbose:
            print("Error in get_wordlist for ", examplename)
        nonefailed=False
    if  splitsfreely!=F.splits_freely_rel(W.wordlist):
        if verbose:
            print("Error in splits_freely for ", examplename)
        nonefailed=False
    if  iscircle!=W.is_circle():
        if verbose:
            print("Error in is_circle for ", examplename)
        nonefailed=False
    if isrigid!=F.is_rigid_rel(W,maxlength):
        if verbose:
            print("Error in is_rigid_rel for ", examplename)
        nonefailed=False
    if not F.are_equivalent_wordlists(newcutpoints,wg.find_cut_points(F,W)):
        if verbose:
            print("Error in split.find_cut_points for ", examplename)
        nonefailed=False
    #cuts=split.find_cut_pairs(F,W,newwordlist,maxlength)[0]
    cuts=wg.find_universal_splitting_words(F,W,newwordlist,maxlength)[0]
    if not F.are_equivalent_wordlists(list(cuts['cutpoints']),newcutpoints):
        if verbose:
            print("Error finding cut points in split.findCutPairs for ", examplename)
        nonefailed=False
    if not F.simplify_wordlist(list(cuts['uncrossed']),newuncrossed)==[]:
        if verbose:
            print("Error too many uncrossed cut pairs in split.findCutPairs for ", examplename)
        nonefailed=False
    if not F.simplify_wordlist(newuncrossed,list(set.union(set(cuts['uncrossed']),set(cuts['othercuts']))))==[]:
        if verbose:
            print("Error didn't find all uncrossed cut pairs in split.findCutPairs for ", examplename)
        nonefailed=False        

    # test some random words to see if splitting info is correct
    for i in range(1,maxlength):
        w=F.randomword(i)
        w=F.cyclic_reduce(w)
        if len(w)>0:
            if iscircle:
                if not wg.gives_cut(F,W,w)!=bool(F.is_conjugate_into(w,*newwordlist)):
                    if verbose:
                        print("Error: W is a circle, so ",w," should be a cut pair in ", examplename)
                    nonefailed=False
                    break
            else:
                if not wg.gives_cut(F,W,w)==bool(F.is_conjugate_into(w,*set.union(set(newuncrossed),set(newcutpoints)))):
                    if verbose:
                        print("Warning",w()," gives a cut but wasn't found in ", examplename)
                        print("It may be that ",w()," is a crossed cut pair and everything is ok. Check by hand.")
                    #nonefailed=False
                    break
    if verbose and nonefailed:
        print("All tests passed for ",examplename,".")
    return nonefailed
        
def testall(maxlength=30, randomautomorphismlength=0,verbose=False,debug=False):
    
    if all([cutpairtest(maxlength,verbose,debug,randomautomorphismlength,k,knownexamples[k]['freegroup'],knownexamples[k]['wordlist'],knownexamples[k]['splitsfreely'],knownexamples[k]['iscircle'],knownexamples[k]['isrigid'],knownexamples[k]['cutpoints'],knownexamples[k]['uncrossed']) for k in knownexamples]):
        print("Found all expected cut points/pairs")
    else:
        print("Failed to find expected cut points/pairs")
    


