import freegroups.freegroup as freegroup
import freegroups.AutF as aut
import grouptheory.group as group
import copy
import freegroups.whiteheadgraph as wg
from knownexamples import *



def rjsjtest(maxlength,verbose,debug,randomautomorphismlength,examplename,freegroup,wordlist,splitsfreely):
    nonefailed=True
    if splitsfreely:
        return nonfailed
    # take a known example and mix it up with an automorphism alpha
    F=freegroup
    rank=F.rank
    alpha,alphainv=aut.random_automorphism_pair(F,randomautomorphismlength)

    if verbose:
        print("Trying example ", examplename, " changed by automorphism:\n", alpha)
    newwordlist=[alpha(w) for w in wordlist]
    gamma, wordmap=F.get_RJSJ(newwordlist, withmap=True)
    minimizedwordlist=[]
    for (v,w,p) in wordmap:
        minimizedwordlist.append(gamma.localgroup(v).get_inclusion(F)(w))
    if not F.is_RJSJ(wordmap,gamma,verbose=verbose):
        if verbose:
            print("Error computing RJSJ for", examplename,".")
        nonefailed=False
        if debug:
            print(gamma, wordmap)
            print('********************************')

    if verbose and nonefailed:
        print("Correctly found RJSJ for",examplename,".")
    return nonefailed
        
def testall(maxlength=30, randomautomorphismlength=0,verbose=False, debug=False):
    if all([rjsjtest(maxlength,verbose,debug,randomautomorphismlength,k,knownexamples[k]['freegroup'],knownexamples[k]['wordlist'],knownexamples[k]['splitsfreely']) for k in knownexamples]):
        print("Found expected rJSJ's.")
    else:
        print("Some rJSJ test failed.")
