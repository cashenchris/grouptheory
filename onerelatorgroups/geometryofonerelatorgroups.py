import networkx as nx
import grouptheory.automaticgroups as ag
import grouptheory.smallcancellation as sc
import grouptheory.freegroups.freegroup as fg
import grouptheory.freegroups.whiteheadgraph as wg
from grouptheory.onerelatorgroups.IvanovSchupp import IvanovSchupp
from grouptheory.onerelatorgroups.BlufsteinMinian import BlufsteinMinian
import pexpect
from fractions import Fraction
import itertools
import os
from collections import deque
import grouptheory.freegroups.enumeratefreegroupwords as enum
import glob
#import subprocess32
import subprocess

# One of the checks that certify_hyperbolicity does is to use the walrus package in GAP. The GAP startup takes a while. If multiple checks are to be run it is faster to spawn GAP only once and reuse it as follows:
# words = list of words in free group of rank r, each defining a 1=relator quotient
# >>> gapinstance=spawngapforwalrus(r)
# >>> for word in words:
#             certify_hyperbolicity(word, gap=gapinstance)
def is_hyperbolic(relator,reportreason=False,**kwargs):
    """
    Check if one-relator group defined by input relator is hyperbolic.
    Returns True if group is known hyperbolic, False if known not hyperbolic, None if inconclusive.

    If reportreason=True output additionally contains string with name of the method certifying hyperbolicity, or None is all tests are inconclusive.
    tryhard=2 will make kbmag retry several times with random orderings of the generators if the first attempt fails.

    pass no_minimization=True to check the relator as given. Otherwise the relator will first be Whitehead minimized.


    >>> is_hyperbolic([1,2,3],reportreason=True)
    (True, 'free')
    >>> is_hyperbolic('bbb',reportreason=True)
    (True, 'torsion')
    >>> is_hyperbolic([1,2,-1,-2,3,4,-3,-4],reportreason=True)
    (True, 'cyclically pinched')
    >>> is_hyperbolic('abcabcdeDEdeDE',reportreason=True)
    (False, 'cyclically pinched')
    >>> is_hyperbolic('ababcabccabcbcbcbcabcbcbcbcbc',reportreason=True)
    (True, 'Ivanov Schupp')
    >>> is_hyperbolic('cacbcbcbcabacbcaba',reportreason=True,no_minimization=True) # this example is small cancellation as written but not after minimzation 
    (True, 'small cancellation')
    >>> is_hyperbolic('CCBBCAAbbcaa',reportreason=True)
    (True, 'walrus')
    """
    def format_return(result,reason):
        if reportreason:
            return result,reason
        else:
            return result
    if not relator:
        return format_return(True,'free')
    F,r1=fg.parseinputword(relator)
    if 'no_minimization' in kwargs and kwargs['no_minimization']: # do not replace word with Whitehead minimal word in same automorphic orbit
        r2=r1
    else: # find a Whitehead minimal element in the same automorphism orbit of the relator
        r2=wg.whitehead_minimal_representative(r1)
    if len(r2)<=1:
        return format_return(True,'free')
    if F.degree(r2)>1:# one-relator groups with torsion are hyperbolic
        return format_return(True,'torsion')
    # Check if relator is cyclically pinched
    (m,n)=is_cyclically_pinched(r2,reportpowers=True)
    if m is not None:
        if m>1 and n>1:
            return format_return(False,'cyclically pinched')
        else:
            return format_return(True, 'cyclically pinched')
    # Try Ivanov-Scupp criteria, always gives an answer if some generator appears at most 3 times
    if 'verbose' in kwargs and kwargs["verbose"]:
        print("Checking Ivanov-Schupp criteria.")
    ISresult=IvanovSchupp(r2,reportreason=False)
    if ISresult is not None:
        return format_return(ISresult,'Ivanov Schupp')
    if 'verbose' in kwargs and kwargs["verbose"]:
        print("Ivanonv-Schupp does not apply.")
    # check if relator defines small cancellation presentation
    Cprimebound=sc.Cprimebound([r2])
    if sc.smallcancellation([r2],Cprimebound):
        return format_return(True,'small cancellation')
    # check Blufstein-Minian condition
    if BlufsteinMinian(r2,Cprimebound):
        return format_return(True,'Blufstein Minian')
    if 'verbose' in kwargs and kwargs["verbose"]:
        print("Not small cancellation.")
    # try GAP with walrus
    if 'walrus' in kwargs and kwargs['walrus']==False:
        pass
    else:
        if 'verbose' in kwargs and kwargs['verbose']:
            print("Trying GAP with walrus.")
        if 'gapprompt' in kwargs:
            gapprompt=kwargs["gapprompt"]
        else:
            gapprompt='gap>'
        if 'gapfreegroupname' in kwargs:
            gapfreegroupname=kwargs["gapfreegroupname"]
        else:
            gapfreegroupname='f'
        if 'walrusparameter' in kwargs:
            walrusparameter=kwargs["walrusparamter"]
        else:
            walrusparameter='1/100'
        if 'gap' in kwargs:
            walrus=checkhyperbolicitywithwalrus(r2.letters,gap=kwargs["gap"],gapfreegroupname=gapfreegroupname,gapprompt=gapprompt,theparameter=walrusparameter)
        else:
            if 'pathtogap' in kwargs:
                walrus=checkhyperbolicitywithwalrus(r2.letters,gap=None,gapfreegroupname=gapfreegroupname,gapprompt=gapprompt,theparameter=walrusparameter,pathtogap=kwargs['pathtogap'])
            else:
                walrus=checkhyperbolicitywithwalrus(r2.letters,gap=None,gapfreegroupname=gapfreegroupname,gapprompt=gapprompt,theparameter=walrusparameter)
        if walrus:
            return format_return(True,'walrus')
        if 'verbose' in kwargs and kwargs["verbose"]:
            print("walrus failed.")
    # try kbmag
    if 'kb' in kwargs and kwargs['kb']==False:
        return format_return(None,None)
    if 'verbose' in kwargs and kwargs["verbose"]:
        print("Trying kbmag")
    autandhyp=ag.certify_hyperbolicity(r2(),**kwargs)
    if autandhyp:
        return format_return(True,'kbmag')
    else: # all tests inconclusive
        return format_return(None, None)


def is_cyclically_pinched(relator,reportwords=False,reportpowers=False):
    F,rels=fg.parseinputwords([relator])
    r=rels[0].letters
    rr=r+r
    for (startingindex,L) in itertools.product(range(len(r)),range(1,len(r))):
        prefix=rr[startingindex:startingindex+L]
        prefixletters={abs(x) for x in prefix}
        suffix=rr[startingindex+L:startingindex+len(r)]
        suffixletters={abs(x) for x in suffix}
        if not prefixletters&suffixletters:
            result=True
            words=(F.word(prefix),F.word(suffix))
            if reportpowers:
                powers=(F.degree(words[0]),F.degree(words[1]))
            break
    else:
        result=False
        words=(None,None)
        powers=(None,None)
    if reportwords:
        return words
    if reportpowers:
        return powers
    return result
            
            



def checkhyperbolicitywithwalrus(theword,theparameter='1/100',gap=None,gapfreegroupname='f',gapprompt='gap>',fulloutput=False,pathtogap=None):
    """
    Check hyperbolicity of the one relator group with relator defined by theword using the walrus package in GAP.

    Input 'gap' should be a pexpect process of GAP with walrus loaded, a free group defined, ready for input. If None, process will be spawned. Spawning takes time, so if this function will be run more than once is is better to spawn one gap process and reuse it for each function call. 
    """ 
    if gap is None:
        gap=spawngapforwalrus(max(abs(x) for x in theword),gapfreegroupname=gapfreegroupname,gapprompt=gapprompt,pathtogap=pathtogap)
    walrusstring='IsHyperbolic(PregroupPresentationFromFp('+gapfreegroupname+',[],['+converttogapword(theword,gapfreegroupname)+']),'+theparameter+');'
    gap.sendline(walrusstring)
    gap.expect(gapprompt)
    output=gap.before
    if 'true' in output:
        result=True
    elif 'fail' in output:
        result=False
    else:
        raise ValueError("Did not find 'true' or 'fail' in output:",output)
    if fulloutput:
        return result,output
    else:
        return result

def spawngapforwalrus(rank,gapfreegroupname='f',gapprompt='gap>',pathtogap=None):
    """
    Sapwns a GAP process, loads package walrus, and defines free group of specified rank.
    """
    # Sometimes pexpect can find gap on its own. If not, specify the full path to gap.
    if pathtogap is None:
        gap=pexpect.spawn('gap -b', encoding='utf-8')
    else:
        gap=pexpect.spawn(pathtogap+' -b',encoding='utf-8')
    gap.expect(gapprompt)
    gap.sendline('walrusloaded:=LoadPackage("walrus");')
    gap.expect(gapprompt)
    gap.sendline('walrusloaded;')
    gap.expect(gapprompt)
    output=gap.before
    if not 'true' in output:
        raise ValueError("Unable to load package 'walrus' in this version of gap.")
    gap.sendline(gapfreegroupname+':=FreeGroup('+str(rank)+');;')
    gap.expect(gapprompt)
    return gap

def converttogapword(listsofints,gapfreegroupname='f'):
    """
    Convert list of nonzero intergers representing a word in a free group to a string readable by GAP.

    >>> converttogapword([1,2,-1,-2,])
    'f.1*f.2*f.1^-1*f.2^-1'
    """
    outputstring=list()
    for i in range(len(listsofints)):
        x=listsofints[i]
        if x>0:
            outputstring.append(gapfreegroupname+'.'+str(x))
        else:
            outputstring.append(gapfreegroupname+'.'+str(abs(x))+"^-1")
        if i<len(listsofints)-1:
            outputstring.append("*")
    return ''.join(outputstring)



#----------------------------------

def girth(relator,verbose=False,**kwargs):
    """
    Given a word defining a shortlex automatic one-relator group, return the girth of the Cayley graph for the corresponding one-relator presentation and a pair of equivalent words witnessing shortest relation.

    Default generator ordering is [...,'B','A','a','b',...]. Use keyword 'generators' to specify different ordering of generators and inverses.
    Use filename=kbmagbasefilename to specify existing kbmag files defining the automatric structure.
    Use cleanup=True to delete all auxilliary kbmag files or cleanup=False to keep them. Default is to cleanup if this function creates the files and not if they already exist.
    Set timeout=n to set time limit of n seconds on attempt to generate automatric group structure.
    """
    if not relator:
        return float('inf')
    F,r=fg.parseinputword(relator)
    relatorasstring=r()
    relatoraslist=r.letters
    if 'filename' in kwargs:
        thefilename=kwargs['filename']
    else:
        thefilename="OneRelatorGroup-"+relatorasstring
    if 'cleanup' in kwargs:
        cleanup=kwargs['cleanup']
    else:
        cleanup=not os.path.isfile(thefilename) # default is to cleanup if the files are not already existing
    if not os.path.isfile(thefilename):
        if 'generators' in kwargs:
            generators=kwargs['generators']
        else:
            generators=[x.upper() for x in F.gens[::-1]]+F.gens
        ag.writetokbmagfile(thefilename,generators,[relatorasstring])
    if 'timeout' in kwargs:
        timeout=kwargs['timeout']
    else:
        timeout=10 # default timeout 10s       
    if not os.path.isfile(thefilename+'.diff1'):
        try:
            #subprocess32.run(['autgroup','-silent',thefilename],check=True,timeout=timeout)
            subprocess.run(['autgroup','-silent',thefilename],check=True,timeout=timeout)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
            #except (subprocess32.TimeoutExpired, subprocess32.CalledProcessError) as e:
            if cleanup:
                files = glob.glob('./'+thefilename+"*")
                for file in files:
                    os.remove(file)
            raise e
    foundrelation=False
    currentlength=0
    while not foundrelation:
        currentlength+=1
        if 2*currentlength==len(relatorasstring) or 2*currentlength-1==len(relatorasstring):
            return len(relatorasstring),relatorasstring,''
        if verbose:
            print("Searching words of length "+str(currentlength))
        wordsofcurrentlength=enum.generate_words(F.rank,currentlength,currentlength)
        for theintlist in wordsofcurrentlength:
            thestring=F.word(theintlist)()
            reducedstring=ag.groupelement(thestring,thefilename).string
            if len(reducedstring)<len(thestring):
                assert(len(reducedstring)+1==len(thestring))
                if 'cleanup' in kwargs and kwargs['cleanup']==True:
                    files = glob.glob('./'+thefilename+"*")
                    for file in files:
                        try:
                            os.remove(file)
                        except:
                            pass
                return 2*currentlength-1,thestring,reducedstring
            if thestring!=reducedstring: # found a relation of lengh currentlength*2, but keep searching for one of length 2*currentlength-1
                foundrelation=True
                u=thestring
                v=reducedstring
    if 'cleanup' in kwargs and kwargs['cleanup']==True:
        files = glob.glob('./'+thefilename+"*")
        for file in files:
            try:
                os.remove(file)
            except:
                pass
    return currentlength*2,u,v
    

def ball(F,relator,thefilename,verbose=False,**kwargs):
    NotImplemented
    try:
        d=deque(relator.letters)
    except AttributeError:
        d=deque(relator)
    cyclicpermutations=[]
    for i in range(len(d)):
        cyclicpermutations.append(list(d))
        d.rotate()
    B=nx.MultiDiGraph()
    B.add_node((0,0))
    for i in range(len(relator)):
        if cyclicpermutations[i][0]>0:
            B.add_edge((0,0),(i,1),label=cyclicpermutations[i][0])
        else:
            B.add_edge((i,1),(0,0),label=-cyclicpermutations[i][0])
        if cyclicpermutations[i][-1]>0:
            B.add_edge((i,len(relator)-1),(0,0),label=cyclicpermutations[i][-1])
        else:
            B.add_edge((0,0),(i,len(relator)-1),label=-cyclicpermutations[i][-1])
        for j in range(1,len(relator)-1):
            if cyclicpermutations[i][j]>0:
                B.add_edge((i,j),(i,j+1),label=cyclicpermutations[i][j])
            else:
                B.add_edge((i,j+1),(i,j),label=-cyclicpermutations[i][j])
    sf.fold(B,(0,0))
    # now use automatric structure to check that all vertices are really distinct
    T=nx.maximum_spanning_tree(nx.MultiGraph(B))
    for n in B:
        theword=F.word(edgeword(B,T,(0,0),n))
        minword=ag.groupelement(theword(),thefilename).string
        B.nodes[n]['representative']=minword
    del T
    originalnodelist=[n for n in B.nodes() if n!=(0,0)]
    for n in originalnodelist:
        if n not in B:
            continue
        for m in [x for x in B if x!=n]:
            if B.nodes[n]['representative']==B.nodes[m]['representative']:
                sf.wedge(B,n,m)
                break
    sf.fold(B)
    return B

def edgeword(B,T,v,w):
    path=nx.shortest_path(T,v,w)
    theword=[]
    for i in range(len(path)-1):
        for e in T.edges(path[i],keys=True,data=True):
            if e[1]==path[i+1]:
                break
        else:
            raise RuntimeError
        if e in B.edges(data=True,keys=True):
            theword.append(e[3]['label'])
        elif (e[1],e[0],e[2],e[3]) in B.edges(data=True,keys=True):
            theword.append(-e[3]['label'])
        else:
            raise RuntimeError
    return theword
            













# in terminal, do
# python geometryofonerelatorgroups.py
# to run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()
