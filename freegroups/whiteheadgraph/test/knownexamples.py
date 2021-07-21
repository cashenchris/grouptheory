import grouptheory.group as group
import grouptheory.freegroups.freegroup as freegroup
import grouptheory.freegroups.whiteheadgraph as wg

F2=freegroup.FGFreeGroup(numgens=2)
F3=freegroup.FGFreeGroup(numgens=3)
# a few words
com=F2.word([1,2,-1,-2])
a=F2.word([1])
b=F2.word([2])
bs12=F2.word([1,2,-1,-1,-2])
bs13=F2.word([1,2,-1,-1,-1,-2])
tk=[com,F2.word([1,2,1,-2])]
bs12a=[bs12,a]
k6=[F3.word([1]),F3.word([2]),F3.word([3]),F3.word([1,2]),F3.word([1,3]),F3.word([2,3]),F3.word([1,-2]),F3.word([1,-3]),F3.word([2,-3])]
baumslag=F2.word([-1,-1,-2,-1,2,1,-2,1,2])
k33=F3.word([2,2,1,1,3,3,1,2,3])
bs23=F2.word([1,1,2,-1,-1,-1,-2])


# their corresponding Whitehead graphs
Com=wg.WGraph([com])
#A=wg.WGraph(a)
#B=wg.WGraph(b)
BS12=wg.WGraph([bs12])
BS13=wg.WGraph([bs13])
TK=wg.WGraph(tk)
K4=wg.WGraph([a,b,com]) # rigid
K4HNN=wg.WGraph([a,F2.word([1,-2,1,2,-1,-2,-1,2])])
BS12A=wg.WGraph([bs12,a])
K6=wg.WGraph(k6) # rigid
AS=wg.WGraph([a,b,F2.word([1,2,1,-2,-1,2,-1,-2])]) # rigid 'almost splits'
NonTree=wg.WGraph([F3.word([-1,2,3]),F3.word([-1,3,2]),F3.word([-1,2,2,2]),F3.word([-1,3,3,3])]) # rigid but rigid cube complex is not a tree
K33=wg.WGraph([k33])
Baumslag=wg.WGraph([baumslag])
BS23=wg.WGraph([bs23])

knownexamples={
    #'Name':{freegroup,whiteheadgraph,wordlist,splitsfreely,iscircle,isrigid,cutpoints,uncrossed, virtuallygeometric}
'BS23':{'freegroup':F2,'whiteheadgraph':BS23,'wordlist':[bs23],'splitsfreely':False,'iscircle':False,'isrigid':False, 'cutpoints':[],'uncrossed':[a], 'virtuallygeometric':True},
'Com':{'freegroup':F2,'whiteheadgraph':Com,'wordlist':[com],'splitsfreely':False,'iscircle':True,'isrigid':False,'cutpoints':[],'uncrossed':[], 'virtuallygeometric':True},
'TK':{'freegroup':F2,'whiteheadgraph':TK,'wordlist':tk,'splitsfreely':False,'iscircle':False,'isrigid':False,'cutpoints':[],'uncrossed':[a], 'virtuallygeometric':None},
'K4':{'freegroup':F2,'whiteheadgraph':K4,'wordlist':[a,b,com],'splitsfreely':False,'iscircle':False,'isrigid':True,'cutpoints':[],'uncrossed':[], 'virtuallygeometric':True},
'K4HNN':{'freegroup':F2,'whiteheadgraph':K4HNN,'wordlist':[a,F2.word('aBabABAb')],'splitsfreely':False,'iscircle':False,'isrigid':False,'cutpoints':[a],'uncrossed':[], 'virtuallygeometric':True},
'BS12':{'freegroup':F2,'whiteheadgraph':BS12,'wordlist':[bs12],'splitsfreely':False,'iscircle':False,'isrigid':False,'cutpoints':[],'uncrossed':[a], 'virtuallygeometric':True},
'BS12A':{'freegroup':F2,'whiteheadgraph':BS12A,'wordlist':bs12a,'splitsfreely':False,'iscircle':False,'isrigid':False,'cutpoints':[a],'uncrossed':[], 'virtuallygeometric':True},
'K6':{'freegroup':F3,'whiteheadgraph':K6,'wordlist':k6,'splitsfreely':False,'iscircle':False,'isrigid':True,'cutpoints':[],'uncrossed':[], 'virtuallygeometric':False},
'AS':{'freegroup':F2,'whiteheadgraph':AS,'wordlist':[a,b,F2.word([1,2,1,-2,-1,2,-1,-2],)],'splitsfreely':False,'iscircle':False,'isrigid':True,'cutpoints':[],'uncrossed':[], 'virtuallygeometric':None},
'NonTree':{'freegroup':F3,'whiteheadgraph':NonTree,'wordlist':[F3.word([-1,2,3]),F3.word([-1,3,2]),F3.word([-1,2,2,2]),F3.word([-1,3,3,3])],'splitsfreely':False,'iscircle':False,'isrigid':True,'cutpoints':[],'uncrossed':[], 'virtuallygeometric':None},
'K33':{'freegroup':F3,'whiteheadgraph':K33,'wordlist':[k33],'splitsfreely':False,'iscircle':False,'isrigid':True,'cutpoints':[],'uncrossed':[],'virtuallygeometric':False},
'Baumslag':{'freegroup':F2,'whiteheadgraph':Baumslag,'wordlist':[baumslag],'splitsfreely':False,'iscircle':False,'isrigid':False,'cutpoints':[],'uncrossed':[a], 'virtuallygeometric':True},
    }
