import grouptheory.freegroups.freegroup as fg


                   
    
def IvanovSchupp(relator,reportreason=False):
    """
    Check if one-relator group defined by input relator is hyperbolic according to Ivanov and Schupp criteria.
    Returns True if group is known hyperbolic, False if known not hyperbolic, None if inconclusive

    >>> IvanovSchupp('',reportreason=True)
    (True, 'free')
    >>> IvanovSchupp('a',reportreason=True)
    (True, 'free')
    >>> IvanovSchupp('aa',reportreason=True)
    (True, 'torsion')
    >>> IvanovSchupp('abaCBC',reportreason=True)
    (False, 'Thm3(1)')
    >>> IvanovSchupp('ababcbc',reportreason=True)
    (True, 'Thm3(1)')
    >>> IvanovSchupp('abbAbcBC',reportreason=True)
    (True, 'Thm3(2)')
    >>> IvanovSchupp('bcBCCBAcbCBcbCBabc',reportreason=True)
    (False, 'Thm3(2a)')
    >>> IvanovSchupp('acaacacaacacaacbAAB',reportreason=True)
    (False, 'Thm3(2b)')
    >>> IvanovSchupp('ababbacb',reportreason=True)
    (True, 'Thm3(3)')
    >>> IvanovSchupp('ababaccBCbccBCbb',reportreason=True)
    (False, 'Thm3(3a)')
    >>> IvanovSchupp('abaBcbCCBcbCCbaccBCbccBCbb',reportreason=True)
    (False, 'Thm3(3b)')
    >>> IvanovSchupp('abaccBCbbaBcbCCb',reportreason=True)
    (False, 'Thm3(3c)')
    >>> IvanovSchupp('abaBcbCCbaBcbCCBcbCCb',reportreason=True)
    (False, 'Thm3(3d)')
    >>> IvanovSchupp('ababbcBCbcBCbcBCBAcbCB',reportreason=True)
    (True, 'Thm3(4)')
    >>> IvanovSchupp('ababbcBCBAcbCB',reportreason=True)
    (False, 'Thm3(4a)')
    >>> IvanovSchupp('ababbcBCbcBCBAcbCB',reportreason=True)
    (False, 'Thm3(4b)')
    >>> IvanovSchupp('ababcabccabcbcbcbcabcbcbcbcbc',reportreason=True)
    (True, 'Thm4')
    >>> IvanovSchupp('abacabcabCbc',reportreason=True)
    (False, 'Thm4(3)')
    >>> IvanovSchupp('ababcBCababccc',reportreason=True)
    (None, None)
    >>> IvanovSchupp('aaaaabbbbbccccc',reportreason=True)
    (None, None)
    """
    def format_return(result,reason):
        if reportreason:
            return result,reason
        else:
            return result
    if not relator:
        return format_return(True,'free')
    F,r1=fg.parseinputword(relator)
    r2=F.cyclic_reduce(r1)
    if r1!=r2:
        raise ValueError("Input relator is not cyclically reduced.")
    if F.degree(r2)>1:# one-relator groups with torsion are hyperbolic
        return format_return(True,'torsion')
    # Try Ivanov-Schupp criteria
    for a in range(1,1+F.rank):
        acount=[abs(x) for x in r2.letters].count(a)
        if acount==0: # neither theorem applies, check next a
            continue
        elif acount==1: # relator is primitive
            return format_return(True, 'free')
        # Ivanov-Schupp Theorem 3 applies when acount=2 or 3
        elif acount==2: 
            if r2.letters.count(a)==0:
                therel=(r2**(-1)).letters
            else:
                therel=r2.letters
            firsta=therel.index(a)
            therel=therel[firsta:]+therel[:firsta]
            if therel[1:].count(a)>0: # Ivanov-Schupp Thm 3 case (1)
                nexta=1+therel[1:].index(a)
                B=F.word(therel[1:nexta])
                if nexta<len(therel)-1:
                    C=F.word(therel[nexta+1:])
                else:
                    C=F.word([])
                if F.degree(B*C**(-1))>1:
                    return format_return(False,'Thm3(1)')
                else:
                    return format_return(True,'Thm3(1)')
            else: # Ivanov-Schupp Thm 3 case (2)
                nexta=1+therel[1:].index(-a)
                B=F.word(therel[1:nexta])
                C=F.word(therel[nexta+1:]) # no index error here since word is cyclically reduced
                if F.degree(B)>1 and F.degree(C)>1:
                    return format_return(False,'Thm3(2b)')
                elif F.is_conjugate_into(B,C) or F.is_conjugate_into(C,B):# group is not hyperbolic by case (2a), simplified because by previous case we know either B or C is not a proper power
                    return format_return(False,'Thm3(2a)')
                else:
                    return format_return(True,'Thm3(2)')
        elif acount==3:
            if r2.letters.count(-a)>r2.letters.count(a):
                therel=(r2**(-1)).letters
            else:
                therel=r2.letters
            if therel.count(-a)>0: # Theorem 3 case (4)
                firsta=therel.index(a)
                seconda=1+firsta+therel[1+firsta:].index(a)
                nega=therel.index(-a)
                if nega>seconda or nega<firsta:
                    therel=therel[firsta:]+therel[:firsta]
                else:
                    therel=therel[seconda:]+therel[:seconda]
                firsta=therel.index(a)
                seconda=1+firsta+therel[1+firsta:].index(a)
                nega=therel.index(-a)
                B=F.word(therel[1+firsta:seconda])
                C=F.word(therel[1+seconda:nega])
                D=F.word(therel[1+nega:])
                Z1,n1=F.max_root(B**(-1)*C*B)
                Z2,n2=F.max_root(D)
                if Z2==Z1**(-1):
                    Z2=Z1
                    n2=-n2
                if Z1==Z2 or Z1==F.word([]) or Z2==F.word([]):
                    if abs(n1)==abs(n2):
                        return format_return(False,'Thm3(4a)')
                    elif n1==-2*n2 or n2==-2*n1:
                        return format_return(False,'Thm3(4b)')
                    else:
                        return format_return(True,'Thm3(4)')
                else:
                    return format_return(True,'Thm3(4)')
            else: # Theorem 3 case (3)
                therel=therel[therel.index(a):]+therel[:therel.index(a)]
                firsta=0
                seconda=1+therel[1:].index(a)
                thirda=1+seconda+therel[1+seconda:].index(a)
                B=F.word(therel[1+firsta:seconda])
                C=F.word(therel[1+seconda:thirda])
                try:
                    D=F.word(therel[1+thirda:])
                except IndexError:
                    D=F.word([])
                Z1,n1=F.max_root(C*B**(-1))
                Z2,n2=F.max_root(D*B**(-1))
                if Z2==Z1**(-1):
                    Z2=Z1
                    n2=-n2
                if Z1==Z2 or Z1==F.word([]) or Z2==F.word([]):
                    if (n1==0 and abs(n2)>1) or (n2==0 and abs(n1)>1):
                        return format_return(False,'Thm3(3a)')
                    elif abs(n1)==abs(n2) and abs(n1)>1:
                        return format_return(False,'Thm3(3b)')
                    elif n1!=0 and n1==-n2:
                        return format_return(False,'Thm3(3c)')
                    elif n1!=0 and (n1==2*n2 or n2==2*n1):
                        return format_return(False,'Thm3(3d)')
                    else:
                        return format_return(True,'Thm3(3)')
                else:
                    return format_return(True,'Thm3(3)')
        else: # Theorem 3 does not apply. Theorem 4 applies to some words with higher acount. Need all a's same sign plus some other conditions.
            therel=r2.letters
            if therel.count(a)==0: # there are no a's, but we know there are some a's or a^-1's, so the inverse must have some a's
                therel=(r2**(-1)).letters
            if therel.count(a)>0 and therel.count(-a)==0: # only positive a's appear. Now check words between a's are all distinct.
                firstaindex=therel.index(a)
                R=F.word(therel[firstaindex:]+therel[:firstaindex]) #cylically permute r so that a is first letter
                T=[]
                currentword=[]
                for i in range(1,len(R)):
                    if R.letters[i]==a:
                        T.append(F.word(currentword))
                        currentword=[]
                    else:
                        currentword.append(R.letters[i])
                    if i==len(R)-1:
                        T.append(F.word(currentword))
                if len(T)==len(set(T)): #intermediate words are distinct, so theorem 4 applies
                    if R.letters.count(a)>4:
                        return format_return(True,'Thm4')
                    else: # check Ivanov-Schupp Thm 4 case (3)
                        assert(R.letters.count(a)==4) # if it were fewer we should have caught this in the theorem 3 part above
                        for i in range(4):
                            if T[i%4]*T[(i+1)%4]**(-1)*T[(i+2)%4]*T[(i+3)%4]**(-1)==F.word([]): 
                                return format_return(False,'Thm4(3)')
                        else: 
                            return format_return(True,'Thm4')
            else: # Theorem 4 does not apply with this a.
                pass
    # Ivanov-Schupp does not apply
    return format_return(None,None)


   








# in terminal, do
# python IvanovSchupp.py
# to run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()
