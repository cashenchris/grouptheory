#!/usr/bin/env python3

## usage:  ./irank.py word
## compute the imprimitivity rank of the given word in a free group
## word is a string where a-z represent generators of free group and A-Z their inverses

## examples:
## $ ./irank.py aaaa
## 1
## $ ./irank.py xyXY
## 2
## $ ./irank.py aabbcc
## 3

import grouptheory.freegroups.freegroup as freegroup
import grouptheory.freegroups.imprimitivity_rank as ir
import sys

####

F=freegroup.FGFreeGroup(numgens=26)
if len(sys.argv)>1:
    w=F.word(sys.argv[1])
    print(ir.imprimitivityrank(w))
else:
    print("Usage: ./irank.py word")
