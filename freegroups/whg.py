# imports Whitehead graph stuff, sets up some basic examples, and runs some tests.
# in ipython do:
# %run grouptheory/freegroups/whg.py
import grouptheory.group
import grouptheory.freegroups.freegroup 
import networkx as nx
import grouptheory.freegroups.AutF as aut
import grouptheory.freegroups.graphofgroups as gog
import grouptheory.freegroups.whiteheadgraph as wg
from grouptheory.freegroups.whiteheadgraph.test.knownexamples import *
import grouptheory.freegroups.whiteheadgraph.test.checkcutpair as checkcutpair
import grouptheory.freegroups.whiteheadgraph.test.rjsj as rjsj

checkcutpair.testall()
rjsj.testall()
