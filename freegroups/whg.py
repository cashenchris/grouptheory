# imports Whitehead graph stuff and sets up some basic examples.
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
