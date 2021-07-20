# imports Whitehead graph stuff and sets up some basic examples.
import grouptheory.group
import freegroup 
import networkx as nx
import AutF as aut
import graphofgroups as gog
import whiteheadgraph as wg
from whiteheadgraph.test.knownexamples import *
import whiteheadgraph.test.checkcutpair
import whiteheadgraph.test.rjsj

whiteheadgraph.test.checkcutpair.testall()
whiteheadgraph.test.rjsj.testall()
