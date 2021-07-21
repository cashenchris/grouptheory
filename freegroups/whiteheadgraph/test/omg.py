import networkx as nx
from networkx.classes.graph import Graph
from networkx import NetworkXException, NetworkXError
import networkx.convert as convert
from copy import deepcopy
import grouptheory.freegroups.whiteheadgraph.orderedmultigraph as omg
import copy

g=omg.OrderedMultiGraph()
g.add_edge(1,2,'e')
g.add_edge(1,2,'f1')
g.add_edge(1,2,'f2')
g.add_edge(2,3,'g')
g.add_vertex(4)
g.add_edge(2,3)
g.add_edge(2,3)
g.add_edge(2,3, 'h', 1, 1)
g.remove_edge('f1')

circle=omg.OrderedMultiGraph()
circle.add_edge(1,2,'e1')
circle.add_edge(2,3,'e2')
circle.add_edge(3,4,'e3')
circle.add_edge(4,1,'e4')

circle.is_connected()
circle.is_circle()

circle.find_cut_vertex()

h=omg.splice(g,circle,1,1,[0,1],('g',),('c',),lookforisolatedvertices=True)
#k,foundfree=omg.selfsplice(circle,1,2,[0,1],lookforisolatedvertices=True,reportfreeedges=True)
