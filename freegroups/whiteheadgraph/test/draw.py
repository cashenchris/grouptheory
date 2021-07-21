import networkx as nx
import matplotlib as plt
print("networkx drawtest")
g=nx.dodecahedral_graph()
try: # draw
        nx.draw(g)
except: # matplotlib not available
        pass
