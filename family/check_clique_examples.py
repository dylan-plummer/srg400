"""Read three named group-50 examples and compute maximum cliques in ordinary graphs."""
from pathlib import Path
import subprocess
import sys

import igraph as ig
import networkx as nx
import pynauty

HERE=Path(__file__).resolve().parent


def main():
    for number,want,name in ((384,10,'srg400_order50'),(386,8,'srg400_clique8'),(387,9,'srg400_clique9')):
        data=subprocess.check_output([sys.executable,str(HERE/'export.py'),str(number)])
        (HERE/(name+'.g6')).write_bytes(data)
        g=nx.from_graph6_bytes(data.strip())
        h=ig.Graph(400,edges=list(g.edges()))
        aut=pynauty.autgrp(pynauty.Graph(400,adjacency_dict={v:list(g[v]) for v in g}))
        assert aut[1]*10**aut[2]==50==h.count_automorphisms()
        omega=max(1+h.induced_subgraph(list(g[v])).clique_number() for v in set(aut[3]))
        assert omega==want
        print('class',number,'order',50,'clique number',omega)


if __name__=='__main__':
    main()
