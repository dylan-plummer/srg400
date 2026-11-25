"""Compute group orders with nauty and bliss, count four-cliques, and find clique numbers."""
from collections import Counter
from pathlib import Path

import igraph as ig
import networkx as nx
import pynauty

HERE = Path(__file__).resolve().parent
GRAPHS = HERE.parent / 'graphs'


def four_cliques(g):
    neighbours = [sum(1 << j for j in g[i]) for i in g]
    edge_counts = Counter()
    for i,j in g.edges():
        common = neighbours[i] & neighbours[j]
        count = 0
        while common:
            bit = common & -common
            common ^= bit
            count += (neighbours[bit.bit_length()-1] & common).bit_count()
        edge_counts[count] += 1
    total = sum(count * freq for count,freq in edge_counts.items())
    assert total % 6 == 0
    return total//6, edge_counts


def main():
    for suffix in ('', '_swap'):
        name = f'srg_400_168_68_72{suffix}.g6'
        g = nx.from_graph6_bytes((GRAPHS/name).read_bytes().strip())
        graph = ig.Graph(n=400, edges=list(g.edges()), directed=False)
        h = pynauty.Graph(400, adjacency_dict={i:list(g[i]) for i in g})
        _, mantissa, exponent, orbits, count = pynauty.autgrp(h)
        order = round(mantissa * 10**exponent)
        assert order == graph.count_automorphisms()
        roots = sorted(set(orbits))
        print(name, 'group order from nauty and bliss', order, 'vertex orbits', count, flush=True)
        print('orbit size distribution', dict(Counter(orbits.count(r) for r in roots)), flush=True)
        total, edge_counts = four_cliques(g)
        print('four-cliques', total, 'edge count range', min(edge_counts), max(edge_counts), flush=True)
        # Every clique containing a vertex consists of that vertex plus a clique in its neighbourhood.
        numbers = Counter(1 + graph.induced_subgraph(list(g[r])).clique_number() for r in roots)
        print('largest clique through a representative, counted by orbits', dict(sorted(numbers.items())),
              'clique number', max(numbers), flush=True)


if __name__ == '__main__':
    main()
