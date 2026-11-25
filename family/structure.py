"""Certify fibre recognition, clique projection bounds, and the cover automorphism group."""
from collections import Counter
from itertools import combinations
from math import comb

import igraph as ig
import pynauty

from family import tile_data


def recover(ports):
    # Entries of column[q][p] identify the tile carrying the block (p,q).
    column = [[sum(1<<t for t,b in enumerate(neighbours) if b>>p&1)
               for p in range(50)] for neighbours in ports]
    for n in (2,4):
        tested = 0
        for odd in combinations(range(50),n):
            for col in column:
                used = 0
                for p in odd:
                    used ^= col[p]
                if used.bit_count() in (1,2):
                    break
            else:
                raise AssertionError(('unresolved odd fibres', odd))
            tested += 1
        assert tested == comb(50,n)
        print('odd-fibre sets with an obstruction',n,tested)
    supports = [sum(s) for row in ports for s in combinations(row,4)]
    assert len(supports) == len(set(supports)) == 1750
    print('different four-direction neighbourhood supports', len(set(supports)))


def group(tiles):
    # The three colours distinguish the 50 points, 35 tiles, and 1050 edges.
    graph = {i:[] for i in range(1135)}
    edges = []
    node = 85
    for t,tile in enumerate(tiles):
        for p,q,_ in tile:
            for v in (p,q,50+t):
                graph[node].append(v)
                graph[v].append(node)
                edges.append((v,node))
            node += 1
    assert node == 1135
    colors = [set(range(50)), set(range(50,85)), set(range(85,1135))]
    g = pynauty.Graph(node,adjacency_dict=graph,vertex_coloring=colors)
    generators,size,exp,orbits,_ = pynauty.autgrp(g)
    order = round(size*10**exp)
    assert order == ig.Graph(node,edges=edges).count_automorphisms(color=[0]*50+[1]*35+[2]*1050) == 50
    counts = Counter(orbits[:50])
    assert sorted(counts.values()) == [25,25]
    print('cover automorphism group order',order,'point orbits',sorted(counts.values()))
    from family import construct
    translations = {tuple(construct.shift(p,d) for p in range(50)) for d in construct.TRANSLATIONS}
    outer = next(g for g in generators if tuple(g[:50]) not in translations)
    r,c = divmod(outer[0],5)
    lin = [construct.shift(outer[i],(-r,-c)) for i in range(25)]
    assert all(lin[construct.shift(u,divmod(v,5))] == construct.shift(lin[u],divmod(lin[v],5))
               for u in range(25) for v in range(25))
    assert sum(lin[i] == i for i in range(25)) == 5
    assert sum(construct.shift(i,divmod(lin[i],5)) == 0 for i in range(25)) == 5
    print('the involution on the translation plane has one fixed and one negated dimension')
    return outer


def clique_bound(tiles, ports):
    adjacent = [sum(s) for s in ports]
    forbidden = [[0]*50 for _ in ports]
    stars = 0
    for tile in tiles:
        support = {(p,q) for p,q,_ in tile}
        nodes = sorted({k for pair in support for k in pair})
        for four in combinations(nodes,4):
            if not all(pair in support for pair in combinations(four,2)):
                continue
            stars += 1
            for i,j in combinations(four,2):
                forbidden[i][j] |= sum(1<<v for v in four if v not in (i,j))
                forbidden[j][i] = forbidden[i][j]
    assert stars == 175
    # Branch in increasing order, deleting nonneighbours and forbidden star triples.
    counts = Counter()
    def visit(chosen,candidates):
        counts[len(chosen)] += 1
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            p = bit.bit_length()-1
            rest = candidates & adjacent[p]
            for q in chosen:
                rest &= ~forbidden[q][p]
            visit(chosen+[p],rest)
    visit([], (1<<50)-1)
    print('admissible projected cliques by size',sorted(counts.items()))
    assert max(counts) == 10 and counts[10] == 95


def main():
    tiles, _ = tile_data()
    ports = [[] for _ in range(50)]
    for tile in tiles:
        neighbours = {}
        for p,q,_ in tile:
            neighbours[p] = neighbours.get(p,0) | (1<<q)
            neighbours[q] = neighbours.get(q,0) | (1<<p)
        for p,s in neighbours.items():
            ports[p].append(s)
    assert all(len(ns)==7 and sum(ns).bit_count()==42 for ns in ports)
    recover(ports)
    group(tiles)
    clique_bound(tiles,ports)
    bound = 30**50 * 2**165 // 50
    print('lower bound on unlabelled SRGs',bound)


if __name__ == '__main__':
    main()
