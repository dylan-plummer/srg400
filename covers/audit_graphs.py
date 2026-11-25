"""Recover the cover and the companion using only ordinary graph6 input."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import random

import igraph as ig
import networkx as nx

HERE = Path(__file__).resolve().parent


class Components:
    def __init__(self, vertices):
        self.parent = {v: v for v in vertices}

    def root(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def join(self, x, y):
        self.parent[self.root(x)] = self.root(y)

    def parts(self):
        groups = defaultdict(list)
        for x in self.parent:
            groups[self.root(x)].append(x)
        return sorted(sorted(group) for group in groups.values())


def read_graph(path):
    g = nx.from_graph6_bytes(path.read_bytes().strip())
    assert len(g) == 400
    return [sum(1 << j for j in g[i]) for i in range(400)]


def check_parameters(rows, degree):
    edges = [(i, j) for i, j in combinations(range(400), 2) if rows[i] >> j & 1]
    assert all(r.bit_count() == degree and not r >> i & 1 for i, r in enumerate(rows))
    counts = Counter(((rows[i] >> j) & 1, (rows[i] & rows[j]).bit_count())
                     for i, j in combinations(range(400), 2))
    assert len(edges) == 200 * degree
    assert counts == {(1, 68 if degree == 168 else 74): len(edges),
                      (0, 72): 79800 - len(edges)}
    plain = ig.Graph(400, edges=edges)
    return plain.count_automorphisms()


def planes(rows, degree):
    words = [r | ((1 << i) if degree == 171 else 0) for i, r in enumerate(rows)]
    assert len(set(words)) == 400
    first = {}
    collisions = defaultdict(list)
    for a, b in combinations(range(400), 2):
        z = words[a] ^ words[b]
        if z in first:
            collisions[z].append((a, b))
        else:
            first[z] = (a, b)
    fibres = Components(range(400))
    count = 0
    for z, other in collisions.items():
        pairs = [first[z], *other]
        for (a, b), (c, d) in combinations(pairs, 2):
            assert len({a, b, c, d}) == 4
            fibres.join(a, b)
            fibres.join(a, c)
            fibres.join(a, d)
            count += 1
    assert count % 3 == 0
    return count // 3, fibres.parts()


def cover_from_blocks(rows, fibres):
    masks = [sum(1 << x for x in f) for f in fibres]
    edge_ports = defaultdict(list)
    edges = []
    for p, q in combinations(range(50), 2):
        parts = {rows[x] & masks[q] for x in fibres[p]}
        if parts == {0}:
            continue
        edge = (p, q)
        edges.append(edge)
        for a, b in (edge, edge[::-1]):
            split = {rows[x] & masks[b] for x in fibres[a]}
            assert len(split) == 2 and all(x.bit_count() == 4 for x in split)
            assert sum(split) == masks[b]
            edge_ports[b, min(split)].append(edge)
    assert len(edges) == 1050 and len(edge_ports) == 350
    tiles = Components(edges)
    for group in edge_ports.values():
        assert len(group) == 6
        for edge in group[1:]:
            tiles.join(group[0], edge)
    tiles = tiles.parts()
    assert len(tiles) == 35 and all(len(t) == 30 for t in tiles)
    return tiles, masks, edge_ports


def cover_signature(tiles):
    edges = []
    for t, tile in enumerate(tiles):
        for e, (p, q) in enumerate(tile, 85 + 30*t):
            edges.extend([(p, e), (q, e), (50+t, e)])
    g = ig.Graph(1135, edges=edges)
    colours = [0]*50 + [1]*35 + [2]*1050
    canon = g.permute_vertices(g.canonical_permutation(color=colours))
    signature = sha256(repr(sorted(canon.get_edgelist())).encode()).hexdigest()
    return signature, g.count_automorphisms(color=colours)


def mate(rows, degree, tiles, masks, edge_ports, fibres):
    supports = [{p for edge in tile for p in edge} for tile in tiles]
    assert all(len(s) == 10 for s in supports)
    def partitions(remaining, chosen):
        if not remaining:
            yield chosen
            return
        p = min(remaining)
        for i, s in enumerate(supports):
            if p in s and s <= remaining:
                yield from partitions(remaining - s, chosen | {i})
    parallel, = partitions(set(range(50)), set())
    assert len(parallel) == 5
    result = list(rows)
    for t, tile in enumerate(tiles):
        if t in parallel:
            for p in supports[t]:
                halves = [m for (q, m), edges in edge_ports.items()
                          if p == q and edges[0] in tile]
                a, = halves
                for x in fibres[p]:
                    expected = (a if a >> x & 1 else masks[p] ^ a) ^ (1 << x)
                    assert rows[x] & masks[p] == (expected if degree == 171 else 0)
                    result[x] ^= expected
        else:
            for p, q in tile:
                for x in fibres[p]:
                    result[x] ^= masks[q]
                for y in fibres[q]:
                    result[y] ^= masks[p]
    return result


def encode(rows):
    graph = nx.Graph()
    graph.add_nodes_from(range(400))
    graph.add_edges_from((i,j) for i,j in combinations(range(400),2) if rows[i] >> j & 1)
    return nx.to_graph6_bytes(graph, header=False)


def main():
    census = json.loads((HERE/'equivariant.json').read_text())
    signatures = set()
    print('name\tdegree\tcommon_adjacent\tcommon_nonadjacent\tplanes\tgraph_aut\tcover_aut\tcover_bliss_sha256\tgraph6_sha256', flush=True)
    for rec in census:
        i = rec['cover']
        cert = None
        for kind in ('cover', 'sym'):
            stem = HERE / 'graphs' / f'{kind}_{i:02}_400'
            lower = read_graph(Path(str(stem) + '_168.g6'))
            upper_path = Path(str(stem) + '_171.g6')
            for degree in (168, 171):
                path = Path(str(stem) + f'_{degree}.g6')
                rows = read_graph(path)
                group = check_parameters(rows, degree)
                if kind == 'sym':
                    assert group == rec['max_group']
                count, fibres = planes(rows, degree)
                assert count == 700 and sorted(map(len, fibres)) == [8]*50
                tiles, masks, ports = cover_from_blocks(rows, fibres)
                sig, order = cover_signature(tiles)
                if cert is None:
                    cert = sig
                    assert sig not in signatures
                    signatures.add(sig)
                assert cert == sig
                if i == 10 and kind == 'sym' and degree == 171:
                    shuffled = list(range(400))
                    random.Random(20261125).shuffle(shuffled)
                    permuted = [sum(1 << j for j, y in enumerate(shuffled) if rows[x] >> y & 1)
                                for x in shuffled]
                    n, randomized_fibres = planes(permuted, degree)
                    assert n == 700
                    t, _, _ = cover_from_blocks(permuted, randomized_fibres)
                    assert cover_signature(t) == (sig, order)
                companion = mate(rows, degree, tiles, masks, ports, fibres)
                if degree == 168:
                    if kind == 'sym' and not upper_path.exists():
                        upper_path.write_bytes(encode(companion))
                    assert companion == read_graph(upper_path)
                    lower_order = group
                else:
                    assert companion == lower and group == lower_order
                checksum = sha256(path.read_bytes()).hexdigest()
                print(path.stem, degree, 68 if degree == 168 else 74, 72, count, group, order, sig, checksum, sep='\t', flush=True)
    root = HERE.parents[1]
    public = ['koval/new_srgs_multiorbit_data/out/new_srgs_multiorbit/srg400_171_74_72.g6',   # Kovaľ preprint (Zenodo 23296158)
              'koval/srg-open-cases/data/srg400_171_74_72_2f5e8c0d6d124206.g6']                # github.com/tonykoval/srg-open-cases
    for path in [root / q for q in public if (root / q).exists()]:
        rows = read_graph(path)
        group = check_parameters(rows, 171)
        count, _ = planes(rows, 171)
        print(path.relative_to(root), 171, 74, 72, count, group, 'external','external', sha256(path.read_bytes()).hexdigest(), sep='\t', flush=True)
    assert len(signatures) == 12


if __name__ == '__main__':
    main()
