"""Independent check of the twelve C5 x C5 covers and the degree-171 companion (written from the text only).

For each cover (three ten-tuples in representatives.json):
  * the 35 developed tiles partition the 1050 edges of B = complement of Hoffman-Singleton, and every point
    lies in 7 tiles;
  * nauty on the point/tile/star incidence graph: the twelve covers are pairwise non-isomorphic, with groups
    of order 50 or 100;
  * there is exactly one set P of five tiles whose supports partition the 50 points;
  * a random decoration gives srg(400,168,68,72); keeping X on P, replacing X by -X on the other 30 tiles and
    adding the two K4s split by the P-character in each fiber gives srg(400,171,74,72);
  * the zero-sum 4-sets (rows of A for degree 168, of A+I for degree 171) are exactly 700 affine planes.
Optionally compares with graph6 files of srg(400,171,74,72) given on the command line (zero-sum 4-set count).
Usage: independent_covers.py [other_400_171.g6 ...]
"""
import json
import random
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
import pynauty

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / 'srg496' / 'family'))
K5 = list(combinations(range(5), 2))


def zero_sum_quads(A):
    """4-sets of vertices whose rows of A+I sum to 0 or 1 mod 2 (A+I is formed here)."""
    n = len(A)
    M = (A + np.eye(n, dtype=np.int64)) % 2
    rows = [int(''.join(map(str, r)), 2) for r in M]
    full = (1 << n) - 1
    buckets = defaultdict(list)
    for a, b in combinations(range(n), 2):
        x = rows[a] ^ rows[b]
        buckets[min(x, x ^ full)].append((a, b))
    quads = set()
    for prs in buckets.values():
        for (a, b), (c, d) in combinations(prs, 2):
            s = {a, b, c, d}
            if len(s) == 4:
                quads.add(tuple(sorted(s)))
    return quads


def pt(s, x, y):
    return 25 * s + 5 * (x % 5) + (y % 5)


def hoffman_singleton():
    H = np.zeros((50, 50), dtype=np.int64)
    for x in range(5):
        for y in range(5):
            for dy, s in ((1, 0), (2, 1)):
                a, b = pt(s, x, y), pt(s, x, y + dy)
                H[a, b] = H[b, a] = 1
            for u in range(5):
                v = (y - 3 * (u - x) ** 2) % 5
                a, b = pt(0, x, y), pt(1, u, v)
                H[a, b] = H[b, a] = 1
    assert np.all(H.sum(1) == 7)
    H2 = H @ H
    assert np.all(H2[H == 1] == 0) and np.all(H2[(H == 0) & ~np.eye(50, dtype=bool)] == 1)
    return H


def translate(p, a, b):
    s, r = divmod(p, 25)
    return pt(s, r // 5 + a, r % 5 + b)


def signing():
    D = np.zeros((5, 10), dtype=np.int64)
    for j, (u, v) in enumerate(K5):
        i, k = (u, v) if (v - u) % 5 in (1, 2) else (v, u)
        D[i, j], D[k, j] = 1, -1
    return 2 * np.eye(10, dtype=np.int64) - D.T @ D


def tiles_of(starters, B):
    X = signing()
    tiles = {}
    for st in starters:
        for a in range(5):
            for b in range(5):
                pts = [translate(p, a, b) for p in st]
                key = frozenset(pts)
                if key not in tiles:
                    tiles[key] = pts                    # point list in K5-edge order; X is the signing
    tiles = list(tiles.values())
    assert len(tiles) == 35
    cover = np.zeros((50, 50), dtype=np.int64)
    for pts in tiles:
        for i, j in combinations(range(10), 2):
            if X[i, j]:
                cover[pts[i], pts[j]] += 1
                cover[pts[j], pts[i]] += 1
    assert np.array_equal(cover, B), 'tiles do not partition B'
    assert all(sum(p in t for t in tiles) == 7 for p in range(50))
    return tiles, X


def incidence_cert(tiles):
    """Point / tile / star incidence graph (stars = the five symbol stars of each tile)."""
    adj = defaultdict(set)
    n = 50 + 35
    stars = []
    for t, pts in enumerate(tiles):
        for u in range(5):
            s = n + len(stars)
            stars.append(s)
            adj[s] |= {50 + t} | {pts[j] for j, e in enumerate(K5) if u in e}
    N = n + len(stars)
    for s in list(adj):
        for v in adj[s]:
            adj[v].add(s)
    g = pynauty.Graph(N, adjacency_dict={v: sorted(w) for v, w in adj.items()},
                      vertex_coloring=[set(range(50)), set(range(50, 85)), set(range(85, N))])
    _, size, exp, _, _ = pynauty.autgrp(g)
    return pynauty.certificate(g), round(size * 10 ** exp)


def partitions_P(tiles):
    sup = [frozenset(t) for t in tiles]
    out = []
    for c in combinations(range(35), 5):
        if len(frozenset().union(*(sup[i] for i in c))) == 50:
            out.append(c)
    return out


def chi(c, x):
    return -1 if bin(c & x).count('1') % 2 else 1


def lift(tiles, X, rng, P=None):
    """Random decoration; with P given, the degree-171 companion."""
    lab, sgn, mine = {}, {}, defaultdict(list)
    for t, pts in enumerate(tiles):
        for p in pts:
            mine[p].append(t)
    for p in range(50):
        cs = list(range(1, 8))
        rng.shuffle(cs)
        for t, c in zip(mine[p], cs):
            lab[p, t], sgn[p, t] = c, rng.choice((1, -1))
    A = np.zeros((400, 400), dtype=np.int64)
    for t, pts in enumerate(tiles):
        f = 1 if (P is None or t in P) else -1
        for i, j in combinations(range(10), 2):
            if not X[i, j]:
                continue
            p, q = pts[i], pts[j]
            for x in range(8):
                for y in range(8):
                    if sgn[p, t] * chi(lab[p, t], x) * sgn[q, t] * chi(lab[q, t], y) == f * X[i, j]:
                        A[8 * p + x, 8 * q + y] = A[8 * q + y, 8 * p + x] = 1
    if P is not None:
        for p in range(50):
            tp = next(t for t in mine[p] if t in P)
            for x, y in combinations(range(8), 2):
                if chi(lab[p, tp], x) == chi(lab[p, tp], y):
                    A[8 * p + x, 8 * p + y] = A[8 * p + y, 8 * p + x] = 1
    return A


def srg_ok(A, k, lam, mu):
    assert np.array_equal(A, A.T) and not A.diagonal().any() and np.all(A.sum(1) == k)
    assert np.array_equal(A @ A, (k - mu) * np.eye(400, dtype=np.int64) + (lam - mu) * A + mu)


def read_g6(path):
    data = Path(path).read_bytes().split(b'\n')[0].strip()
    if data.startswith(b'>>graph6<<'):
        data = data[10:]
    if data[0] == 126:
        n, body = ((data[1] - 63) << 12) | ((data[2] - 63) << 6) | (data[3] - 63), data[4:]
    else:
        n, body = data[0] - 63, data[1:]
    bits = [((c - 63) >> s) & 1 for c in body for s in range(5, -1, -1)]
    A = np.zeros((n, n), dtype=np.int64)
    k = 0
    for j in range(n):
        for i in range(j):
            A[i, j] = A[j, i] = bits[k]
            k += 1
    return A


def main():
    rng = random.Random(400171)
    B = 1 - np.eye(50, dtype=np.int64) - hoffman_singleton()
    reps = json.loads((HERE / 'representatives.json').read_text())
    assert len(reps) == 12
    certs = []
    for c, starters in enumerate(reps):
        tiles, X = tiles_of(starters, B)
        cert, order = incidence_cert(tiles)
        certs.append(cert)
        Ps = partitions_P(tiles)
        assert len(Ps) == 1
        A168 = lift(tiles, X, rng)
        srg_ok(A168, 168, 68, 72)
        A171 = lift(tiles, X, rng, P=set(Ps[0]))
        srg_ok(A171, 171, 74, 72)
        q168 = len(zero_sum_quads(A168 + np.eye(400, dtype=np.int64)))   # rows of A
        q171 = len(zero_sum_quads(A171))                                   # rows of A + I
        assert q168 == q171 == 700
        print(f'cover {c:2d}: tiles partition B, |group| = {order}, unique P, srg 168 and 171 ok, '
              f'zero-sum 4-sets {q168}/{q171}', flush=True)
    assert len(set(certs)) == 12
    print('the twelve covers are pairwise non-isomorphic')
    for path in sys.argv[1:]:
        A = read_g6(path)
        srg_ok(A, 171, 74, 72)
        print(f'{Path(path).name}: srg(400,171,74,72), zero-sum 4-sets of A+I: {len(zero_sum_quads(A))}')


if __name__ == '__main__':
    main()
