"""Shared helpers for the 400 audit (independent of the original constructor)."""
import itertools, random
import numpy as np
import os
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/'
STARTERS = [list(map(int, l.split())) for l in open(SRC + 'starters.txt') if l.strip()]
P = lambda i, j: 5 * (i % 5) + (j % 5)
Q = lambda i, j: 25 + 5 * (i % 5) + (j % 5)
def hs():
    H = np.zeros((50, 50), int)
    for i in range(5):
        for j in range(5):
            for d, f in ((1, P), (2, Q)):
                H[f(i, j), f(i, j + d)] = H[f(i, j + d), f(i, j)] = 1
            for u in range(5):
                v = (j - 3 * (u - i) ** 2) % 5
                H[P(i, j), Q(u, v)] = H[Q(u, v), P(i, j)] = 1
    return H
def tr(lab, a, b):
    s, x = divmod(lab, 25); i, j = divmod(x, 5); return 25 * s + 5 * ((i + a) % 5) + (j + b) % 5
K5 = list(itertools.combinations(range(5), 2))
D = np.zeros((5, 10), int)
for e, (i, j) in enumerate(K5):
    t, h = (i, j) if (j - i) % 5 in (1, 2) else (j, i)
    D[t, e] = -1; D[h, e] = 1
X = 2 * np.eye(10, dtype=int) - D.T @ D
def tiles(starters=STARTERS):
    out = {}
    for row in starters:
        for a in range(5):
            for b in range(5):
                lab = [tr(p, a, b) for p in row]
                M = np.zeros((50, 50), int)
                for x in range(10):
                    for y in range(10): M[lab[x], lab[y]] = X[x, y]
                key = frozenset(lab)
                if key in out: assert (out[key] == M).all()
                else: out[key] = M
    return out
T = tiles(); TL = sorted(T, key=lambda k: sorted(k)); TI = {k: i for i, k in enumerate(TL)}
INC = {p: [k for k in TL if p in k] for p in range(50)}
TILE_OF = {}
for k in TL:
    for p in k:
        for q in k:
            if T[k][p, q]: TILE_OF[p, q] = k
CHI = np.array([[(-1) ** bin(c & x).count('1') for x in range(8)] for c in range(8)])
def graph(assign, sign=None):
    A = np.zeros((400, 400), np.int8)
    for (p, q), k in TILE_OF.items():
        fp = CHI[assign[p, k]] * (sign[p, k] if sign else 1)
        fq = CHI[assign[q, k]] * (sign[q, k] if sign else 1)
        A[8*p:8*p+8, 8*q:8*q+8] = (np.outer(fp, fq) == T[k][p, q])
    return A
def is_srg(A, k=168, l=68, m=72):
    A = A.astype(np.int64); n = len(A); A2 = A @ A; I = np.eye(n, dtype=np.int64)
    return (A == A.T).all() and not A.diagonal().any() and (A2 == k * I + l * A + m * (1 - I - A)).all()
def origin_assign(perm0, perm25, sign0=None, sign25=None):
    """Translation-invariant data: perm_o[c-1] = character of the c-th tile at origin o (tiles sorted as TL order)."""
    asg, sg = {}, {}
    for o, perm, sgn in ((0, perm0, sign0), (25, perm25, sign25)):
        inc0 = INC[o]
        for a in range(5):
            for b in range(5):
                for c, k in enumerate(inc0):
                    kk = frozenset(tr(p, a, b) for p in k)
                    asg[tr(o, a, b), kk] = perm[c]; sg[tr(o, a, b), kk] = sgn[c] if sgn else 1
    return asg, sg
def random_assign(rng, signs=True):
    asg, sg = {}, {}
    for p in range(50):
        for c, k in zip(rng.sample(range(1, 8), 7), INC[p]):
            asg[p, k] = c; sg[p, k] = rng.choice((1, -1)) if signs else 1
    return asg, sg
def certificate(A):
    import pynauty
    g = pynauty.Graph(len(A), adjacency_dict={i: list(np.flatnonzero(A[i])) for i in range(len(A))})
    return pynauty.certificate(g)
def autorder(A):
    import pynauty
    g = pynauty.Graph(len(A), adjacency_dict={i: list(np.flatnonzero(A[i])) for i in range(len(A))})
    _, m, e, orbs, no = pynauty.autgrp(g); return round(m * 10 ** e), no
# ---- 15-cliques via F2 systems over the 100 HS 15-cocliques ----
def hs_cocliques():
    import igraph as ig
    H = hs(); g = ig.Graph.Adjacency(H.tolist(), mode='undirected')
    return [tuple(sorted(s)) for s in g.largest_independent_vertex_sets()]
def bits(c): return [(c >> i) & 1 for i in range(3)]
def f2_rank_solve(rows, rhs):
    """Gaussian elimination over F2. rows: list of int bitmasks; return (consistent, rank)."""
    piv = {}
    for r, b in zip(rows, rhs):
        while r:
            h = r.bit_length() - 1
            if h in piv: pr, pb = piv[h]; r ^= pr; b ^= pb
            else: piv[h] = (r, b); break
        else:
            if b: return False, None
    return True, len(piv)
def fifteen_cliques(S, asg, sg):
    """Number of 15-cliques of the lift lying over the B-clique S (15 points)."""
    idx = {p: i for i, p in enumerate(S)}
    rows, rhs = [], []
    for p, q in itertools.combinations(S, 2):
        k = TILE_OF[p, q]
        r = 0
        for t, bt in enumerate(bits(asg[p, k])): r |= bt << (3 * idx[p] + t)
        for t, bt in enumerate(bits(asg[q, k])): r |= bt << (3 * idx[q] + t)
        # adjacency: chi_p(x) chi_q(y) s_p s_q = X  <=>  c_p.x + c_q.y = [s_p s_q X == -1]
        rows.append(r); rhs.append(int(sg[p, k] * sg[q, k] * T[k][p, q] == -1))
    ok, rk = f2_rank_solve(rows, rhs)
    return 2 ** (45 - rk) if ok else 0
