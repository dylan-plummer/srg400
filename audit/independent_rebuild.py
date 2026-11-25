"""Independent rebuild of the srg(400,168,68,72) graphs from the written construction only (NumPy).
Rebuilds Hoffman-Singleton, the 35 tiles, both exported graphs (exact comparison), and 20 random members."""
import sys, itertools, random
import numpy as np
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
def tiles(starters=None):
    starters = starters or globals()['starters']
    """List of (labels tuple, X) for distinct developed tiles; also check sign consistency."""
    out = {}
    for row in starters:
        for a in range(5):
            for b in range(5):
                lab = [tr(p, a, b) for p in row]
                M = np.zeros((50, 50), int)
                for x in range(10):
                    for y in range(10): M[lab[x], lab[y]] = X[x, y]
                key = frozenset(lab)
                if key in out: assert (out[key] == M).all(), 'sign mismatch'
                else: out[key] = M
    return out
def graph(T, assign, sign=None):
    """assign[(p, key)] = character c in 1..7; sign[(p,key)] = +-1."""
    xs = np.arange(8)
    chi = np.array([[(-1) ** bin(c & x).count('1') for x in xs] for c in range(8)])
    A = np.zeros((400, 400), int)
    for key, M in T.items():
        for p in key:
            for q in key:
                if M[p, q] == 0: continue
                fp = chi[assign[p, key]] * (sign[p, key] if sign else 1)
                fq = chi[assign[q, key]] * (sign[q, key] if sign else 1)
                A[8*p:8*p+8, 8*q:8*q+8] = (np.outer(fp, fq) == M[p, q])
    return A
def is_srg(A, k=168, l=68, m=72):
    n = len(A); A2 = A @ A; I = np.eye(n, dtype=int); J = np.ones((n, n), int)
    return (A == A.T).all() and not A.diagonal().any() and (A2 == k * I + l * A + m * (J - I - A)).all()
import os
HERE = os.path.dirname(os.path.abspath(__file__))
starters = [list(map(int, l.split())) for l in open(os.path.join(HERE, '..', 'data', 'starters.txt')) if l.strip()]
if __name__ == '__main__':
    H = hs(); I = np.eye(50, dtype=int); J = np.ones((50, 50), int)
    print('HS srg(50,7,0,1):', (H.sum(1) == 7).all() and (H @ H == 7 * I + 0 * H + 1 * (J - I - H)).all())
    B = J - I - H
    T = tiles(starters)
    S = sum(abs(M) for M in T.values())
    print('tiles', len(T), 'partition of B:', (S == B).all(), 'X spectrum', np.round(np.linalg.eigvalsh(X)).astype(int).tolist())
    inc = {p: [k for k in T if p in k] for p in range(50)}
    print('tiles per point', set(len(v) for v in inc.values()))
    # README labelling
    sup = lambda k: sorted(tuple(sorted((p, q))) for p in k for q in k if p < q and T[k][p, q] != 0)
    assign = {}
    for o in (0, 25):
        inc0 = sorted(inc[o], key=sup)
        for a in range(5):
            for b in range(5):
                for c, k in enumerate(inc0, 1):
                    assign[tr(o, a, b), frozenset(tr(p, a, b) for p in k)] = c
    import networkx as nx
    base = os.path.join(HERE, '..', 'graphs')
    for suf in ('', '_swap'):
        asg = dict(assign)
        if suf:
            for k in inc[0]:
                if asg[0, k] in (1, 2): asg[0, k] = 3 - asg[0, k]
        A = graph(T, asg)
        F = nx.to_numpy_array(nx.from_graph6_bytes(open(f'{base}/srg_400_168_68_72{suf}.g6', 'rb').read().strip()), dtype=int)
        print(suf or 'base', 'srg:', is_srg(A), 'equals export:', (A == F).all())
    # random assignments + signs
    rng = random.Random(1)
    ok = 0
    for trial in range(20):
        asg = {}; sg = {}
        for p in range(50):
            perm = rng.sample(range(1, 8), 7)
            for c, k in zip(perm, inc[p]): asg[p, k] = c; sg[p, k] = rng.choice((1, -1))
        ok += is_srg(graph(T, asg, sg))
    print('random bijections + random signs: srg in', ok, 'of 20')
