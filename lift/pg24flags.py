"""Flag graph of PG(2,4): (p,L) ~ (p',L') iff exactly one of p' in L, p in L' holds (p != p', L != L')."""
import itertools, numpy as np, networkx as nx, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from gf import GF
F = GF(4)
vecs = [v for v in itertools.product(range(4), repeat=3) if any(v)]
def norm(v):
    i = next(i for i in range(3) if v[i]); inv = next(x for x in range(1, 4) if F.mul[v[i], x] == 1)
    return tuple(int(F.mul[inv, c]) for c in v)
pts = sorted({norm(v) for v in vecs}); assert len(pts) == 21
dot = lambda a, b: int(F.add[F.add[F.mul[a[0], b[0]], F.mul[a[1], b[1]]], F.mul[a[2], b[2]]])
flags = [(i, j) for i in range(21) for j in range(21) if dot(pts[i], pts[j]) == 0]
assert len(flags) == 105
inc = lambda i, j: dot(pts[i], pts[j]) == 0
A = np.zeros((105, 105), int)
for a, (p, L) in enumerate(flags):
    for b, (q, M) in enumerate(flags):
        if p != q and L != M and (inc(q, L) != inc(p, M)): A[a, b] = 1
A2 = A @ A; k = A[0].sum()
print('degree', k, 'lambda', set(A2[A == 1]), 'mu', set(A2[(A == 0) & ~np.eye(105, dtype=bool)]))
nx.write_graph6(nx.from_numpy_array(A), 'pg24flags.g6', header=False)
