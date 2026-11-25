"""F=2, complete-fibre case of the signed-tile lift: a signing X of a base srg B with restricted eigenvalues
r_B, s_B such that X^2 = (rho+sigma) X - rho*sigma I with rho = r_B+2, sigma = s_B+2 (and B's degree = -rho*sigma).
Then  A = I (x) (J2-I2) + lift  is srg on 2v vertices. Encodes the condition as CP-SAT on edge sign bits.
usage: signing.py BASE {bh81 | g6file} time workers"""
import sys, itertools, time
import numpy as np
from ortools.sat.python import cp_model
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
def brouwer_haemers():
    from gf import GF
    F = GF(81)
    gen = next(x for x in range(1, 81) if F.order[x] == 80)
    D = set(); y = 1
    for i in range(80):
        if i % 4 == 0: D.add(y)
        y = int(F.mul[y, gen])
    A = np.zeros((81, 81), int)
    for a in range(81):
        for d in D: A[a, F.add[a, d]] = 1
    return A
def eig(A):
    n = len(A); k = A[0].sum(); A2 = A @ A
    lam = A2[0][np.flatnonzero(A[0])[0]]; mu = A2[0][[j for j in range(1, n) if not A[0, j]][0]]
    d = (lam - mu) ** 2 + 4 * (k - mu); q = int(round(d ** .5)); return k, lam, mu, (lam - mu + q) // 2, (lam - mu - q) // 2
def solve(A, limit, workers):
    n = len(A); k, lam, mu, r, s = eig(A); rho, sig = r + 2, s + 2
    assert k == -rho * sig, (k, rho, sig)
    a1, a0 = rho + sig, -rho * sig             # X^2 = a1 X + a0 I
    md = cp_model.CpModel(); E = {}
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j]: E[i, j] = md.NewBoolVar('')     # 1 = negative sign
    e = lambda i, j: E[min(i, j), max(i, j)]
    # switching: fix a BFS spanning tree to positive
    seen = {0}; frontier = [0]
    while frontier:
        nxt = []
        for u in frontier:
            for w in np.flatnonzero(A[u]):
                w = int(w)
                if w not in seen: seen.add(w); nxt.append(w); md.Add(e(u, w) == 0)
        frontier = nxt
    par = {}
    def P(i, w, j):
        key = (min(i, j), max(i, j), w)
        if key not in par:
            b = md.NewBoolVar(''); md.AddBoolXOr([e(i, w), e(w, j), b.Not()]); par[key] = b   # b = e_iw xor e_wj
        return par[key]
    for i in range(n):
        for j in range(i + 1, n):
            W = [int(w) for w in np.flatnonzero(A[i] & A[j])]
            # sum over w of X_iw X_wj = sum (1 - 2 p_w) must equal a1 * X_ij (X_ij = 1-2e or 0)
            ps = [P(i, w, j) for w in W]
            if A[i, j]: md.Add(len(W) - 2 * sum(ps) == a1 * (1 - 2 * e(i, j)))
            else: md.Add(len(W) - 2 * sum(ps) == 0)
    sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = limit; sol.parameters.num_workers = workers
    st = sol.Solve(md); name = sol.StatusName(st); X = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        X = np.zeros((n, n), int)
        for (i, j), b in E.items(): X[i, j] = X[j, i] = 1 - 2 * sol.Value(b)
        assert (X @ X == a1 * X + a0 * np.eye(n, dtype=int)).all()
    return name, (rho, sig), X
if __name__ == '__main__':
    A = brouwer_haemers() if sys.argv[1] == 'bh81' else None
    if A is None:
        import networkx as nx; A = nx.to_numpy_array(nx.read_graph6(sys.argv[1]), dtype=int)
    print('base', eig(A), flush=True)
    t0 = time.time(); name, sp, X = solve(A, float(sys.argv[2]), int(sys.argv[3]))
    print(sys.argv[1], 'tile spectrum', sp, name, f'{time.time()-t0:.0f}s', flush=True)
    if X is not None: np.save(sys.argv[1].split('/')[-1] + '_signing.npy', X)
