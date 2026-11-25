"""The signed-tile lift (Section 2 of the note), for any base graph, tiles and Hadamard fibre.

Input: a base graph B on v points; tiles = list of (points, X) with X a symmetric 0/+-1 matrix on those points, the
tiles edge-disjoint with union B; a normalized Hadamard matrix of order F; at each point a bijection from its incident
tiles to the F-1 nonconstant rows of the Hadamard matrix; optional signs. Output: the adjacency matrix of the lift on
F*v vertices, vertex (p, x) = F*p + x. Theorem 2.1: if B is strongly regular with eigenvalues k, r, s, every tile has
spectrum in {r, s} and every point lies in exactly F-1 tiles, the lift is srg with eigenvalues Fk/2, Fr/2, Fs/2.
usage: python general_lift.py    (runs the controls)
"""
import itertools, random
import numpy as np

def sylvester(m):
    H = np.ones((1, 1), int)
    for _ in range(m): H = np.block([[H, H], [H, -H]])
    return H

def lift(v, tiles, Had, assign, signs=None):
    F = len(Had); A = np.zeros((F * v, F * v), np.int8)
    for t, (pts, X) in enumerate(tiles):
        for i, p in enumerate(pts):
            for j, q in enumerate(pts):
                if X[i][j] == 0: continue
                fp = Had[assign[p, t]] * (signs[p, t] if signs else 1)
                fq = Had[assign[q, t]] * (signs[q, t] if signs else 1)
                A[F * p:F * p + F, F * q:F * q + F] = (np.outer(fp, fq) == X[i][j])
    return A

def srg_parameters(A):
    """Return (v,k,lambda,mu) if A is strongly regular (exact integer check on all pairs), else None."""
    A = A.astype(np.int64); n = len(A); k = int(A[0].sum())
    if not (A == A.T).all() or A.diagonal().any() or not (A.sum(1) == k).all(): return None
    A2 = A @ A; off = ~np.eye(n, dtype=bool)
    lam = set(A2[(A == 1) & off].tolist()); mu = set(A2[(A == 0) & off].tolist())
    if len(lam) > 1 or len(mu) > 1: return None
    return (n, k, lam.pop() if lam else 0, mu.pop() if mu else 0)

def random_assignment(v, tiles, F, rng, signs=True):
    inc = {p: [t for t, (pts, _) in enumerate(tiles) if p in pts] for p in range(v)}
    assert all(len(ts) == F - 1 for ts in inc.values()), 'every point must lie in F-1 tiles'
    asg, sg = {}, {}
    for p, ts in inc.items():
        for row, t in zip(rng.sample(range(1, F), F - 1), ts):
            asg[p, t] = row; sg[p, t] = rng.choice((1, -1)) if signs else 1
    return asg, sg

def control_fano():
    """K7 split into the 7 lines of the Fano plane (all-plus triangles, spectrum {2,-1}); F = 4: srg(28,12,6,4)."""
    lines = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]
    tiles = [(L, [[0 if i == j else 1 for j in range(3)] for i in range(3)]) for L in lines]
    return 7, tiles, sylvester(2)

def paley_hadamard(q=11):
    """Normalized Hadamard matrix of order q+1 (Paley I, q = 3 mod 4 prime)."""
    chi = lambda a: 0 if a % q == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)
    Q = np.array([[chi(j - i) for j in range(q)] for i in range(q)])
    H = np.ones((q + 1, q + 1), int); H[1:, 1:] = Q - np.eye(q, dtype=int)
    H = H * H[:, :1]; H = H * H[:1, :]          # normalize the first row and column to all-ones
    assert (H @ H.T == (q + 1) * np.eye(q + 1, dtype=int)).all()
    return H

def control_k12():
    """K12 tiled by its 66 edges (K2, spectrum {1,-1}); Hadamard order F = 12 (not a power of 2): srg(144,66,30,30)."""
    tiles = [((a, b), [[0, 1], [1, 0]]) for a, b in itertools.combinations(range(12), 2)]
    return 12, tiles, paley_hadamard(11)

def control_400():
    """The 35 signed T(5) tiles of the complement of Hoffman-Singleton; F = 8: srg(400,168,68,72)."""
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'audit'))
    from independent_rebuild import tiles as hs_tiles
    T = hs_tiles(); out = []
    for key, M in T.items():
        pts = sorted(key); out.append((pts, [[int(M[p, q]) for q in pts] for p in pts]))
    return 50, out, sylvester(3)

if __name__ == '__main__':
    rng = random.Random(2026)
    for name, ctl in (('Fano/K7, F=4', control_fano), ('K12 by edges, F=12', control_k12), ('co-Hoffman-Singleton, F=8', control_400)):
        v, tiles, Had = ctl()
        res = set()
        for trial in range(5):
            asg, sg = random_assignment(v, tiles, len(Had), rng)
            res.add(srg_parameters(lift(v, tiles, Had, asg, sg)))
        print(name, '-> parameters over 5 random assignments with random signs:', res)
