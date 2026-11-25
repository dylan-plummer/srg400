"""Independent count of base cliques that can carry a clique of some member: B-cliques S such that every triple of S
lying pairwise in one tile T has positive sign product in X_T (balanced). Necessary and, per tile, characters-free."""
import itertools, collections
from lib400 import *
H = hs(); Badj = [set(np.flatnonzero(1 - H[p])) - {p} for p in range(50)]
def ok_add(S, p):
    for q, r in itertools.combinations(S, 2):
        k = TILE_OF[q, r]
        if TILE_OF[p, q] is k and TILE_OF[p, r] is k and T[k][p, q] * T[k][p, r] * T[k][q, r] == -1: return False
    return True
cnt = collections.Counter()
def rec(S, cand):
    cnt[len(S)] += 1
    for i, p in enumerate(cand):
        if ok_add(S, p): rec(S + [p], [q for q in cand[i + 1:] if q in Badj[p]])
rec([], list(range(50)))
print('balanced base cliques by size', sorted(cnt.items()))
