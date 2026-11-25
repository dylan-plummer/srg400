"""Exact count of translation-invariant members (unsigned): 7!^2 data / GL(3,2) per origin = 30^2 classes."""
import itertools, sys, time
import numpy as np
from lib400 import *
# GL(3,2)-orbit representatives of bijections {7 tiles} -> {1..7}: fix images of a basis... just canonicalize.
GL = []
for cols in itertools.permutations(range(1, 8), 3):
    M = lambda c: (cols[0] if c & 1 else 0) ^ (cols[1] if c & 2 else 0) ^ (cols[2] if c & 4 else 0)
    img = [M(c) for c in range(1, 8)]
    if len(set(img)) == 7 and 0 not in img: GL.append({c: M(c) for c in range(1, 8)})
assert len(GL) == 168
# characters transform by the dual (transpose-inverse) action; orbits are the same as for the primal action on 1..7
reps, seen = [], set()
for perm in itertools.permutations(range(1, 8)):
    if perm in seen: continue
    orb = {tuple(g[c] for c in perm) for g in GL}
    seen |= orb; reps.append(perm)
print('classes per origin', len(reps))
certs = {}
t0 = time.time()
for i, p0 in enumerate(reps):
    for p25 in reps:
        asg, sg = origin_assign(p0, p25)
        A = graph(asg, sg)
        c = certificate(A)
        certs.setdefault(c, []).append((p0, p25))
print('translation-invariant (unsigned) members', len(reps) ** 2, 'non-isomorphic', len(certs), f'{time.time()-t0:.0f}s')
import collections
print('class sizes', collections.Counter(len(v) for v in certs.values()))
import pickle; pickle.dump(certs, open('ti_certs.pkl', 'wb'))
# all srg? check a few + aut orders of representatives
for c, v in list(certs.items())[:5]:
    A = graph(*origin_assign(*v[0])); print(v[0], is_srg(A), autorder(A))
