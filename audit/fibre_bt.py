"""Exact backtracking: all 8-cocliques C in which every outside vertex has 0 or 4 neighbours ('fibre-like').
Branch on an outside vertex z whose count is strictly between 0 and 4: it needs more neighbours in C."""
import sys, random, time
import numpy as np
from lib400 import *
def all_fibre_like(A, size=8, deg=4):
    n = len(A); nb = [int(''.join('1' if A[i, j] else '0' for j in range(n))[::-1], 2) for i in range(n)]
    nbl = [np.flatnonzero(A[i]).tolist() for i in range(n)]
    found = set()
    def rec(C, allowed, cnt):
        # allowed: bitmask of vertices still allowed (non-adjacent to C, > min constraint handled by canonical set)
        k = len(C)
        if k == size:
            if all(c in (0, deg) for c in cnt.values()): found.add(tuple(sorted(C)))
            return
        # deficient vertices
        best, bestopts = None, None
        for z, c in cnt.items():
            if 0 < c < deg:
                if deg - c > size - k: return
                opts = nb[z] & allowed
                if bin(opts).count('1') < deg - c: return
                if best is None or bin(opts).count('1') < bin(bestopts).count('1'): best, bestopts = z, opts
        if best is None:
            # no deficient vertex: C must already be complete or every vertex has 0/4 -> pick any allowed
            return  # with C nonempty some vertex is deficient unless all at 0/4; size<8 then impossible for our use
        o = bestopts
        while o:
            b = o & -o; o ^= b; y = b.bit_length() - 1
            cnt2 = dict(cnt)
            ok = True
            for z in nbl[y]:
                cnt2[z] = cnt2.get(z, 0) + 1
                if cnt2[z] > deg: ok = False; break
            if ok: rec(C + [y], allowed & ~nb[y] & ~b, cnt2)
            allowed &= ~b      # later branches exclude y (it was the choice here)... keeps completeness? see note
    return found, rec
def fibre_like_through(A, x):
    n = len(A); nb = [int(''.join('1' if A[i, j] else '0' for j in range(n))[::-1], 2) for i in range(n)]
    nbl = [np.flatnonzero(A[i]).tolist() for i in range(n)]
    found = set()
    def rec(C, allowed, cnt):
        k = len(C)
        if k == 8:
            if all(c in (0, 4) for c in cnt.values()): found.add(tuple(sorted(C)))
            return
        best, bestopts, bc = None, None, 10**9
        for z, c in cnt.items():
            if 0 < c < 4:
                if 4 - c > 8 - k: return
                opts = nb[z] & allowed; m = bin(opts).count('1')
                if m < 4 - c: return
                if m < bc: best, bestopts, bc = z, opts, m
        if best is None: return
        # z=best needs >=1 more neighbour from opts: branch "first chosen neighbour of z is y" (exclude earlier ys)
        o = bestopts
        while o:
            b = o & -o; o ^= b; y = b.bit_length() - 1
            cnt2 = dict(cnt); ok = True
            for z in nbl[y]:
                cnt2[z] = cnt2.get(z, 0) + 1
                if cnt2[z] > 4: ok = False; break
            if ok: rec(C + [y], allowed & ~nb[y] & ~b, cnt2)
            allowed &= ~b
    cnt = {}
    for z in nbl[x]: cnt[z] = 1
    rec([x], ((1 << n) - 1) & ~nb[x] & ~(1 << x), cnt)
    return found
if __name__ == '__main__':
    seed, mode = int(sys.argv[1]), sys.argv[2]
    if mode == 'readme':
        A = graph(*origin_assign(*[tuple(range(1, 8))] * 2))
    else:
        asg, sg = random_assign(random.Random(seed), mode == 'signed'); A = graph(asg, sg)
    assert is_srg(A)
    t0 = time.time(); allsets = set(); step = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    for x in range(0, 400, step):
        allsets |= fibre_like_through(A, x)
    fibres = {tuple(range(8 * p, 8 * p + 8)) for p in range(50)}
    print(f'seed {seed} {mode}: fibre-like 8-sets through checked vertices: {len(allsets)}, all fibres: {allsets <= fibres}, '
          f'missing fibres among checked: {len({f for f in fibres if f[0] // 8 * 8 in range(0, 400, step) or True} - allsets) if step == 1 else "n/a"}, {time.time()-t0:.0f}s', flush=True)
