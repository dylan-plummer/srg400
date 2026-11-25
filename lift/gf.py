"""Small finite fields GF(p^e) by table lookup (used to build the F = 2 base graphs)."""
import itertools
import numpy as np
from sympy import factorint

class GF:
    def __init__(self, q):
        (p, e), = factorint(q).items()
        self.q, self.p, self.e = q, p, e
        if e == 1:
            self.mul = np.fromfunction(lambda a, b: (a * b) % p, (q, q), dtype=np.int64).astype(np.int64)
            self.add = np.fromfunction(lambda a, b: (a + b) % p, (q, q), dtype=np.int64).astype(np.int64)
        else:
            poly = self._irreducible(p, e)
            digits = lambda x: [(x // p**i) % p for i in range(e)]
            enc = lambda d: sum(int(c) * p**i for i, c in enumerate(d))
            self.add = np.zeros((q, q), np.int64); self.mul = np.zeros((q, q), np.int64)
            D = [digits(x) for x in range(q)]
            for a in range(q):
                for b in range(q):
                    self.add[a, b] = enc([(x + y) % p for x, y in zip(D[a], D[b])])
                    prod = [0] * (2 * e - 1)
                    for i, x in enumerate(D[a]):
                        if x:
                            for j, y in enumerate(D[b]): prod[i + j] = (prod[i + j] + x * y) % p
                    for k in range(2 * e - 2, e - 1, -1):          # reduce by monic poly of degree e
                        c = prod[k]
                        if c:
                            for i in range(e + 1): prod[k - e + i] = (prod[k - e + i] - c * poly[i]) % p
                    self.mul[a, b] = enc(prod[:e])
        self.neg = np.array([int(np.where(self.add[x] == 0)[0][0]) for x in range(q)])
        self.order = {}
        for x in range(1, q):
            y, k = x, 1
            while y != 1: y = self.mul[y, x]; k += 1
            self.order[x] = k
    @staticmethod
    def _irreducible(p, e):
        for coeffs in itertools.product(range(p), repeat=e):
            poly = list(coeffs) + [1]
            if coeffs[0] == 0: continue
            # irreducible iff no root-free factor: brute force divisibility by monic polys of degree <= e//2
            ok = True
            for d in range(1, e // 2 + 1):
                for c2 in itertools.product(range(p), repeat=d):
                    f = list(c2) + [1]
                    r = poly[:]
                    for k in range(len(r) - 1, d - 1, -1):
                        c = r[k]
                        if c:
                            for i in range(d + 1): r[k - d + i] = (r[k - d + i] - c * f[i]) % p
                    if not any(r[:d]): ok = False; break
                if not ok: break
            if ok: return poly
        raise ValueError
