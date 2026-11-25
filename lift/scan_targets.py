"""Arithmetic scan for the signed-tile lift (Theorem 2.1 of the note), with empty (c=0) or complete (c=1) fibres.

Base B: a row of Brouwer's table (existence from the row or its complement), or a complete multipartite graph
K_{m x n} (restricted eigenvalues 0, -n). Fibre: Hadamard order F (2 or 4u). Tiles: signed graphs with spectrum
{rho, sigma}, rho > 0 > sigma, degree -rho*sigma, each point in F-1 tiles. Need
    rho = r_B + 2c,  sigma = s_B + 2c,  k_B = (F-1)(-rho sigma).
Result: srg(Fv, c(F-1) + F k_B/2) with restricted eigenvalues F rho/2 - c, F sigma/2 - c.
Prints the targets with v' <= VMAX whose row is open in the snapshot."""
import csv, math, sys
VMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1300
rows = list(csv.DictReader(open(__import__('os').path.join(__import__('os').path.dirname(__file__), '..', 'data', 'brouwer_table_snapshot.csv'))))
tab = {}
for r in rows:
    try: tab[(int(r['v']), int(r['k']), int(r['lam']), int(r['mu']))] = (r['status'], r['comments'])
    except ValueError: pass
comp = lambda v, k, l, m: (v, v - k - 1, v - 2 - 2 * k + m, v - 2 * k + l)
def status(p):
    sts = [tab[x][0] for x in (p, comp(*p)) if x in tab]
    if any(s.endswith('!') or s == '+' for s in sts): return 'exists'
    if '-' in sts: return 'does not exist'
    if '?' in sts: return 'open'
    return 'infeasible/absent'
def eig(v, k, l, m):
    d = (l - m) ** 2 + 4 * (k - m); q = math.isqrt(d)
    return ((l - m + q) // 2, (l - m - q) // 2) if q * q == d and (l - m + q) % 2 == 0 else None
def srg_from(v, k, r, s):
    m = k + r * s; l = m + r + s; return (v, k, l, m)
bases = []
for p in tab:
    e = eig(*p)
    if e and p[1] < p[0] - 1: bases.append((p, e, status(p)))
for m in range(2, 400):                    # complete multipartite K_{m x n}
    for n in range(2, 400):
        if m * n * 2 <= VMAX: bases.append(((m * n, n * (m - 1), n * (m - 2), n * (m - 1)), (0, -n), f'K_{m}x{n}'))
Fs = [2] + [4 * u for u in range(1, 40)]
out = set()
for p, (rB, sB), st in bases:
    v, kB = p[0], p[1]
    for c in (0, 1):
        rho, sig = rB + 2 * c, sB + 2 * c
        if not (rho > 0 > sig): continue
        for F in Fs:
            if F * v > VMAX or kB != (F - 1) * (-rho * sig): continue
            if F * rho % 2 or F * sig % 2: continue
            t = srg_from(F * v, c * (F - 1) + F * kB // 2, F * rho // 2 - c, F * sig // 2 - c)
            stt = status(t)
            if stt == 'open': out.add((t, p, st, c, F, rho, sig))
for t, p, st, c, F, rho, sig in sorted(out):
    print(f'{t}  <- base {p} [{st}], fibres {"complete" if c else "empty"} of size {F}, '
          f'tiles: signed spectrum {{{rho},{sig}}}, degree {-rho*sig}, {F-1} per point')
