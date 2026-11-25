"""Second exhaustive census. No quotient by the normalizer of translations."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import json
import subprocess
import tempfile

import pynauty

HERE = Path(__file__).resolve().parent
FIVE = range(5)
EDGES = list(combinations(FIVE, 2))
TOUCH = [(i,j) for i,j in combinations(range(10),2) if set(EDGES[i]) & set(EDGES[j])]
DELTAS = list(product(FIVE, repeat=2))


def move(x, vector):
    sheet, within = divmod(x, 25)
    a,b = divmod(within,5)
    return 25*sheet + 5*((a+vector[0])%5) + (b+vector[1])%5


def orbit(starter):
    translates = [[move(p,d) for p in starter] for d in DELTAS]
    return frozenset(tuple(sorted(tuple(sorted((f[i],f[j]))) for i,j in TOUCH)) for f in translates)


def main():
    allowed = set()
    for p,q in combinations(range(50),2):
        side_p,side_q = p//25,q//25
        x,y = divmod(p%25,5)
        a,b = divmod(q%25,5)
        if side_p == side_q:
            forbidden = x == a and (b-y)%5 in ((2,3) if side_p else (1,4))
        else:
            forbidden = (b-y+3*(x-a)**2)%5 == 0
        if not forbidden:
            allowed.add((p,q))
    remain = allowed.copy()
    classes = []
    while remain:
        p,q = min(remain)
        equivalence = {tuple(sorted((move(p,d),move(q,d)))) for d in DELTAS}
        assert len(equivalence)==25 and equivalence <= remain
        classes.append(equivalence)
        remain -= equivalence
    edge_id = {e:i for i,s in enumerate(classes) for e in s}
    full = (1<<len(classes))-1
    assert len(classes)==42
    shorts = defaultdict(set)
    for v in ((1,0),(1,1),(1,2),(1,3),(1,4),(0,1)):
        for sheet0,sheet1,offset in product((0,25),(0,25),range(25)):
            row = {}
            for step, start in ((1,sheet0),(2,sheet1+offset)):
                for k in FIVE:
                    row[tuple(sorted((k,(k+step)%5)))] = move(start,(k*v[0],k*v[1]))
            tile = tuple(row[e] for e in EDGES)
            if len(set(tile)) != 10:
                continue
            edges = [tuple(sorted((tile[i],tile[j]))) for i,j in TOUCH]
            if any(e not in allowed for e in edges):
                continue
            hist = Counter(edge_id[e] for e in edges)
            if sorted(hist.values()) != [5]*6:
                continue
            mask = sum(1<<t for t in hist)
            shorts[mask].add(orbit(tile))
    options = list(shorts)
    def pack(free):
        if free==0:
            return 1
        pivot=free&-free
        return sum(pack(free^t) for t in options if t&pivot and free&t==t)
    allshort=pack(full)
    print('base edges, difference classes, short orbits, short masks, all-short covers',
          len(allowed),len(classes),sum(map(len,shorts.values())),len(shorts),allshort,flush=True)
    assert allshort == 0
    patterns = defaultdict(list)
    seen = set()
    for mask in range(1024):
        if mask in seen:
            continue
        marked = [EDGES[i] for i in range(10) if mask >> i & 1]
        for perm in permutations(FIVE):
            seen.add(sum(1 << EDGES.index(tuple(sorted((perm[a],perm[b])))) for a,b in marked))
        hist = Counter(((mask>>i)&1)+((mask>>j)&1) for i,j in TOUCH)
        patterns[tuple(hist[j] for j in range(3))].append(mask)
    pairs = defaultdict(list)
    for m,n in combinations(options,2):
        if not m & n:
            pairs[full^m^n].append((m,n))
    def profile(mask):
        hist=Counter()
        for j,equivalence in enumerate(classes):
            if mask >> j & 1:
                p,q=next(iter(equivalence))
                hist[p//25+q//25]+=1
        return tuple(hist[i] for i in range(3))
    jobs=[(m,t) for m in pairs for t in patterns[profile(m)]]
    rows = [' '.join(str(edge_id.get(tuple(sorted((p,q))),-1)) for q in range(50)) for p in range(50)]
    data='\n'.join(rows + [f'{m} {t}' for m,t in jobs])+'\n'
    with tempfile.TemporaryDirectory() as tmp:
        exe=Path(tmp)/'check'
        subprocess.run(['g++','-O3','-std=c++20',str(HERE/'census_check.cpp'),'-o',str(exe)],check=True)
        result=subprocess.run([str(exe)],input=data,text=True,capture_output=True,check=True)
    print('masks',len({m for m,t in jobs}),result.stderr.strip(),flush=True)
    covers = set()
    for line in result.stdout.splitlines():
        m,pat,*row=map(int,line.split())
        long=orbit(row)
        assert len(long)==25
        for a,b in pairs[m]:
            for left,right in product(shorts[a],shorts[b]):
                c=left|right|long
                if c in covers:
                    continue
                assert len(c)==35 and set().union(*c)==allowed
                covers.add(c)
    print('distinct decompositions on the fixed labelled host',len(covers),flush=True)
    # Canonical graph on 1135 vertices, with one edge node per host edge.
    def certificate(cover):
        graph = {i:[] for i in range(1135)}
        for s,tile in enumerate(sorted(cover)):
            for j,(p,q) in enumerate(tile,85+30*s):
                for other in (p,q,50+s):
                    graph[j].append(other)
                    graph[other].append(j)
        colour=[set(range(50)),set(range(50,85)),set(range(85,1135))]
        g=pynauty.Graph(1135,adjacency_dict=graph,vertex_coloring=colour)
        return sha256(pynauty.certificate(g)).hexdigest()
    original=json.loads((HERE/'representatives.json').read_text())
    reference=[frozenset().union(*(orbit(row) for row in rep)) for rep in original]
    assert all(c in covers for c in reference)
    keys=[certificate(c) for c in reference]
    assert len(set(keys))==12
    census=Counter(keys.index(certificate(c)) for c in covers)
    assert len(census)==12
    print('labelled covers per catalogue class',sorted(census.items()),flush=True)

if __name__ == '__main__':
    main()
