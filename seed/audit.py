"""Read the exports without importing the constructor. Check both integer matrix identities."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from hashlib import sha256
import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
GRAPHS = HERE.parent / 'graphs'


def read_export(name, valency, adjacent, other):
    path = GRAPHS / name
    graph = nx.from_graph6_bytes(path.read_bytes().strip())
    a = nx.to_numpy_array(graph, nodelist=range(400), dtype=np.int32)
    assert graph.number_of_nodes() == 400 and np.array_equal(a, a.T)
    assert np.all(a.diagonal() == 0) and np.all(a.sum(1) == valency)
    squared = a @ a
    expected = (valency-other)*np.eye(400, dtype=np.int32) + (adjacent-other)*a + other
    assert np.array_equal(squared, expected)
    i,j = np.triu_indices(400,1)
    hist = Counter(zip(a[i,j].tolist(), squared[i,j].tolist()))
    print(name, 'edges', graph.number_of_edges(), 'unordered pairs', sorted(hist.items()))
    print('sha256', sha256(path.read_bytes()).hexdigest())
    return a


def cell_data(a, translations):
    coords=[(s,r,c) for s in range(2) for r in range(5) for c in range(5)]
    def hs(p,q):
        s,r,c = p
        t,u,v = q
        if t < s:
            return hs(q,p)
        if s == t:
            return r == u and (v-c)%5 in ((1,4) if s==0 else (2,3))
        return (v-c+3*(u-r)**2)%5 == 0
    h=np.array([[hs(i,j) for j in coords] for i in coords],dtype=np.int32)
    assert np.all(h.sum(1)==7) and np.array_equal(h@h,6*np.eye(50)-h+1)
    A=a.reshape(50,8,50,8)
    Q=4*(1-np.eye(50,dtype=int)-h)
    assert np.array_equal(A.sum(3),np.repeat(Q[:,None,:],8,axis=1))
    seeds=[list(map(int,s.split())) for s in (HERE/'starters.txt').read_text().splitlines()]
    print('induced Hoffman-Singleton edge counts on the three starters',
          [int(h[np.ix_(s,s)].sum()//2) for s in seeds])
    for dr,dc in translations:
        perm=[8*coords.index((s,(r+dr)%5,(c+dc)%5))+x for s,r,c in coords for x in range(8)]
        assert np.array_equal(a[np.ix_(perm,perm)],a)
    return h


def rebuild_from_incidence(a, h, swap_at_zero):
    """Reconstruct the entire matrix using signed incidence and Walsh columns."""
    edges = list(combinations(range(5), 2))
    d = np.zeros((5, 10), dtype=np.int32)
    for e, (u, v) in enumerate(edges):
        d[u, e], d[v, e] = (-1, 1) if v-u in (1, 2) else (1, -1)
    x = 2*np.eye(10, dtype=np.int32) - d.T @ d
    assert np.array_equal(x @ x, 6*np.eye(10)-x)
    edge_indices = list(zip(*np.triu(x != 0, 1).nonzero()))

    def shift(k, u, v):
        return 25*(k//25) + 5*((k//5+u) % 5) + (k+v) % 5

    def translate_key(key, u, v):
        return tuple(sorted(tuple(sorted((shift(i,u,v), shift(j,u,v)))) for i,j in key))

    developed = {}
    weights = {}
    seeds = (HERE/'starters.txt').read_text().splitlines()
    for seed in seeds:
        labels = list(map(int, seed.split()))
        for u, v in product(range(5), repeat=2):
            points = [shift(k, u, v) for k in labels]
            signs = {tuple(sorted((points[i],points[j]))): x[i,j] for i,j in edge_indices}
            key = tuple(sorted(signs))
            if key in developed:
                assert weights[key] == signs
            else:
                developed[key] = points
                weights[key] = signs
    b = 1 - np.eye(50, dtype=np.int32) - h
    assert Counter(e for key in developed for e in key) == Counter(
        (i,j) for i,j in combinations(range(50), 2) if b[i,j])
    ranks = {}
    for k in range(50):
        incident = [key for key, points in developed.items() if k in points]
        incident.sort(key=lambda key: translate_key(key, -(k//5), -(k%5)))
        assert len(incident) == 7
        if k == 0 and swap_at_zero:
            incident[0], incident[1] = incident[1], incident[0]
        ranks.update({(k,key): i for i,key in enumerate(incident,1)})
    w = np.array([[(-1)**((i&j).bit_count()) for j in range(8)] for i in range(8)], dtype=np.int32)
    assert np.array_equal(w.T @ w, 8*np.eye(8))
    twice_adjacency = np.kron(b, np.ones((8,8), dtype=np.int32))
    for key, points in developed.items():
        f = np.zeros((400,10), dtype=np.int32)
        for e, k in enumerate(points):
            f[8*k:8*k+8,e] = w[:,ranks[k,key]]
        twice_adjacency += f @ x @ f.T
    assert np.array_equal(twice_adjacency, 2*a)
    print('incidence-matrix reconstruction agrees with the export in all 160000 entries')


def main():
    for suffix in ('', '_swap'):
        a=read_export(f'srg_400_168_68_72{suffix}.g6',168,68,72)
        b=read_export(f'srg_400_231_134_132{suffix}.g6',231,134,132)
        assert np.array_equal(a+b,1-np.eye(400,dtype=int))
        h = cell_data(a, [] if suffix else [(1,0),(0,1)])
        rebuild_from_incidence(a, h, swap_at_zero=bool(suffix))
    print('all checks passed')


if __name__=='__main__':
    main()
