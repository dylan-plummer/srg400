"""Verify the graph6 exports independently, including canonical isomorphism checks."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pynauty

HERE = Path(__file__).resolve().parent
SEED = HERE.parent/'seed'


def main():
    variants = (HERE/'translation_signs.g6').read_bytes().splitlines()
    assert len(variants) == 32
    assert variants[0] == (HERE.parent/'graphs'/'srg_400_168_68_72.g6').read_bytes().strip()
    swap = (HERE.parent/'graphs'/'srg_400_168_68_72_swap.g6').read_bytes().strip()
    assert swap not in variants
    certificates = set()
    orders = Counter()
    for i,data in enumerate(variants+[swap]):
        g = nx.from_graph6_bytes(data)
        a = nx.to_numpy_array(g,nodelist=range(400),dtype=np.int32)
        assert np.array_equal(a@a,96*np.eye(400)-4*a+72)
        assert np.all(a.sum(0)==168) and np.all(a.diagonal()==0)
        ng = pynauty.Graph(400,adjacency_dict={v:list(g[v]) for v in g})
        cert = pynauty.certificate(ng)
        assert cert not in certificates
        certificates.add(cert)
        _,size,exp,_,_ = pynauty.autgrp(ng)
        order = round(size*10**exp)
        orders[order] += 1
        # Recover the planes just from the exported adjacency rows.
        rows = [sum(1<<v for v in g[u]) for u in range(400)]
        pairs = {}
        planes = set()
        for u,v in combinations(range(400),2):
            key = rows[u]^rows[v]
            for i0,j in pairs.get(key,()):
                planes.add(tuple(sorted((i0,j,u,v))))
            pairs.setdefault(key,[]).append((u,v))
        assert len(planes)==700
        assert all(len({u//8 for u in plane})==1 for plane in planes)
        print(i if i<32 else 'previous-swap',order,sha256(cert).hexdigest(),flush=True)
    print('Full canonical forms all distinct. Group orders',dict(orders))


if __name__ == '__main__':
    main()
