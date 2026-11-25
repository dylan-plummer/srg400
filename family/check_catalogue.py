"""Cross-check catalogue representatives as ordinary 400-vertex graphs."""
from collections import Counter
import csv
from pathlib import Path
import random

import igraph as ig
import networkx as nx
import numpy as np
import pynauty

from family import HERE,realise,tile_data


def main():
    with (HERE/'catalogue.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==14656
    other=set((HERE/'translation_signs.g6').read_bytes().splitlines())
    selected={int(next(r['class'] for r in rows if r['group_order']==k)) for k in ('25','50')}
    selected.update(random.Random(400).sample(range(len(rows)),16))
    tiles,flags=tile_data()
    known_groups=(HERE/'group_samples.g6').read_bytes().splitlines()
    canonical=set()
    orders=Counter()
    for i in sorted(selected):
        r=rows[i];expected=int(r['group_order'])
        labels=tuple(tuple(map(int,r[k])) for k in ('P_characters','Q_characters'))
        data=realise(tiles,flags,labels,int(r['negative_flags'])).strip()
        g=nx.from_graph6_bytes(data)
        a=nx.to_numpy_array(g,nodelist=range(400),dtype=np.int32)
        assert np.array_equal(a@a,96*np.eye(400)-4*a+72)
        full=ig.Graph(400,edges=list(g.edges()))
        assert full.count_automorphisms()==expected
        ng=pynauty.Graph(400,adjacency_dict={v:list(g[v]) for v in g})
        certificate=pynauty.certificate(ng)
        assert certificate not in canonical
        canonical.add(certificate)
        orders[expected]+=1
        if data in known_groups:
            assert (data in other)==(expected==25)
        print(i,expected,'checked')
    print('Distinct canonical forms and exact group orders for the catalogue sample',dict(orders))


if __name__=='__main__':
    main()
