"""Cross-check the other agent's catalogue against our unsigned translation-invariant count (888)."""
import os, sys, pickle, csv, collections
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'family'))
from family import tile_data, realise
import networkx as nx
from lib400 import certificate, autorder
certs = pickle.load(open('ti_certs.pkl', 'rb'))
tiles, flags = tile_data()
rows = list(csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'family', 'catalogue.tsv')), delimiter='\t'))
hit = collections.Counter(); orders = collections.Counter(); inmine = 0
sample = [r for r in rows if r['negative_flags'] == '0']
for r in sample:
    labels = (tuple(map(int, r['P_characters'])), tuple(map(int, r['Q_characters'])))
    g6 = realise(tiles, flags, labels, 0)
    A = nx.to_numpy_array(nx.from_graph6_bytes(g6.strip()), dtype=np.int8)
    c = certificate(A)
    o, _ = autorder(A)
    assert o == int(r['group_order'])
    inmine += c in certs; orders[o] += 1
print('catalogue rows with zero signs', len(sample), 'aut orders', dict(orders), 'found among our 888 unsigned classes', inmine)
mine = collections.Counter()
for c, v in certs.items():
    from lib400 import graph, origin_assign
    mine[autorder(graph(*origin_assign(*v[0])))[0]] += 1
print('aut orders of our 888 unsigned classes', dict(mine))
