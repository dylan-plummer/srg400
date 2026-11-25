"""Classify the unsigned covers, certify extraction from any graph on a cover, and export examples."""
from collections import Counter
import contextlib,csv,hashlib,json
from pathlib import Path
from short_orbits import *
import pynauty
from construct import build,graph6
sys.path.insert(0,str(here.parent/'srg400_family_20261124'))
from structure import recover

def develop(starters):
    tiles=set()
    for row in starters:
        stars=[tuple(row[i] for i,(a,b) in enumerate(EDGES) if s in (a,b)) for s in range(5)]
        orbit={tuple(sorted(tuple(sorted(shift(p,t) for p in st)) for st in stars)) for t in TRANSLATIONS}
        assert not tiles&orbit
        tiles|=orbit
    assert len(tiles)==35
    return sorted(tiles)
def cover_graph(tiles):
    graph={i:set() for i in range(260)}
    def link(i,j):graph[i].add(j);graph[j].add(i)
    for j,tile in enumerate(tiles):
        for h,star in enumerate(tile):
            z=85+5*j+h;link(z,50+j)
            for p in star:link(z,p)
    return pynauty.Graph(260,adjacency_dict=graph,vertex_coloring=[set(range(50)),set(range(50,85)),set(range(85,260))])
def make_ports(tiles):
    ports=[[] for i in range(50)];edges=set();inside=[]
    for tile in tiles:
        support=sorted(set().union(*map(set,tile)))
        assert len(support)==10
        inside.append(sum(not in_host(difference(a,b)) for a,b in combinations(support,2)))
        for p in support:
            neighbours=set(x for star in tile if p in star for x in star if x!=p)
            assert len(neighbours)==6
            ports[p].append(sum(1<<q for q in neighbours))
            for q in neighbours:
                if q>p:
                    assert (p,q) not in edges and in_host(difference(p,q))
                    edges.add((p,q))
    assert edges=={(p,q) for p,q in combinations(range(50),2) if in_host(difference(p,q))}
    return ports,Counter(inside)
def main():
    short={int(k):v for k,v in json.loads((here/'short_orbits.json').read_text()).items()}
    joins={int(k):v for k,v in json.loads((here/'short_pairs.json').read_text()).items()}
    old=[[int(x) for x in line.split()] for line in (here.parent/'srg400_construction_20261124'/'starters.txt').read_text().splitlines()]
    known={};summ=[];start=[];graphs=here/'graphs';graphs.mkdir(exist_ok=True)
    def add(rows):
        tiles=develop(rows);g=cover_graph(tiles);key=pynauty.certificate(g)
        if key in known:return
        index=len(known);known[key]=index
        gen,size,ex,orbits,_=pynauty.autgrp(g);order=round(size*10**ex)
        ports,inside=make_ports(tiles)
        with contextlib.redirect_stdout(sys.stderr):recover(ports)
        adj,_=build(rows);blob=graph6(adj)
        name=f'cover_{index:02d}_400_168.g6';(graphs/name).write_bytes(blob)
        sig=dict(sorted(inside.items()));start.append(rows)
        data={'class':index,'group':order,'point_orbits':sorted(Counter(orbits[:50]).values()),'induced_HS_edges':sig,'sha256':hashlib.sha256(blob).hexdigest()}
        summ.append(data);print(data,flush=True)
    add(old)
    touched=set()
    for line in (here/'embeddings.txt').read_text().splitlines():
        m,pattern,*e=map(int,line.split())
        for a,b in joins[m]:
            for f in short[a]:
                for g in short[b]:
                    rows=[f,g,e];sig=pynauty.certificate(cover_graph(develop(rows)))
                    add(rows);touched.add(known[sig])
    assert touched==set(range(len(known)))
    (here/'representatives.json').write_text(json.dumps(start,indent=1)+'\n')
    (here/'catalogue.json').write_text(json.dumps(summ,indent=1)+'\n')
    print('unsigned isomorphism classes',len(known))
if __name__=='__main__':main()
