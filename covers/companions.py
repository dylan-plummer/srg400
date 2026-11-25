"""Use the intrinsic parallel class of five tiles to construct the degree-171 companions."""
import contextlib
from catalogue import *
from construct import translate
from full_count import free_tile_switches, fixed_under
from family import construct as factory

def signed_data(rows):
    _,original=build(rows)
    tiles=sorted(original,key=lambda t:tuple((a,b) for a,b,_ in t))
    flags={};index={t:j for j,t in enumerate(tiles)}
    for o in (0,25):
        inc=[t for t in tiles if any(o in (a,b) for a,b,_ in t)]
        for delta in TRANSLATIONS:
            for ch,t in enumerate(inc,1):
                flags[shift(o,delta),index[translate(t,delta)]]=ch
    return tiles,flags

def parallel(tiles):
    supports=[sum(1<<p for p in {p for a,b,_ in t for p in (a,b)}) for t in tiles]
    def rec(rem,part):
        if not rem:yield part;return
        first=rem&-rem
        for i,m in enumerate(supports):
            if m&first and rem&m==m:yield from rec(rem^m,part+(i,))
    out=list(rec((1<<50)-1,()));assert len(out)==1
    return set(out[0])

def cover_action(rows,tiles):
    systems=develop(rows);tg=cover_graph(systems);element_gens,size,exp,_,_=pynauty.autgrp(tg)
    base_index={tuple(sorted(tuple(sorted(s)) for s in stars)):i for i,stars in enumerate(systems)}
    supports=[frozenset(p for a,b,_ in t for p in (a,b)) for t in tiles]
    assert len(set(supports))==35
    by_support={v:i for i,v in enumerate(supports)}
    gens=[]
    for perm in element_gens:
        pts=tuple(perm[:50])
        tgts=[50+by_support[frozenset(pts[p] for p in s)] for s in supports]
        gens.append(pts+tuple(tgts))
    els={tuple(range(85))};queue=list(els)
    for e in queue:
        for g in gens:
            h=tuple(g[i] for i in e)
            if h not in els:els.add(h);queue.append(h)
    assert len(els)==round(size*10**exp)
    return els

def companion(adj,tiles,flags):
    pos=parallel(tiles)
    ports=[[] for _ in range(50)];own={}
    result=list(adj)
    for j,t in enumerate(tiles):
        for a,b,_ in t:
            if j not in pos:
                for base,other in ((a,b),(b,a)):
                    for x in range(8):result[8*base+x]^=255<<(8*other)
        for p in {p for a,b,_ in t for p in (a,b)}:
            mask=sum(1<<(q if p!=q else v) for v,q,_ in t if p in (v,q))
            if j in pos:
                own[p]=j;mask|=1<<p
            ports[p].append(mask)
    for p in range(50):
        ch=flags[p,own[p]]
        for a,b in combinations(range(8),2):
            if (ch&(a^b)).bit_count()%2==0:
                result[8*p+a]|=1<<(8*p+b);result[8*p+b]|=1<<(8*p+a)
    for x,r in enumerate(result):
        assert r.bit_count()==171
        for y in range(x):assert (r&result[y]).bit_count()==(74 if (r>>y)&1 else 72)
    with contextlib.redirect_stdout(sys.stderr):recover(ports)
    return result

def main():
    allrows=json.loads((here/'representatives.json').read_text());rec=[]
    for cls,rows in enumerate(allrows):
        tiles,flags=signed_data(rows)
        adj,_=build(rows);b=companion(adj,tiles,flags)
        blob=graph6(b);(here/'graphs'/f'cover_{cls:02d}_400_171.g6').write_bytes(blob)
        with contextlib.redirect_stdout(sys.stderr):free_tile_switches(flags)
        group=cover_action(rows,tiles)
        # Fixed-point counts are constant on a conjugacy class. Evaluate just
        # one representative, while keeping every conjugacy-class size.
        inverses={g:tuple(g.index(i) for i in range(85)) for g in group}
        unseen=set(group);terms=Counter()
        while unseen:
            g=min(unseen)
            orbit={tuple(h[g[j]] for j in inverses[h]) for h in group}
            unseen-=orbit;terms[fixed_under(tiles,flags,g)]+=len(orbit)
        size=sum(k[1]*v for k,v in terms.items())//len(group)
        assert size*len(group)==sum(k[1]*v for k,v in terms.items())
        data={'cover':cls,'parallel':sorted(parallel(tiles)),'cover_order':len(group),'classes':size,
              'burnside':[[o,f,n] for (o,f),n in sorted(terms.items())]}
        rec.append(data);print(cls,len(group),size,flush=True)
    (here/'counts.json').write_text(json.dumps(rec,indent=1)+'\n')
    print('all',sum(x['classes'] for x in rec),'outside the old cover',sum(x['classes'] for x in rec[1:]))
if __name__=='__main__':main()
