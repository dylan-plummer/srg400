"""Independent subgroup count using K5 orientations and quotienting first by translations."""
from collections import Counter
from functools import cache
from itertools import combinations,permutations
from pathlib import Path
from family import SEED,construct,tile_data,add_basis
from full_count import cover_group,free_tile_switches
from symmetric_census import code,normal


@cache
def local_fixed(perm,signs):
    # Solve (P+I)z=signs modulo each invariant simplex code. This counts
    # fixed local cosets independently of full_count's 480-state lookup.
    def move(w):
        return sum(1<<j for i,j in enumerate(perm) if w>>i&1)
    total=0
    for simplex in simplex_codes():
        if any(move(w) not in simplex for w in simplex):
            continue
        space={}
        for w in simplex:
            add_basis(space,w)
        for i,j in enumerate(perm):
            add_basis(space,(1<<i)^(1<<j))
        if not normal(signs,space):
            total+=1<<(7-len(space))
    return total


@cache
def simplex_codes():
    return {code(labels) for labels in permutations(range(1,8))}


def oriented_heads(tiles):
    supports={frozenset((p,q) for p,q,s in tile):t for t,tile in enumerate(tiles)}
    edges=list(combinations(range(5),2))
    heads={}
    for line in (SEED/'starters.txt').read_text().splitlines():
        seed=list(map(int,line.split()))
        for delta in construct.TRANSLATIONS:
            pts=[construct.shift(p,delta) for p in seed]
            support=frozenset(tuple(sorted((pts[i],pts[j]))) for i,j in combinations(range(10),2)
                              if set(edges[i])&set(edges[j]))
            t=supports[support]
            for e,(i,j) in enumerate(edges):
                head=j if j-i in (1,2) else i
                head_star=frozenset(pts[k] for k,edge in enumerate(edges) if head in edge)
                key=(pts[e],t)
                assert heads.setdefault(key,head_star)==head_star
    assert len(heads)==350
    return heads


def all_subgroups(group):
    index={g:i for i,g in enumerate(group)}
    mult=[[index[tuple(g[x] for x in h)] for h in group] for g in group]
    one=index[tuple(range(85))]
    subs={frozenset([one])}
    todo=list(subs)
    for s in todo:
        for x in range(50):
            if x in s:
                continue
            h=set(s)
            queue=list(s)
            for a in queue:
                for y in [x,*s]:
                    k=mult[a][y]
                    if k not in h:
                        queue.append(k)
                        h.add(k)
            h=frozenset(h)
            if h not in subs:
                subs.add(h)
                todo.append(h)
    return one,mult,subs


def odd_orbits(elements,points):
    return sorted({tuple(sorted({g[i] for g in elements})) for i in points})


def fixed_subgroup(sub,group,one,mult,flags,heads):
    # The odd Sylow subgroup acts by translations. An arbitrary representative
    # fixed up to tile reversals can be made invariant by odd-order averaging.
    order2=[i for i in sub if i!=one and mult[i][i]==one]
    odd=[]
    for i in sub:
        square=mult[i][i]
        fourth=mult[square][square]
        if mult[fourth][i]==one:
            odd.append(group[i])
    points=odd_orbits(odd,range(50))
    torb=odd_orbits(odd,range(50,85))
    if not order2:
        assert len(sub)==len(odd)
        return 480**len(points)//(1<<len(torb))
    assert len(sub)==2*len(odd)
    g=group[order2[0]]
    assert all(g[g[i]]==i for i in range(85))
    fixed_points=[o for o in points if g[o[0]] in o]
    fixed_tiles=[o for o in torb if g[o[0]] in o]
    ft={v-50:j for j,o in enumerate(fixed_tiles) for v in o}
    ptpairs=(len(points)-len(fixed_points))//2
    tilecycles=(len(torb)+len(fixed_tiles))//2
    local=[]
    for o in fixed_points:
        p=o[0]
        move=[None]*7
        signs=0
        effect=[0]*len(fixed_tiles)
        for (v,t),c in flags.items():
            if v!=p:
                continue
            t2=g[50+t]-50
            dest=flags[g[p],t2]-1
            move[c-1]=dest
            transported=frozenset(g[i] for i in heads[p,t])
            signs |= (transported!=heads[g[p],t2])<<dest
            if t in ft:
                effect[ft[t]] ^= 1<<dest
        local.append((tuple(move),signs,effect))
    total=0
    for mask in range(1<<len(fixed_tiles)):
        count=1
        for perm,s,e in local:
            for i,flip in enumerate(e):
                if mask>>i&1:
                    s ^= flip
            count*=local_fixed(perm,s)
        total+=count
    f,rem=divmod(480**ptpairs*total,1<<tilecycles)
    assert rem==0
    return f


def main():
    tiles,flags=tile_data()
    free_tile_switches(flags)
    group=cover_group(tiles)
    one,mult,subs=all_subgroups(group)
    heads=oriented_heads(tiles)
    fixes={s:fixed_subgroup(s,group,one,mult,flags,heads) for s in subs}
    exact={}
    orders=Counter()
    for s in sorted(subs,key=len,reverse=True):
        exact[s]=fixes[s]-sum(st for t,st in exact.items() if s<t)
        assert exact[s]>=0
        orders[len(s)]+=exact[s]
    assert sum(orders.values())==fixes[frozenset([one])]
    counts={}
    for n,points in sorted(orders.items()):
        counts[n],rem=divmod(points,50//n)
        assert rem==0
        print('full automorphism order',n,'unlabelled graphs',counts[n])
    total=sum(counts.values())
    print('subgroups of the cover group',len(subs))
    print('exact unlabelled classes',total)
    print('classes with nontrivial automorphism group',total-counts[1])
    assert counts[25]+counts[50]==14656 and counts[50]==512
    direct=next(line.split()[-1] for line in Path(__file__).with_name('full_count.out').read_text().splitlines()
                if line.startswith('exact unlabelled classes'))
    assert str(total)==direct


if __name__=='__main__':
    main()
