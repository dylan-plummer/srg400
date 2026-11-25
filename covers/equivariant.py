"""Classify decorations fixed by the translation group on each unsigned cover."""
from companions import *
from itertools import permutations
from symmetric_census import code,act,normal,tile_masks
from family import add_basis,realise
from functools import cache

def induced(tiles,flags,perm):
    move=[-1]*14;equations=set()
    signs=[{(p,q):s for p,q,s in t} for t in tiles]
    def port(p,t):return 7*(p//25)+flags[p,t]-1
    for i,t in enumerate(tiles):
        j=perm[50+i]-50
        for a,b,s in t:
            aa,bb=perm[a],perm[b]
            for u,v in [(a,aa),(b,bb)]:
                x,y=port(u,i),port(v,j)
                assert move[x] in (-1,y);move[x]=y
            m=(1<<port(aa,j))^(1<<port(bb,j))
            equations.add((m,s!=signs[j][tuple(sorted((aa,bb)))]))
    assert sorted(move)==list(range(14))
    # Translate-invariant sign lifts exist: the tile flips on each translation
    # orbit are the residual freedom. Solve on just fourteen orbit flags.
    gamma=next(z for z in range(1<<14) if all((z&m).bit_count()%2==v for m,v in equations))
    return move,gamma

def main():
    rows0=json.loads((here/'representatives.json').read_text())
    labels={}
    for x in permutations(range(1,8)):labels.setdefault(code(x),x)
    codes=list(labels);idx={x:i for i,x in enumerate(codes)}
    records=[]
    for cover,rows in enumerate(rows0):
        tiles,flags=signed_data(rows);grp=cover_action(rows,tiles);K=[g for g in grp if g[0]==0]
        flip=tile_masks(tiles,flags);assert len(flip)==3
        actions=[]
        for g in K:
            move,gamma=induced(tiles,flags,g)
            table=[[idx[tuple(sorted(act(w,[move[i+7*s]-7*s for i in range(7)]) for w in c))] for c in codes] for s in range(2)]
            actions.append((move,gamma,table))
        basis={}
        for a,b in product(range(30),repeat=2):
            bs={};basis[a,b]=bs
            for m in flip+list(codes[a])+[w<<7 for w in codes[b]]:add_basis(bs,m)
            assert len(bs)==9
        seen=set();count=Counter();examples={}
        for a,b in product(range(30),repeat=2):
            bs=basis[a,b];masks=[0]
            for i in range(14):
                if i not in bs:masks +=[v|(1<<i) for v in masks]
            for z in masks:
                if (a,b,z) in seen:continue
                orb=set()
                for mov,gg,tab in actions:
                    aa,bb=tab[0][a],tab[1][b];zz=normal(act(z,mov)^gg,basis[aa,bb])
                    orb.add((aa,bb,zz))
                assert (a,b,z) in orb and not seen&orb
                seen|=orb;sz=len(grp)//len(orb);count[sz]+=1
                examples.setdefault(sz,(a,b,z))
        assert len(seen)==28800
        best=max(count);a,b,z=examples[best]
        blob=realise(tiles,flags,(labels[codes[a]],labels[codes[b]]),z)
        (here/'graphs'/f'sym_{cover:02d}_400_168.g6').write_bytes(blob)
        rec={'cover':cover,'counts':dict(count),'max_group':best,'P':labels[codes[a]],'Q':labels[codes[b]],'negative_flags':z}
        records.append(rec);print(cover,dict(sorted(count.items())),flush=True)
    (here/'equivariant.json').write_text(json.dumps(records,indent=1)+'\n')
if __name__=='__main__':main()
