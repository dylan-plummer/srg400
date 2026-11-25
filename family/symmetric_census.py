"""Classify decorations that repeat under the fixed translations of the three-starter cover."""
from collections import Counter
import csv
from itertools import permutations,product
from pathlib import Path
from family import add_basis,tile_data,realise
from structure import group

HERE = Path(__file__).resolve().parent


def code(labels):
    return tuple(sorted(sum(1<<i for i,c in enumerate(labels) if (c&x).bit_count()%2) for x in range(8)))


def act(mask, move):
    return sum(1<<j for i,j in enumerate(move) if mask>>i&1)


def normal(mask,basis):
    for i in sorted(basis,reverse=True):
        if (mask>>i)&1:
            mask ^= basis[i]
    return mask


def outer_action(tiles,flags):
    perm = group(tiles)
    move = [-1]*14
    def coord(p,t):
        return 7*(p//25)+flags[p,t]-1
    signs = [{(p,q):s for p,q,s in t} for t in tiles]
    constraints = set()
    for t,tile in enumerate(tiles):
        t2 = perm[50+t]-50
        for p,q,s in tile:
            p2,q2=perm[p],perm[q]
            for v,v2 in ((p,p2),(q,q2)):
                i,j=coord(v,t),coord(v2,t2)
                assert move[i] in (-1,j)
                move[i]=j
            a,b=coord(p2,t2),coord(q2,t2)
            constraints.add(((1<<a)^(1<<b),s!=signs[t2][tuple(sorted((p2,q2)))]))
    assert sorted(move)==list(range(14)) and all(move[move[i]]==i for i in range(14))
    gamma = next(z for z in range(1<<14) if all((m&z).bit_count()%2==b for m,b in constraints))
    return move,gamma


def tile_masks(tiles,flags):
    masks = []
    for t in range(len(tiles)):
        masks.append(sum(1<<i for i in {7*(p//25)+c-1 for (p,s),c in flags.items() if s==t}))
    return sorted(set(masks))


def full_kernel(codes,flags):
    # A local translation fixes all blocks iff its sign on each tile is constant.
    # Use 35 unknown tile constants, with four simplex-code constraints per base point.
    basis = {}
    for p in range(50):
        incident = sorted(((c,t) for (v,t),c in flags.items() if v==p))
        simplex = codes[p//25]
        for word in range(1,128):
            if all((w&word).bit_count()%2==0 for w in simplex):
                add_basis(basis,sum(1<<t for c,t in incident if word>>(c-1)&1))
    return 35-len(basis)


def main():
    tiles,flags=tile_data()
    move,gamma=outer_action(tiles,flags)
    tile_flips=tile_masks(tiles,flags)
    assert len(tile_flips)==3
    planes={}
    for labels in permutations(range(1,8)):
        planes.setdefault(code(labels), labels)
    assert len(planes)==30
    codes=list(planes)
    idx={z:i for i,z in enumerate(codes)}
    on_codes=[[idx[tuple(sorted(act(w,[move[7*s+i]-7*s for i in range(7)]) for w in z))]
               for z in codes] for s in range(2)]
    gauges={}
    for a,b in product(range(30),repeat=2):
        basis={}
        for m in [*tile_flips,*codes[a],*(w<<7 for w in codes[b])]:
            add_basis(basis,m)
        gauges[a,b]=basis
    groups=Counter()
    identity_classes=0
    seen=set()
    samples={}
    catalogue=[]
    burnside_fixed=0
    kernels=Counter()
    for a,b in product(range(30),repeat=2):
        basis=gauges[a,b]
        kernel=full_kernel((codes[a],codes[b]),flags)
        kernels[kernel]+=1
        other=on_codes[0][a],on_codes[1][b]
        if other==(a,b):
            # Fixed cosets solve (1+move)z = gamma modulo the gauge.
            image_basis=basis.copy()
            for i,j in enumerate(move):
                add_basis(image_basis,(1<<i)^(1<<j))
            if not normal(gamma,image_basis):
                burnside_fixed += 1<<(14-len(image_basis))
        # The complement coordinates give one vector for each sign coset.
        reps=[0]
        for i in range(14):
            if i not in basis:
                reps += [r^(1<<i) for r in reps]
        identity_classes+=len(reps)
        for z in reps:
            k=a,b,z
            if k in seen:
                continue
            z2=normal(act(z,move)^gamma,gauges[other])
            image=(*other,z2)
            assert normal(act(z2,move)^gamma,basis)==z
            seen.update((k,image))
            order=(50 if k==image else 25)*(1<<kernel)
            groups[order]+=1
            samples.setdefault(order,k)
            catalogue.append((len(catalogue),order,''.join(map(str,planes[codes[a]])),
                              ''.join(map(str,planes[codes[b]])),z))
    assert identity_classes==len(seen)
    assert (identity_classes+burnside_fixed)//2==len(catalogue)
    assert kernels=={0:900}
    with (HERE/'catalogue.tsv').open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n')
        writer.writerow(('class','group_order','P_characters','Q_characters','negative_flags'))
        writer.writerows(catalogue)
    print('classes over the identity on the base',identity_classes)
    print('fixed sign cosets from the independent rank calculation',burnside_fixed)
    print('dimensions of the fibre-translation kernel',dict(kernels))
    print('full isomorphism classes',sum(groups.values()),'automorphism orders',dict(sorted(groups.items())))
    with (HERE/'group_samples.g6').open('wb') as f:
        for order,(a,b,z) in sorted(samples.items()):
            f.write(realise(tiles,flags,(planes[codes[a]],planes[codes[b]]),z))
            print('group representative',order,'Fano indices',a,b,'signs',z)


if __name__=='__main__':
    main()
