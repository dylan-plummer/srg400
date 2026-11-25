"""Count the whole unlabelled family by Burnside on local affine simplex codes."""
from collections import Counter
from functools import cache
from itertools import permutations, product

from family import construct,tile_data,add_basis
from structure import group
from symmetric_census import code,outer_action,tile_masks


def cycles(perm):
    left=set(range(len(perm)))
    result=[]
    while left:
        i=min(left)
        cycle=[]
        while i in left:
            left.remove(i)
            cycle.append(i)
            i=perm[i]
        assert i==cycle[0]
        result.append(cycle)
    return result


def free_tile_switches(flags):
    rows=[sum(1<<t for q,t in flags if q==p) for p in range(50)]
    basis={}
    for r in rows:
        add_basis(basis,r)
    free=[j for j in range(35) if j not in basis]
    good=[]
    for word in range(1<<len(free)):
        x=sum(1<<j for i,j in enumerate(free) if word>>i&1)
        for j,b in sorted(basis.items()):
            x |= ((x&b).bit_count()%2)<<j
        assert all((x&r).bit_count()%2==0 for r in rows)
        if all((x&r).bit_count() in (0,4) for r in rows):
            good.append(x)
    assert good==[0]
    print('point-tile incidence rank',len(basis),'even sets',1<<len(free),'allowed sets',len(good))


def local_states():
    codes={code(z) for z in permutations(range(1,8))}
    states={tuple(sorted(w^s for w in c)) for c in codes for s in range(128)}
    assert len(states)==480
    return sorted(states)


STATES=local_states()


@cache
def word_action(perm):
    return tuple(sum(1<<v for u,v in enumerate(perm) if w>>u&1) for w in range(128))


@cache
def local_fixed(perm,signs):
    act=word_action(perm)
    return sum(all((act[w]^signs) in state for w in state) for state in STATES)


def cover_group(tiles):
    idx={tile:i for i,tile in enumerate(tiles)}
    outer=group(tiles)
    e=list(range(85))
    full=set()
    for d in construct.TRANSLATIONS:
        trans=[construct.shift(p,d) for p in range(50)]
        trans.extend(50+idx[construct.translate(t,d)] for t in tiles)
        for g in (e,outer):
            full.add(tuple(trans[g[i]] for i in range(85)))
    assert len(full)==50
    return sorted(full)


def phases(tiles,g):
    """Pick either of the two sign lifts on each tile by propagating along its edges."""
    result={}
    for t,tile in enumerate(tiles):
        target={(a,b):s for a,b,s in tiles[g[50+t]-50]}
        constraints=[(a,b,(s!=target[tuple(sorted((g[a],g[b])))])) for a,b,s in tile]
        first=constraints[0][0]
        value={first:0}
        while len(value)<10:
            prev=len(value)
            for a,b,s in constraints:
                if a in value:
                    value[b]=value[a]^s
                elif b in value:
                    value[a]=value[b]^s
            assert len(value)>prev
        assert all(value[a]^value[b]==s for a,b,s in constraints)
        result.update({(p,t):v for p,v in value.items()})
    return result


def fixed_under(tiles,flags,g):
    import math
    base=g[:50]
    perm=tuple(t-50 for t in g[50:])
    point_cycles=cycles(base)
    tile_cycles=cycles(perm)
    order=math.lcm(*(len(c) for c in point_cycles))
    if order==1:
        return order,480**50//2**35
    gamma=phases(tiles,g)
    # Compute the d-th power of this particular lift on one flag per tile orbit.
    allowed=[]
    for orbit in tile_cycles:
        t=orbit[0]
        p=next(p for p,u in flags if u==t)
        s=0
        for _ in range(order):
            s ^= gamma[p,t]
            p,t=base[p],perm[t]
        if (order//len(orbit))%2:
            allowed.append((s,))
        elif s:
            return order,0
        else:
            allowed.append((0,1))
    incident={p:[t for c,t in sorted((c,t) for (q,t),c in flags.items() if q==p)] for p in range(50)}
    # Record the signed action after one circuit of each point orbit. Alongside
    # the reference signs, keep track of how an extra tile flip changes them.
    circuits=[]
    for orbit in point_cycles:
        dest=[]
        signs=0
        sensitivity=[]
        for j,t in enumerate(incident[orbit[0]]):
            s=0
            mask=0
            for p in orbit:
                s ^= gamma[p,t]
                mask ^= 1<<t
                t=perm[t]
            k=flags[orbit[0],t]-1
            dest.append(k)
            signs |= s<<k
            sensitivity.append((k,mask))
        circuits.append((tuple(dest),signs,sensitivity))
    lifted=0
    for parity in product(*allowed):
        z=sum(bit<<orb[0] for orb,bit in zip(tile_cycles,parity))
        choices=1
        for dest,s,sensitivity in circuits:
            sign=s^sum(((z&mask).bit_count()%2)<<k for k,mask in sensitivity)
            choices *= local_fixed(dest,sign)
            if not choices:
                break
        lifted+=choices
    # Conjugating by tile flips accounts for 2^(35 - #tile orbits) lifts
    # per parity vector. Divide their total by the free group of 2^35 flips.
    quotient,rem=divmod(lifted,1<<len(tile_cycles))
    assert rem==0
    if order==5:
        assert quotient==480**10//(1<<len(tile_cycles))
    return order,quotient


def main():
    tiles,flags=tile_data()
    free_tile_switches(flags)
    # Independent control on the two origins and the three translated tile classes.
    move,gamma=outer_action(tiles,flags)
    flip=tile_masks(tiles,flags)
    total=0
    for bits in product(range(2),repeat=3):
        z=gamma
        for m,b in zip(flip,bits):
            if b:z^=m
        total+=local_fixed(tuple(move[:7]),z&127)*local_fixed(tuple(i-7 for i in move[7:]),z>>7)
    assert total==8*512
    print('two-origin control: fixed invariant classes',total//8)
    terms=Counter(fixed_under(tiles,flags,g) for g in cover_group(tiles))
    numerator=sum(n*fixed for (order,fixed),n in terms.items())
    count,remainder=divmod(numerator,50)
    assert remainder==0
    for (order,fixed),n in sorted(terms.items()):
        print('element order',order,'elements',n,'fixed classes',fixed)
    print('exact unlabelled classes',count)
    print('excess over earlier lower bound',count-30**50*2**165//50)


if __name__=='__main__':
    main()
