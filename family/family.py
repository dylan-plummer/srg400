"""Data of the 35 tiles and explicit translation-invariant sign variants."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SEED = HERE.parent / 'seed'
sys.path.insert(0, str(SEED))
import construct


def tile_data():
    seeds = [list(map(int,l.split())) for l in (SEED/'starters.txt').read_text().splitlines()]
    _, signed_tiles = construct.build(seeds)
    tiles = sorted(signed_tiles, key=lambda t: tuple((p,q) for p,q,_ in t))
    index = {t:i for i,t in enumerate(tiles)}
    flags = {}
    for origin in (0,25):
        incident = [t for t in tiles if any(origin in (p,q) for p,q,_ in t)]
        for delta in construct.TRANSLATIONS:
            p = construct.shift(origin,delta)
            for char,t in enumerate(incident,1):
                flags[p,index[construct.translate(t,delta)]] = char
    assert len(flags) == 350
    return tiles, flags


def add_basis(basis, vector):
    while vector:
        bit = vector.bit_length()-1
        if bit not in basis:
            basis[bit] = vector
            return True
        vector ^= basis[bit]
    return False


def representatives(tiles, flags):
    # Bit 7*p+c-1 switches the flag carrying character c at point p.
    tile_switch = [0]*len(tiles)
    translations = [0]*150
    invariant_signs = [0]*14
    for (p,t),c in flags.items():
        bit = 1 << (7*p+c-1)
        tile_switch[t] |= bit
        invariant_signs[7*(p//25)+c-1] |= bit
        for j in range(3):
            if c >> j & 1:
                translations[3*p+j] |= bit
    basis = {}
    for mask in tile_switch + translations:
        add_basis(basis,mask)
    assert len(basis) == 185
    masks = [v for v in invariant_signs if add_basis(basis,v)]
    print('independent translation-invariant sign coordinates', len(masks))
    assert len(masks) == 5
    return masks


def graph(tiles, flags, flipped):
    rows = [0]*400
    for t,tile in enumerate(tiles):
        for p,q,sign in tile:
            c,d = flags[p,t],flags[q,t]
            flip = ((flipped >> (7*p+c-1)) ^ (flipped >> (7*q+d-1))) & 1
            for x in range(8):
                for y in range(8):
                    if ((x&c).bit_count()+(y&d).bit_count()+flip)%2 == (sign<0):
                        rows[8*p+x] |= 1 << (8*q+y)
                        rows[8*q+y] |= 1 << (8*p+x)
    assert all(x.bit_count()==168 for x in rows)
    counts = Counter(((rows[i] >> j)&1,(rows[i]&rows[j]).bit_count()) for i,j in combinations(range(400),2))
    assert counts == {(1,68):33600,(0,72):46200}
    return construct.graph6(rows)


def realise(tiles,flags,labels,mask):
    fg={}
    flip=0
    for (p,t),c in flags.items():
        sheet=p//25
        char=labels[sheet][c-1]
        fg[p,t]=char
        if mask>>(7*sheet+c-1)&1:
            flip |= 1<<(7*p+char-1)
    return graph(tiles,fg,flip)


def main():
    tiles,flags = tile_data()
    masks = representatives(tiles,flags)
    with (HERE/'translation_signs.g6').open('wb') as f:
        for i in range(1<<len(masks)):
            # These particular masks have disjoint supports.
            flip = sum(m for j,m in enumerate(masks) if (i>>j)&1)
            data = graph(tiles,flags,flip)
            if i == 0:
                assert data == (HERE.parent/'graphs'/'srg_400_168_68_72.g6').read_bytes()
            f.write(data)
    print('32 exact pair-count checks passed')


if __name__ == '__main__':
    main()
