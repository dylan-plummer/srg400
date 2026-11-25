"""Standalone constructor. Write python from_starters.py 10 171 --symmetric > example.g6."""
import argparse
from itertools import combinations, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def offset(p, delta):
    sheet,r=divmod(p,25)
    x,y=divmod(r,5)
    a,b=delta
    return 25*sheet+5*((x+a)%5)+(y+b)%5


def unsigned(tile):
    return tuple((a,b) for a,b,s in tile)


def develop(starters):
    symbols=list(combinations(range(5),2))
    oriented=[(a,b) if (b-a)%5 in (1,2) else (b,a) for a,b in symbols]
    result=set()
    for row in starters:
        for delta in product(range(5),repeat=2):
            labels=[offset(p,delta) for p in row]
            signed=[]
            for i,j in combinations(range(10),2):
                same=set(oriented[i])&set(oriented[j])
                if same:
                    v,=same
                    negative=(oriented[i].index(v)==oriented[j].index(v))
                    a,b=sorted((labels[i],labels[j]))
                    signed.append((a,b,(-1 if negative else 1)))
            result.add(tuple(sorted(signed)))
    assert len(result)==35
    return sorted(result,key=unsigned)


def partition(support):
    def visit(free):
        if not free:
            yield []
            return
        p=min(free)
        for i,s in enumerate(support):
            if p in s and s<=free:
                for tail in visit(free-s):
                    yield [i]+tail
    answer,=visit(set(range(50)))
    return set(answer)


def build(starters,degree,decoration):
    tiles=develop(starters)
    support=[{x for a,b,s in tile for x in (a,b)} for tile in tiles]
    special=partition(support)
    index={s:i for i,s in enumerate(map(frozenset,support))}
    assert len(index)==35
    flags={}
    for p in range(50):
        sheet,loc=divmod(p,25)
        delta=divmod(loc,5)
        origins=[s for s in support if 25*sheet in s]
        labels=decoration['Q' if sheet else 'P']
        for j,s in enumerate(origins):
            i=index[frozenset(offset(v,delta) for v in s)]
            flags[p,i]=(labels[j],(decoration['negative_flags']>>(7*sheet+j))&1)
    rows=[0]*400
    def link(x,y):
        rows[x]|=1<<y
        rows[y]|=1<<x
    for i,tile in enumerate(tiles):
        reverse=degree==171 and i not in special
        for p,q,s in tile:
            c,alpha=flags[p,i]
            d,beta=flags[q,i]
            for x,y in product(range(8),repeat=2):
                parity=((x&c).bit_count()+(y&d).bit_count()+alpha+beta+reverse)%2
                if parity==int(s<0):
                    link(8*p+x,8*q+y)
        if degree==171 and i in special:
            for p in support[i]:
                c,alpha=flags[p,i]
                for x,y in combinations(range(8),2):
                    if ((x^y)&c).bit_count()%2==0:
                        link(8*p+x,8*p+y)
    assert all(x.bit_count()==degree for x in rows)
    return rows


def graph6(rows):
    data=bytearray(b'~?EO')
    assert len(rows)==400
    value=0
    bits=0
    for y in range(400):
        for x in range(y):
            value=value*2+(rows[y]>>x&1)
            bits+=1
            if bits==6:
                data.append(value+63)
                value=bits=0
    if bits:
        data.append((value<<(6-bits))+63)
    return bytes(data)+b'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cover',type=int,choices=range(12))
    parser.add_argument('degree',type=int,choices=[168,171])
    parser.add_argument('--symmetric',action='store_true')
    args=parser.parse_args()
    starters=json.loads((HERE/'representatives.json').read_text())[args.cover]
    decoration=(json.loads((HERE/'equivariant.json').read_text())[args.cover] if args.symmetric else
                {'P':range(1,8),'Q':range(1,8),'negative_flags':0})
    sys.stdout.buffer.write(graph6(build(starters,args.degree,decoration)))


if __name__=='__main__':
    main()
