"""Enumerate the order-five triangular tile orbits in the biaffine Hoffman--Singleton model."""
import sys,json
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations,product
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parent/'srg400_construction_20261124'))
from construct import EDGES,shift,signed_edges,translate,TRANSLATIONS
adjacent=[(i,j) for i,j in combinations(range(10),2) if set(EDGES[i])&set(EDGES[j])]
def difference(a,b):
    a,b=sorted((a,b));s=a//25+b//25
    h,k=divmod(a%25,5);i,j=divmod(b%25,5)
    x=5*((i-h)%5)+(j-k)%5
    if s!=1:x=min(x,5*((h-i)%5)+(k-j)%5)
    return s,x
def in_host(pair):
    s,x=pair;a,b=divmod(x,5)
    return ((b+3*a*a)%5!=0 if s==1 else x!=0 and (a!=0 or b not in ((1,4) if s==0 else (2,3))))
classes=sorted({difference(x,y) for x,y in combinations(range(50),2) if in_host(difference(x,y))})
assert len(classes)==42
def orbit_key(row):
    support=tuple(sorted(tuple(sorted((row[i],row[j]))) for i,j in adjacent))
    return min(tuple(sorted(tuple(sorted((shift(a,t),shift(b,t)))) for a,b in support)) for t in TRANSLATIONS)

def short_orbits():
    seen={};by_mask=defaultdict(list)
    for first,second in product((0,25),repeat=2):
        for v in TRANSLATIONS[1:]:
            for offset in range(25):
                row=[None]*10
                for k,base in [(1,first),(2,second+offset)]:
                    for i in range(5):row[EDGES.index(tuple(sorted((i,(i+k)%5))))]=shift(base,(i*v[0],i*v[1]))
                if len(set(row))!=10:continue
                hist=Counter(difference(row[i],row[j]) for i,j in adjacent)
                if len(hist)!=6 or set(hist.values())!={5} or not all(in_host(p) for p in hist):continue
                key=orbit_key(row)
                if key in seen:continue
                mask=sum(1<<classes.index(c) for c in hist)
                seen[key]=row;by_mask[mask].append(row)
    return by_mask
def main():
    by_mask=short_orbits();types=Counter()
    for m,rows in by_mask.items():
        s=sum(bool(m>>i&1) for i,p in enumerate(classes) if p[0]==1);types[s]+=len(rows)
    print('short orbits',sum(map(len,by_mask.values())),'masks',len(by_mask),'by cross multiplicity',dict(types))
    universe=(1<<42)-1
    options=[[m for m in by_mask if m>>j&1] for j in range(42)]
    covers=[]
    def rec(rem,l):
        if not rem:covers.append(l);return
        j=min((i for i in range(42) if rem>>i&1),key=lambda i:sum(m&rem==m for m in options[i]))
        for m in options[j]:
            if m&rem==m:rec(rem^m,l+[m])
    rec(universe,[])
    print('all-short orbit mask covers',len(covers))
    (here/'short_orbits.json').write_text(json.dumps(by_mask)+'\n')
    (here/'seven_short.json').write_text(json.dumps(covers)+'\n')
if __name__=='__main__':main()
