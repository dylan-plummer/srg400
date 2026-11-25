"""Sign and dimension bounds for cliques, uniform over decorations of a fixed cover."""
from collections import Counter
from itertools import combinations
import json

from companions import here, signed_data, parallel


def bounds(tiles, special, degree):
    n=50
    adj=[0]*n
    forbidden=[[0]*n for _ in range(n)]
    ports=[[] for _ in range(n)]
    own=[0]*n
    for t,tile in enumerate(tiles):
        mult=-1 if degree==171 and t not in special else 1
        signs={(p,q):mult*s for p,q,s in tile}
        support=sorted({v for e in signs for v in e})
        for p,q in signs:
            adj[p]|=1<<q
            adj[q]|=1<<p
        for p in support:
            ns=sum(1<<(v if u==p else u) for u,v in signs if p in (u,v))
            if t in special:
                own[p]=ns
            ports[p].append(ns)
        for triple in combinations(support,3):
            a,b,c=triple
            if signs.get((a,b),0)*signs.get((a,c),0)*signs.get((b,c),0)==-1:
                for p,q,r in ((a,b,c),(b,c,a),(c,a,b)):
                    forbidden[p][q]|=1<<r
                    forbidden[q][p]|=1<<r
    maximal=Counter()
    def visit(chosen, word, available):
        if chosen:
            if degree==168:
                bound=len(chosen)
            else:
                # The internal half-space of the distinguished character always
                # contains the selected vertices. Each external tile meeting the
                # projected clique supplies another affine hyperplane constraint.
                bound=0
                for p in chosen:
                    chars=sum(bool(ns&word) or ns==own[p] for ns in ports[p])
                    bound+=8 >> chars.bit_length()
            maximal[bound]+=1
        while available:
            bit=available&-available
            available^=bit
            p=bit.bit_length()-1
            follow=available&adj[p]
            for q in chosen:
                follow&=~forbidden[p][q]
            visit(chosen+[p],word|bit,follow)
    visit([],0,(1<<n)-1)
    best=max(maximal)
    return best,maximal[best],sum(maximal.values())


if __name__=='__main__':
    for cls,rows in enumerate(json.loads((here/'representatives.json').read_text())):
        ts,fs=signed_data(rows)
        part=parallel(ts)
        print(cls,bounds(ts,part,168),bounds(ts,part,171),flush=True)
