"""Reduce all long-tile problems modulo the normalizer of the chosen translations."""
from short_orbits import *
from itertools import permutations
def normalizer():
    group=[]
    for a,b,c,d in product(range(5),repeat=4):
        if (a*d-b*c)%5==0:continue
        for swap in (0,1):
            # The images of the two within-sheet forbidden differences must stay forbidden.
            if b or d not in ((2,3) if swap else (1,4)):continue
            for u,v in TRANSLATIONS:
                perm=[25*(s^swap)+5*((a*x+b*y+s*u)%5)+(c*x+d*y+s*v)%5 for s in (0,1) for x,y in TRANSLATIONS]
                if any(in_host(difference(perm[0],perm[25+5*x+(-3*x*x)%5])) for x in range(5)):continue
                lift=[]
                for s,h in classes:
                    i,j=(0,h) if s==0 else ((25,25+h) if s==2 else (0,25+h))
                    lift.append(classes.index(difference(perm[i],perm[j])))
                group.append((perm,lift))
    return group
def image(m,g):return sum(1<<g[i] for i in range(len(g)) if m>>i&1)
def profile(m):return tuple(sum(1 for i,(s,_) in enumerate(classes) if s==j and m>>i&1) for j in range(3))
def patterns():
    syms=[]
    for g in permutations(range(5)):
        syms.append([EDGES.index(tuple(sorted((g[a],g[b])))) for a,b in EDGES])
    by_profile=defaultdict(list)
    for s in range(1024):
        if s!=min(image(s,g) for g in syms):continue
        h=tuple(sum(((s>>i&1)+(s>>j&1))==k for i,j in adjacent) for k in range(3))
        by_profile[h].append(s)
    return by_profile
def main():
    by_mask={int(k):v for k,v in json.loads((here/'short_orbits.json').read_text()).items()}
    full=(1<<42)-1; colors=patterns();masks=defaultdict(list)
    for a,b in combinations(by_mask,2):
        m=full^a^b
        if not a&b and profile(m) in colors:masks[m].append([a,b])
    norm=normalizer()
    remaining=set(masks);representatives=[]
    while remaining:
        m=min(remaining);orbit={image(m,g) for _,g in norm};assert orbit<=remaining
        remaining-=orbit;representatives.append(m)
    with (here/'problems.in').open('w') as f:
        for i in range(50):
            print(*[classes.index(difference(i,j)) if i!=j and in_host(difference(i,j)) else -1 for j in range(50)],file=f)
        for m in representatives:
            for pat in colors[profile(m)]:print(m,pat,file=f)
    (here/'short_pairs.json').write_text(json.dumps({m:masks[m] for m in representatives})+'\n')
    print('normalizer modulo translations',len(norm),'support masks',len(masks),'orbits',len(representatives),
          'embedding problems',sum(len(colors[profile(m)]) for m in representatives))
if __name__=='__main__':main()
