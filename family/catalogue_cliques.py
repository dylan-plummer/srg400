"""Determine clique numbers by solving the binary lift equations on projected cliques."""
from collections import Counter, defaultdict
from itertools import combinations
import csv
from family import HERE,tile_data,construct


def supports(tiles):
    pairs={}
    forbidden=defaultdict(set)
    for t,tile in enumerate(tiles):
        support={(p,q) for p,q,_ in tile}
        nodes=sorted({p for e in support for p in e})
        for p,q,s in tile:
            pairs[p,q]=t,int(s<0)
        for four in combinations(nodes,4):
            if all(e in support for e in combinations(four,2)):
                for p,q in combinations(four,2):
                    forbidden[p,q].update(set(four)-{p,q})
    reps=defaultdict(list)
    for_size=set()
    def visit(chosen,future):
        if len(chosen)>=8 and chosen not in for_size:
            reps[len(chosen)].append(chosen)
            for d in construct.TRANSLATIONS:
                for_size.add(tuple(sorted(construct.shift(p,d) for p in chosen)))
        for i,p in enumerate(future):
            rest=tuple(q for q in future[i+1:] if (p,q) in pairs and all(q not in forbidden[v,p] for v in chosen))
            visit(chosen+(p,),rest)
    visit((),tuple(range(50)))
    return reps,pairs


def equations(s,pairs,flags):
    rows=[]
    for i,j in combinations(range(len(s)),2):
        p,q=s[i],s[j]
        t,negative=pairs[p,q]
        u,v=7*(p//25)+flags[p,t]-1,7*(q//25)+flags[q,t]-1
        rhs=(negative<<14)^(1<<u)^(1<<v)
        rows.append((15+3*i,15+3*j,u,v,rhs))
    return rows


def conditions(rows,chars):
    # Eliminate the fibre coordinates first. Remaining rows constrain only signs.
    basis={}
    checks=[]
    for i,j,c,d,rhs in rows:
        w=(chars[c]<<i)^(chars[d]<<j)^rhs
        while w:
            k=w.bit_length()-1
            if k in basis:
                w^=basis[k]
                continue
            basis[k]=w
            if k<15:
                checks.append(w)
            break
    return checks


def main():
    tiles,flags=tile_data()
    reps,pairs=supports(tiles)
    lifts={n:[equations(s,pairs,flags) for s in subsets] for n,subsets in reps.items()}
    print('translation orbits of projected sets', {n:len(v) for n,v in sorted(lifts.items())},flush=True)
    families=defaultdict(list)
    with (HERE/'catalogue.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            chars=tuple(map(int,row['P_characters']+row['Q_characters']))
            families[chars].append((int(row['class']),int(row['negative_flags']) | (1<<14)))
    sizes={}
    for chars,variants in families.items():
        for n in sorted(lifts,reverse=True):
            for rows in lifts[n]:
                if not variants:
                    break
                checks=conditions(rows,chars)
                pending=[]
                for i,z in variants:
                    if all((w&z).bit_count()%2==0 for w in checks):
                        sizes[i]=n
                    else:
                        pending.append((i,z))
                variants=pending
        assert not variants
    assert len(sizes)==14656
    with (HERE/'clique_numbers.tsv').open('w') as f:
        print('class','clique_number',sep='\t',file=f)
        for i,n in sorted(sizes.items()):
            print(i,n,sep='\t',file=f)
    print('clique number distribution',dict(sorted(Counter(sizes.values()).items())),flush=True)


if __name__=='__main__':
    main()
