"""Check the two proposed signed components, their colours, and their relative differences."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def components():
    d=np.zeros((7,21),dtype=int)
    for i,(p,q) in enumerate(combinations(range(7),2)):
        d[p,i],d[q,i]=-1,1
    r=np.loadtxt(HERE/'roots616.txt',dtype=int)
    assert r.shape==(28,8) and np.array_equal(r.T@r,28*np.eye(8))
    assert np.all(np.diag(r@r.T)==8) and np.all((r@r.T)%4==0)
    blocks=[2*np.eye(21,dtype=int)-d.T@d,2*np.eye(28,dtype=int)-(r@r.T)//4]
    edges=[]
    for start,x in zip((0,21),blocks):
        assert set(np.unique(x))=={-1,0,1} and np.all(x.diagonal()==0)
        assert np.array_equal(x@x,10*np.eye(len(x))-3*x)
        assert np.all(np.sum(x!=0,axis=1)==10)
        edges.append([(i+start,j+start) for i,j in combinations(range(len(x)),2) if x[i,j]])
    return edges


def difference(c,g,i,j):
    x,y=c[i],c[j]
    z=(g[j]-g[i])%11
    if x>y:
        x,y,z=y,x,(-z)%11
    if x==y:
        z=min(z,(-z)%11)
    return x,y,z


def main():
    tiles=components()
    colors=list(map(int,(HERE/'colors616.txt').read_text().split()))
    assert Counter(colors)=={i:7 for i in range(7)}
    counts=Counter(tuple(sorted((colors[i],colors[j]))) for edges in tiles for i,j in edges)
    assert all(counts[i,j]==(5 if i==j else 10) for i in range(7) for j in range(i,7))
    print('signed tiles',list(map(len,tiles)),'edges; 7 diagonal colour pairs occur 5 times, 21 off-diagonal pairs 10 times')
    labels=list(combinations(range(7),2))
    triangle=[labels.index(e) for e in ((0,1),(0,3),(1,3))]
    assert [colors[i] for i in triangle]==[5,0,6]
    d=np.zeros((7,21),dtype=int)
    for i,(u,v) in enumerate(labels):
        d[u,i],d[v,i]=-1,1
    cycle=np.zeros(21,dtype=int)
    cycle[triangle]=[1,-1,1]
    assert np.all(d@cycle==0)
    print('switching triangle',triangle,'has distinct inner labels 5, 0, 6')
    if (HERE/'positions616.txt').exists():
        pos=list(map(int,(HERE/'positions616.txt').read_text().split()))
        assert len(pos)==49 and all(0<=k<11 for k in pos)
        covers=[]
        for nodes,edges in zip((range(21),range(21,49)),tiles):
            assert len({(pos[i],colors[i]) for i in nodes})==len(nodes)
            proj={difference(colors,pos,i,j) for i,j in edges}
            assert all(d!=0 for x,y,d in proj) and len(proj)==len(edges)
            covers.append(proj)
        print('valid individual translation orbits with',11*len(covers[0]),'and',11*len(covers[1]),'base edges')
        overlaps=[len({(x,y,(min((u*d)%11,(-u*d)%11) if x==y else (u*d)%11))
                           for x,y,d in covers[0]} & covers[1]) for u in range(1,11)]
        print('smallest overlap in difference classes after rescaling one starter',min(overlaps),'of 245')


if __name__=='__main__':
    main()
