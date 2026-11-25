"""Identify the full group of each cover by its action on the translations."""
from companions import here, signed_data, cover_action, json

for i,rows in enumerate(json.loads((here/'representatives.json').read_text())):
    tiles,_=signed_data(rows)
    group=cover_action(rows,tiles)
    stabilizer=[g for g in group if g[0]==0]
    for g in stabilizer:
        if all(g[p]==p for p in range(25)):
            continue
        def period(x):
            z=g[x]
            count=1
            while z!=x:
                count+=1
                z=g[z]
            return count
        if max(map(period,range(25)))!=len(stabilizer):
            continue
        u,v=divmod(g[5],5)
        w,z=divmod(g[1],5)
        linear=lambda a,b:(5*((a*u+b*w)%5)+(a*v+b*z)%5)
        assert all(g[k]==linear(*divmod(k,5)) for k in range(25))
        eigen=[l for l in range(5) if ((l-u)*(l-z)-v*w)%5==0]
        assert ((len(group)==50 and eigen==[1,4]) or
                (len(group)==100 and eigen in ([2,4],[3,4])))
        print(i,len(group),'cyclic point stabilizer',len(stabilizer),'eigenvalues mod 5',eigen)
        break
    else:
        raise AssertionError(i)
