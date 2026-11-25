"""All one-size simplex-boundary tiles in Brouwer's local snapshot of open parameters."""
import csv
from math import comb, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    with (ROOT/'current_table.csv').open() as f:
        rows = list(csv.DictReader(f))
    def params(row):
        return tuple(int(row[key]) for key in ('v','k','lam','mu'))
    table = {params(row):row['status'] for row in rows}
    for row in rows:
        if row['status'] != '?':
            continue
        n,k,l,mu = params(row)
        delta = (l-mu)**2+4*(k-mu)
        gap = isqrt(delta)
        if gap*gap != delta or (l-mu+gap)%2:
            continue
        r,s = (l-mu+gap)//2,(l-mu-gap)//2
        for h in [2,*range(4,n+1,4)]:
            if n%h:
                continue
            t,a = n//h,h//2
            for c in (0,1):
                inside = c*(h-1)
                if any(x%a for x in (r+c,s+c,k-inside)):
                    continue
                u,v = (r+c)//a,-(s+c)//a
                if min(u,v)<1:
                    continue
                b = u*v*(h-1)
                if b>=t or k!=inside+a*b:
                    continue
                size = comb(u+v,u)
                if size>t or t*(h-1)%size:
                    continue
                p,q = (r-inside)//a,(s-inside)//a
                mb = b+p*q
                lb = mb+p+q
                if b==t-1:
                    host=f'complete K{t}'
                elif (not 0<=lb<b or not 0<=mb<=b or t-2-2*b+mb<0 or t-2*b+lb<0
                      or b*(b-1-lb)!=(t-1-b)*mb):
                    host='infeasible base pair counts'
                else:
                    host=f'base {(t,b,lb,mb)}, table status {table.get((t,b,lb,mb),"absent")}'
                print((n,k,l,mu),'h',h,'internal clique',bool(c),'J',u+v,min(u,v),
                      'tiles',t*(h-1)//size,host)


if __name__=='__main__':
    main()
