"""Search for positions mod 11 with the stored inner labels on both signed tiles."""
import argparse
from collections import defaultdict
from ortools.sat.python import cp_model
from transfer616 import HERE,components,difference


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds',type=float,default=50)
    args=parser.parse_args()
    c=list(map(int,(HERE/'colors616.txt').read_text().split()))
    separate=list(map(int,(HERE/'positions616.txt').read_text().split()))
    tiles=components()
    model=cp_model.CpModel()
    g=[model.new_int_var(0,10,f'g{i}') for i in range(49)]
    for nodes in (range(21),range(21,49)):
        for y in range(7):
            model.add_all_different([g[i] for i in nodes if c[i]==y])
    by_colours=defaultdict(list)
    for edges in tiles:
        for i,j in edges:
            if c[i]>c[j]:
                i,j=j,i
            d=model.new_int_var(1,10,'directed_difference')
            model.add_modulo_equality(d,11+g[j]-g[i],11)
            if c[i]==c[j]:
                absolute=model.new_int_var(1,5,'undirected_difference')
                model.add_allowed_assignments([d,absolute],[(z,min(z,11-z)) for z in range(1,11)])
                d=absolute
            by_colours[c[i],c[j]].append(d)
    for (u,v),bucket in by_colours.items():
        model.add_all_different(bucket)
        if u!=v:
            model.add(sum(bucket)==55)
    i,j=tiles[0][0]
    # Independent cyclic shifts of the two tiles and one simultaneous dilation.
    model.add(g[i]==0)
    model.add(g[j]==1)
    model.add(g[21]==0)
    for i,pos in enumerate(separate):
        model.add_hint(g[i],pos)
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=args.seconds
    solver.parameters.num_search_workers=4
    status=solver.solve(model)
    print('joint model on 49 vertices, 245 differences',solver.status_name(status),flush=True)
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        z=[solver.value(x) for x in g]
        used=[difference(c,z,i,j) for edges in tiles for i,j in edges]
        assert len(set(used))==245 and all(d for _,_,d in used)
        print('verified full relative difference cover',*z)


if __name__=='__main__':
    main()
