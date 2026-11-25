"""Check that the submitted roots generate E8 and have the required frame operator."""
from pathlib import Path
from sympy import Matrix,eye
from sympy.matrices.normalforms import hermite_normal_form


def main():
    p=Path(__file__).with_name('roots616.txt')
    r=Matrix([list(map(int,line.split())) for line in p.read_text().splitlines()])
    assert r.shape==(28,8) and r.T*r==28*eye(8)
    for i in range(28):
        v=list(r.row(i))
        assert len({int(x)%2 for x in v})==1 and sum(v)%4==0 and sum(x*x for x in v)==8
    determinant=abs(hermite_normal_form(r.T).det())
    assert determinant==2**8
    print('The doubled roots span a lattice of covolume',determinant)
    print('After division by 2 the lattice is E8 and the frame operator is 7I.')


if __name__=='__main__':
    main()
