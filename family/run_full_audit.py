"""Build the independent nauty checker, check every graph, and compare the canonical digests."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import argparse
import gzip
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nauty',type=Path,default=Path(sys.prefix)/'share/gap/pkg/grape/nauty2_8_6')
    parser.add_argument('--workers',type=int,default=6)
    args=parser.parse_args()
    env=dict(os.environ,CCACHE_DISABLE='1')
    with TemporaryDirectory() as work:
        build=Path(work)
        def cc(name):
            out=build/(name+'.o')
            subprocess.run([env.get('CC','gcc'),'-O3','-march=native','-DWORDSIZE=64','-c',
                            str(args.nauty/(name+'.c')),'-o',str(out)],env=env,check=True)
            return str(out)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            objects=list(executor.map(cc,['nauty','nautil','naugraph','nausparse','schreier','naurng']))
            binary=build/'audit'
            subprocess.run([env.get('CXX','g++'),'-O3','-march=native','-std=c++17','-DWORDSIZE=64',
                            '-I'+str(args.nauty),str(HERE/'audit_catalogue.cpp'),*objects,'-lcrypto','-o',str(binary)],
                           env=env,check=True)
            seed=HERE.parent/'graphs'/'srg_400_168_68_72.g6'
            def check(first):
                cmd=[str(binary),str(seed),str(HERE/'catalogue.tsv'),str(first),str(min(first+256,14656))]
                return subprocess.check_output(cmd,env=env)
            data=b''.join(executor.map(check,range(0,14656,256)))
    rows=[line.split() for line in data.splitlines()]
    assert [int(r[0]) for r in rows]==list(range(14656))
    assert len({r[2] for r in rows})==14656
    counts=Counter(int(r[1]) for r in rows)
    assert counts=={25:14144,50:512}
    table = b'class\tgroup_order\tcanonical_sha256\n' + b''.join(b'\t'.join(r)+b'\n' for r in rows)
    (HERE/'canonical.tsv.gz').write_bytes(gzip.compress(table,mtime=0))
    print('Distinct ordinary 400-vertex canonical adjacency digests',len(rows))
    print('Automorphism group order counts',dict(sorted(counts.items())))
    print('Unordered off-diagonal pairs checked',14656*79800)


if __name__=='__main__':
    main()
