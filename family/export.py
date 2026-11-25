"""Export one numbered class from catalogue.tsv to graph6 on standard output."""
import argparse
import csv
import sys
from family import HERE,realise,tile_data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('number',type=int)
    args=parser.parse_args()
    with (HERE/'catalogue.tsv').open() as f:
        row=next((row for row in csv.DictReader(f,delimiter='\t') if int(row['class'])==args.number),None)
    if row is None:
        parser.error('class number must be between 0 and 14655')
    tiles,flags=tile_data()
    labels=tuple(tuple(map(int,row[k])) for k in ('P_characters','Q_characters'))
    sys.stdout.buffer.write(realise(tiles,flags,labels,int(row['negative_flags'])))


if __name__=='__main__':
    main()
