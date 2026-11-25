"""Construct the 400-vertex graph from three ordered ten-tuples. Standard library only."""
from itertools import combinations, product
from pathlib import Path
from collections import Counter

HERE = Path(__file__).resolve().parent
GRAPHS = HERE.parent / 'graphs'
EDGES = list(combinations(range(5), 2))
TRANSLATIONS = list(product(range(5), repeat=2))


def shift(label, shift_by):
    sheet, position = divmod(label, 25)
    row, col = divmod(position, 5)
    a, b = shift_by
    return 25 * sheet + 5 * ((row + a) % 5) + (col + b) % 5


def directed_incidence(edge):
    i, j = edge
    if (j - i) % 5 in (1, 2):
        return {i: -1, j: 1}
    return {j: -1, i: 1}


def signed_edges(labels):
    assert len(set(labels)) == 10
    roots = [directed_incidence(e) for e in EDGES]
    result = []
    for i, j in combinations(range(10), 2):
        common = roots[i].keys() & roots[j].keys()
        if common:
            x, = common
            a, b = sorted((labels[i], labels[j]))
            result.append((a, b, -roots[i][x] * roots[j][x]))
    assert len(result) == 30
    return tuple(sorted(result))


def translate(tile, translation):
    shifted = []
    for a, b, sign in tile:
        u, v = sorted((shift(a, translation), shift(b, translation)))
        shifted.append((u, v, sign))
    return tuple(sorted(shifted))


def build(starters, swap_at_zero=False):
    tiles = set()
    lengths = []
    for row in starters:
        orbit = {translate(signed_edges(row), t) for t in TRANSLATIONS}
        lengths.append(len(orbit))
        assert tiles.isdisjoint(orbit)
        tiles.update(orbit)
    assert lengths == [5, 5, 25]
    # Sort by unsigned support. Transport seven characters from the origin
    # of each sheet, so every translation also preserves the resulting graph.
    def support(tile):
        return tuple((a, b) for a, b, sign in tile)
    assert len({support(t) for t in tiles}) == 35
    flags = {}
    for origin in (0, 25):
        incident = sorted((tile for tile in tiles if any(origin in (a, b) for a, b, _ in tile)), key=support)
        assert len(incident) == 7
        for t in TRANSLATIONS:
            label = shift(origin, t)
            for char, tile in enumerate(incident, 1):
                flags[label, translate(tile, t)] = char
    assert len(flags) == 350
    if swap_at_zero:
        for tile in tiles:
            char = flags.get((0, tile))
            if char in (1, 2):
                flags[0, tile] = 3-char
    rows = [0] * 400
    pairs = set()
    for tile in tiles:
        for a, b, sign in tile:
            assert (a, b) not in pairs
            pairs.add((a, b))
            for x, y in product(range(8), repeat=2):
                parity = ((x & flags[a, tile]).bit_count() + (y & flags[b, tile]).bit_count()) % 2
                if (parity == 0) == (sign == 1):
                    i, j = 8*a + x, 8*b + y
                    rows[i] |= 1 << j
                    rows[j] |= 1 << i
    assert len(pairs) == 1050
    return rows, tiles


def graph6(rows):
    n = len(rows)
    assert 63 <= n <= 258047
    out = bytearray([126] + [63 + ((n >> s) & 63) for s in (12, 6, 0)])
    bits = 0
    count = 0
    for col in range(n):
        for row in range(col):
            bits = 2 * bits + ((rows[row] >> col) & 1)
            count += 1
            if count == 6:
                out.append(63 + bits)
                bits = count = 0
    if count:
        out.append(63 + (bits << (6-count)))
    return bytes(out) + b'\n'


def main():
    starters = [list(map(int, s.split())) for s in (HERE/'starters.txt').read_text().splitlines()]
    complete = (1 << 400) - 1
    for suffix in ('', '_swap'):
        rows, _ = build(starters, swap_at_zero=bool(suffix))
        assert all(r.bit_count() == 168 and not (r >> i & 1) for i, r in enumerate(rows))
        histogram = Counter(((rows[i] >> j) & 1, (rows[i] & rows[j]).bit_count())
                            for i, j in combinations(range(400), 2))
        assert histogram == {(1, 68): 33600, (0, 72): 46200}
        complement = [complete ^ r ^ (1 << i) for i, r in enumerate(rows)]
        (GRAPHS/f'srg_400_168_68_72{suffix}.g6').write_bytes(graph6(rows))
        (GRAPHS/f'srg_400_231_134_132{suffix}.g6').write_bytes(graph6(complement))
        print(suffix or 'translation-invariant', '35 tiles, 1050 cell pairs, 350 flags', sorted(histogram.items()))


if __name__ == '__main__':
    main()
