# Strongly regular graphs with parameters (400, 168, 68, 72)

Graphs, code and a short note (`note/srg400.pdf`) for the results below (Dylan Plummer, 24 November 2026).

## Results

| parameters | status in Brouwer's table (checked 24 Nov 2026) | here |
|---|---|---|
| srg(400, 168, 68, 72) and complement (400, 231, 134, 132) | open (`?`, "pg(14,11,6)?") | about 6.7149 · 10^121 pairwise non-isomorphic graphs from one tiling (exact count in the note) |

- **Construction.**
  - The complement B of the Hoffman–Singleton graph splits into 35 edge-disjoint copies of T(5): three base blocks
    under C5 × C5, with orbit lengths 5, 5 and 25.
  - Each copy is signed by X = 2I − DᵀD, where D is an oriented incidence matrix of K5. X has eigenvalues 2^6 and
    (−3)^4, the restricted eigenvalues of B.
  - Every point is replaced by 𝔽₂³, and the 7 tiles at a point are matched with the 7 nonzero characters.
  - Vertices (p,x) and (q,y) are adjacent iff pq lies in a tile T and χ_{p,T}(x) χ_{q,T}(y) = X_T(p,q).
  - Spectrum: 168^1, 8^231, (−12)^168.
- **Credit.** Anton Kovaľ's search ledger records the co-Hoffman–Singleton cell graph, the eight-point slab
  construction and signed T(5) components for this row. It left the free C5 × C5 case unresolved, and the three
  starters here realize that case.
- **Lifting lemma (note, Lemma 1).**
  - Setting: any Hadamard fiber of order F, any signed tiles with two eigenvalues {ρ, σ}, and fibers that are either
    empty (c = 0) or complete (c = 1).
  - Condition: B regular with restricted eigenvalues ρ − 2c and σ − 2c, and every point in F − 1 tiles.
  - Then the lift is srg(Fv, c(F−1) + F k_B/2).
  - Controls: Fano plane, F = 4 → (28,12,6,4). K12 by edges with a Paley Hadamard matrix of order 12 → (144,66,30,30).
    This construction, F = 8 → (400,168,68,72).
- **The family.**
  - Characters and flag signs are free.
  - The fibers are intrinsic: zero-sum 4-sets of rows over 𝔽₂ are exactly the affine planes inside fibers.
  - The cover group is C5 × D10, so no member is vertex-transitive.
  - The tiling alone gives exactly N = (2^215·15^50 + 5·2^121·15^20·61² + 8·2^39·15^10 + 16·2^43·15^10 +
    20·2^25·15^4)/50 ≈ 6.7149 · 10^121 classes, by Burnside's lemma over the cover group (`family/COUNT.txt`).
  - Full automorphism groups have orders 1, 2, 5, 10, 25 and 50. Only about 10^-59 of the graphs have a nontrivial
    automorphism (`family/stabilizers.out`).
- **Translation-invariant members.** There are exactly 14,656 classes: 14,144 with |Aut| = 25 and 512 with
  |Aut| = 50. All 14,656 graphs were built, canonically labeled with nauty and checked on all pairs.
- **Cliques.**
  - Every member has clique number ≤ 10.
  - Three K5-edges through a common vertex give a negative triangle of X, and no lifted clique can project onto one.
  - Only 95 base 10-sets survive, and no 11-set.
  - So no member is the point graph of a pg(14,11,6).
- **The two exported graphs.**
  - Translation-invariant: |Aut| = 25, 4,958,400 four-cliques.
  - Swap at point 0: |Aut| = 1, 4,958,208 four-cliques.
  - Both have clique number 10.

## Other parameters

The lemma reaches several open rows up to 1300 vertices arithmetically. Each needs a tiling of an existing base by
signed graphs with a prescribed spectrum: 340/108, 400/126, 400/147, 448/162, 640/210, 640/213, 676/216, 676/243,
700/216, 868/300, 880/270, 1000/336 and 196/75. Others are blocked:
- 288 and 1200: affine planes of orders 6 and 10.
- 736/364: Fisher's inequality.
- 960/364: clique size.
- 448/162 with 28-point tiles: GQ(3,9) has no ovoids.
- Several rows whose base does not exist.

The F = 2 rows ask for a single signing of a base graph:
- (210,33,0,6) from the flag graph of PG(2,4), the unique srg(105,32,4,12). CP-SAT says INFEASIBLE in 656 s, so this
  route is closed.
- (162,21,0,3) from the Brouwer–Haemers graph srg(81,20,1,6). The general run was undecided after 30 minutes, and no
  translation-invariant signing exists.

`lift/signing.py` runs these searches, and its positive control signs K4,4 and lifts it to Clebsch (16,5,0,2). The
outputs are in `lift/signing.out`.

A mixed T(7)/E8 tiling of K_{11×7} would give (616,287,126,140) and,
by switching, (616,300,152,140). Its separate translation orbits are valid, but no joint decomposition has been
found yet. See `lift/scan_targets.out` and `family/TRANSFERS.txt`.

## Files

- `data/starters.txt`: the three base blocks. `data/brouwer_table_snapshot.csv`: the table snapshot used by the
  scans.
- `graphs/`: graph6 files, vertex 8p + x.
  - `srg_400_168_68_72.g6`: the translation-invariant graph.
  - `srg_400_168_68_72_swap.g6`: the asymmetric graph.
  - `srg_400_231_134_132*.g6`: their complements.
  - `SHA256SUMS`.
- `seed/`: the original constructor (standard library only), all-pairs audit and invariants (nauty, bliss,
  four-cliques, clique numbers), with outputs.
- `audit/`: independent checks.
  - `independent_rebuild.py`: rebuild from the written description, with exact comparison against the exports and
    random decorations with signs.
  - `cliques_check.py`: second enumeration of admissible base cliques.
  - `count_ti.py`: the 888 sign-free translation-invariant classes by nauty.
  - `cross_catalogue.py`: comparison with the catalogue.
  - `fibre_bt.py`: an exhaustive search for 8-sets that behave like fibers, run on sample members.
- `lift/`:
  - `general_lift.py`: the lemma for arbitrary base, tiles and Hadamard matrix, with controls.
  - `scan_targets.py`: target rows for c ∈ {0,1} and Hadamard orders.
  - `signing.py`, `pg24flags.py`, `gf.py`: the F = 2 signing searches.
- `family/`: structure and counting (snapshot of the companion analysis, 24 Nov 2026 18:10).
  - `structure.py`: fiber-recovery certificate, cover group, clique bound.
  - `symmetric_census.py`, `catalogue.tsv`: the 14,656 translation-invariant classes. `export.py N` writes class N
    as graph6.
  - `full_count.py`, `COUNT.txt`: exact Burnside count. `stabilizers.py`: count by subgroups, automorphism-group distribution.
  - `run_full_audit.py`, `audit_catalogue.cpp`, `canonical.tsv.gz`: all 14,656 graphs on all pairs, with canonical digests.
  - `catalogue_cliques.py`: clique numbers.
  - `TRANSFERS.txt`, `transfers.py`, `*616*`: the transfers to other rows.
  - Sample graphs: `translation_signs.g6` (32 sign variants), `group_samples.g6`, `srg400_order50.g6`,
    `srg400_clique8.g6`, `srg400_clique9.g6`.
  - The detailed write-up is `family/README.txt`.

## Verification

    ./verify_all.sh            # about 1.5 minutes; PYTHON=... to choose the interpreter

This script does the following:
- rebuilds both graphs from the starters twice (`seed/construct.py`, and `audit/independent_rebuild.py` written only
  from the description), compares them with the stored files byte for byte, and checks
  A² = 96I − 4A + 72J on all pairs;
- runs the lifting lemma on random decorations and the three controls;
- reruns the fiber-recovery certificate, the cover group and the clique enumeration;
- checks the 32 sign variants and 18 catalogue representatives.

The full census, exact count and catalogue audit take longer and are listed at the end of `verify_all.sh`.
Requirements: Python ≥ 3.10 with numpy, networkx, python-igraph, pynauty, sympy and ortools.
