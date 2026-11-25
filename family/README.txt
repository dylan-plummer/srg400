The size and geometry of the (400,168,68,72) family
24 November 2026

On the three-starter decomposition alone there are approximately
6.71493814 x 10^121 isomorphism classes of strongly regular graphs.
The exact integer, its formula and two counts are in COUNT.txt.
Exactly 14,656 classes arise from
the decorations that repeat under the specified 25 translations. Their
full automorphism groups have orders 25 (14,144 classes) and 50 (512).
All the graphs from this fixed decomposition have clique number at most
10 and at least two vertex orbits. These claims concern this cover.

The original three starters and the proof of the graph parameters are
in ../seed/. Anton Kovaľ's co-Hoffman–Singleton
slab and signed triangular graph framework is credited there. For the
different open parameters on 616 points, see TRANSFERS.txt.

1. What may vary

Write B for the complement of the Hoffman–Singleton graph, partitioned
into 35 edge-disjoint triangular graphs T(5). At each of its 50 points
put a fiber F2^3. There are seven incident tiles at each point and ten
points in each tile. A flag is an incident pair (point,tile). To its
seven flags assign the seven nonzero binary characters in any order.
Also choose, independently for each flag, a multiplier +1 or -1 for
its character vector. On each tile this merely conjugates the signed
10 by 10 adjacency matrix by a diagonal sign matrix. Each tile keeps
the restricted eigenvalues 2 and -3. Character orthogonality and the
fiber quotient are unchanged. The original proof therefore gives
A^2 = 96I - 4A + 72J for every choice.

2. Recovering the fibers from the abstract graph

Use adjacency rows over F2. On a target fiber q the row of (p,x) is
zero for an inactive base pair. For an active pair it is an affine
function of y whose linear part is the nonzero character assigned at
q to the tile carrying pq. The constant depends on x and the signs.

Consider four distinct vertices. Let S consist of the base points used
an odd number of times. If |S|=2 or 4, the linear part at q is the XOR
of the characters on the tiles at q that meet S in an odd number of
neighbors of q. structure.py checks the three starters and finds, for
each of the 1,225 pairs and 230,300 four-sets S, at least one q with
exactly one or two such tiles. Distinct characters cannot cancel in
either case. This certificate is independent of every character order
and every flag sign.

When S is empty, two vertices at p differ by a nonzero vector. Their
row difference is identically one on the neighboring fibers reached
through exactly four of the seven directions at p and zero elsewhere.
For all 50 points and all 35 choices of four directions, the 1,750 sets
of 24 neighboring base points are pairwise distinct. Thus two rows
at p and two at a different point cannot sum to zero. For four vertices
over a single point, their rows sum to zero exactly when their four
coordinates sum to zero. This follows since all seven characters occur.

Consequently the four-sets with zero row sum are exactly the 14 affine
planes in each fiber, giving 700 altogether. Two distinct vertices are
in the same fiber exactly when they occur together in one of these
four-sets. Any isomorphism therefore carries fibers to fibers and acts
affinely on each of them. Indeed, after fixing the origin a map that
preserves the four-sets {0,x,y,x+y} is linear.

The abstract graph also recovers B by its nonempty blocks. At a point,
two neighboring fibers use the same tile exactly when they induce the
same unordered partition into two halves at that point. Splits for
different characters are distinct. Link two base edges through their
common endpoint when their splits agree. The connected components of
this relation are the 35 tiles because each T(5) is connected.

3. Counting isomorphism classes

Encode the recovered cover as a three-colored graph. Its vertices are
the 50 base points, the 35 tiles and the 1,050 base edges. An edge-node
is adjacent to its two endpoints and its tile. The action on the points
is faithful. The tile and its member edges determine each other, and
an edge-node is uniquely determined by its endpoints. nauty and bliss
both return group order 50 with two point orbits of length 25.

The translations form a normal subgroup of index two. An element outside
it can be composed with a translation to fix point 0. Its square is a
translation fixing 0, hence the identity. Conjugation on the translation
plane then has one fixed and one negated dimension, as verified on the
25 points of the first sheet. Thus the cover group is C5 x D10, where
D10 is the dihedral group of order ten. No graph from this cover can be
vertex-transitive since every isomorphism preserves the recovered cover.

With the vertices labeled, the character assignments give (7!)^50
choices. The seven unordered splits in each fiber recover the characters
uniquely in the given coordinates. Once these are known, evaluating a
block at (0,0) recovers the product of its two flag signs. On a connected
tile all signs are recovered up to reversing all ten of them. The number
of distinct labeled graphs in the construction is therefore

                         (7!)^50 * 2^(350-35).

Each unlabeled class accounts for at most 50 * (8*168)^50 labeled
graphs here. An isomorphism must come from a cover automorphism and a
choice in AGL(3,2) at each fiber. Dividing proves the asserted bound.
This division only uses an upper bound for the size of each class.

4. Counting the translation-invariant subfamily exactly

Transport the assignments and signs from the origin of each sheet by
the fixed translation action, as in the original graph. The seven ports
at an origin are ordered as in the seed. The seven characters give a
binary length-seven simplex code: the eight vectors of sign changes
arising by translating the local F2^3. Two ordered assignments give the
same code exactly when their column labels differ by GL(3,2), since
their full-rank binary generator matrices have the same row space. This
leaves 7!/168 = 30 codes per origin.

There are fourteen invariant flag signs, seven per sheet. Once the two
codes have been chosen, the signs over the identity base map are taken
modulo local translations and reversing all the signs on a tile. For
translation-invariant differences it suffices to use invariant sign
changes of both kinds: sum any solution of the binary equations over
the group of odd order 25. The invariant local translations supply the
two length-seven simplex codes and invariant tile reversals supply one
vector for each of the three tile orbits. The resulting subspace U has
dimension nine in all 900 cases. Thus there are 30^2 * 2^5 = 28,800
classes over the identity on the base.

Translations already act trivially on these classes. The remaining
cover coset acts by permuting the two sets of seven ports internally
and adding a sign vector gamma accounting for the orientation of K5.
symmetric_census.py derives that vector from the signed edge equations,
and enumerates the orbits of this involution. It also counts the fixed
points by a second calculation: for a preserved pair of codes and port
permutation P, solve (I+P)z = gamma modulo U. If gamma lies in
U + image(I+P), there are 2^(14 - dim(U + image(I+P))) fixed cosets.
There are 512 fixed points in total, giving (28,800+512)/2 = 14,656.

The kernel on the fibers is measured separately. A vertex permutation
over the identity on the base fixes each port, so its linear part is
the identity in every fiber. A local translation preserves the signed
edges exactly when its flag signs are constant on each connected tile.
The script solves for 35 tile constants with the seven entries at every
point constrained to its simplex code. The kernel has dimension zero
for all 900 pairs. Hence the fixed classes have the full cover group
C5 x D10. The other classes have group C5 x C5. The catalogue classifies
decorations repeating under the prescribed translations, up to arbitrary
graph isomorphism. The lower bound in section 3 allows arbitrary choices.

5. Why no member of the family is geometric

A clique projects to a clique in B and contains at most one point from
any fiber. In a tile, three K5 edges incident with one K5 vertex form
a triangle with sign product -1. This follows by multiplying the three
negative incidence products there. A lifted triangle would have sign
product +1 because each of its three character values occurs twice.
Switching signs at flags preserves this obstruction.

There are five such four-edge stars on each of the 35 tiles. A projected
clique therefore meets each of the 175 stars in at most two points.
structure.py enumerates sets in increasing vertex order, rejecting a
new base point if it is a nonneighbor or makes a forbidden triple.
The largest surviving sets are the 95 ten-sets. The bound of ten
is attained in the previously exported graphs. A line of pg(14,11,6)
would contain 15 points. This excludes that geometry throughout the
family and leaves open its existence by other constructions.

The distribution in the translation-invariant catalogue is

    clique number       group order 25     group order 50       total
         8                       12                 17             29
         9                      180                 47            227
        10                   13,952                448         14,400

To test whether a projected set lifts to a clique, give each of its
points three unknown fiber coordinates over F2. Every pair imposes one
affine linear equation on these coordinates, read directly from the
signed block. A simultaneous solution gives a clique. The admissible
sets of sizes 8, 9 and 10 form 1,577, 114 and 7 translation orbits,
respectively. catalogue_cliques.py performs binary elimination for
each required orbit and the corresponding character assignments.
The output covers every class. Classes 384, 386 and 387 all have
the same character assignments and full group of order 50. Their
clique numbers are 10, 8 and 9 respectively. check_clique_examples.py
also verifies these three numbers by ordinary graph clique search.

6. Data and reproduction

catalogue.tsv contains one row per class in the symmetric subfamily.
P_characters and Q_characters give the character numbers in port order
at the two sheet origins. Bit j of negative_flags reverses the j-th
port in zero-based order on P followed by Q. Class numbers start at 0.
To export class 384, a graph with group C5 x D10, from this directory:

    python export.py 384 > graph384.g6

This exporter requires only standard Python 3.10+ and the sibling seed
directory. The full construction and checks can also be run as follows.

    python family.py
    python structure.py
    python symmetric_census.py
    python full_count.py
    python stabilizers.py
    python verify_examples.py
    python check_catalogue.py
    python catalogue_cliques.py
    python check_clique_examples.py
    python run_full_audit.py

translation_signs.g6 contains 32 graphs, one per line. They represent
all sign cosets for the original pair of local codes before the cover
involution. Full 400-vertex nauty canonical forms distinguish all 32
from one another and from the earlier asymmetric graph. The first is
the earlier symmetric graph, so the file gives 31 additional examples.
group_samples.g6 contains class 0 and class 384. Class 384 is also new
relative to all 32 because its group has order 50. The exact integer
matrix check was run on every graph in translation_signs.g6.

Separately, audit_catalogue.cpp decoded the original graph6 file to
recover the block data, without using the starter code. It rebuilt all
14,656 catalogue entries as ordinary graphs, checked every degree and
all 1,169,548,800 unordered distinct vertex pairs across them, computed
their full nauty groups, and obtained 14,656 different canonical matrix
digests. full_audit.out and canonical.tsv.gz record the results. The
digest input is the canonical row-major dense nauty matrix with 64-bit
setwords in native byte order as used on the x86-64 verification host.
run_full_audit.py requires C and C++ compilers, OpenSSL headers and a
configured nauty source directory. The default is the GAP/grape copy
in the srg conda environment. check_catalogue.py also checks a sample
with bliss. The family reconstruction itself uses standard Python.

Related literature on recovering structure and counting graphs in
another prolific strongly regular graph construction includes Cai, Guo,
Gavrilyuk and Ponomarenko, A large family of strongly regular graphs
with small Weisfeiler-Leman dimension. https://arxiv.org/abs/2312.00460
