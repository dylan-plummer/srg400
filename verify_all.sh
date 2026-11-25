#!/bin/sh
# Rebuild and verify the srg(400,168,68,72) graphs, the lifting lemma controls, and the family certificates.
# Needs Python 3.10+ with numpy, networkx, igraph, pynauty, sympy and ortools.
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
echo "== seed: rebuild both graphs from the starters, all-pairs audit, invariants"
(cd seed && $PY -B construct.py && $PY -B audit.py && $PY -B invariants.py)
echo "== audit: independent rebuild from the written description (exact comparison with the exports)"
(cd audit && $PY -B independent_rebuild.py && $PY -B cliques_check.py)
echo "== lift: general lemma on random decorations and signs, with controls (Fano F=4, K12 F=12 Paley)"
(cd lift && $PY -B general_lift.py)
echo "== family: fibre-recovery certificate, cover group, clique bound; 32 sign variants; catalogue samples"
(cd family && $PY -B structure.py && $PY -B family.py && $PY -B verify_examples.py && $PY -B check_catalogue.py)
echo "all checks passed"
# Longer runs (not executed here): family/symmetric_census.py (14,656 classes), family/full_count.py (exact
# Burnside count), family/run_full_audit.py (all 14,656 graphs on all pairs), audit/count_ti.py (888 sign-free
# translation-invariant classes), lift/scan_targets.py (target rows), lift/signing.py (F = 2 signing searches).
