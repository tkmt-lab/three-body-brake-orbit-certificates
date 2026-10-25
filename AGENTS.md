# Computational reading guide

Read README.md, FORMAT.md, claims.json, VERIFICATION-RECORD.md and provenance.json.
Use `python3 -B check_all.py` from the repository root for the same bounded
checks as CI. Python 3.11+ and Poppler are required; no Python package install
is needed. Keep terminal algebra, saved status records, source audits, full
trace replay and a new ODE integration separate.

This is version 1.0.0. Preserve the 13 Arb rows and ten same-zero comparisons.
Do not add second-route holds 11/16 or Arb-only row15 to the paired set.
Do not call numerical residual overlap a uniqueness proof. The inequalities
must be recomputed, not accepted from a status flag. Refuse optimized Python.
Do not alter hashes to hide an unexplained data mismatch.

Optional integrators run only with `optional/run.py --run-heavy` and an
explicit bounded output directory; the default checks never import them.
Audit programs in audit/ are historical mathematical evidence, not default
verification entrypoints. Raw large traces are external to this Git packet;
do not invent assets or URLs. No network, cloud, publication or credentials
are needed for checking. Do not infer least period, stability, geometric
row distinctness, catalogue completeness, novelty priority, human peer review
or a proof-assistant theorem from a passing check.
