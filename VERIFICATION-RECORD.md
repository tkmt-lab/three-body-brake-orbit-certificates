# Verification record — version 1.0.0

The article gives a computer-assisted prose argument for 13 row-associated
boxes, with ten same-zero comparisons. `verify.py` checks the finite rational
implications in 23 saved certificates. The matrices F/J and saved distance
enclosures are inputs, not conclusions of a fresh integration in this packet.

## Layers of evidence

Arb certifies rows 1–6 and 10–16. Handwritten MPFR certifies 1, 2, 3, 12;
CAPD certifies 4, 5, 6, 10, 13, 14. Exact centers and nested radii identify
the same zero under outer-box uniqueness. The three implementations share
the model and symmetry proof; Arb and handwritten MPFR also share the Taylor
sufficient conditions. MPFR and CAPD use MPFR scalar arithmetic.

The CAPD audits recomputed saved local inclusions, time choices, reachable
family/tangent connections, derivative reductions and terminal data for the
six selected paths. They did not independently regenerate Taylor coefficients,
rebuild cloud binaries or prove the software. Fixed-source Taylor semantics,
MPFR/GMP, CAPD interval/AD/set operations, compiler, hardware, and correspondence
of fixed sources to execution records remain trusted. Coverage summaries and
archival checking code are in `audit/`. Default packet verification does not
replay those large streams. Missing full-trace assets cannot be replaced by
rechecking a PASS flag.

A general downward-rounding claim in CAPD's shrinking-step formula has a
counterexample: at 320 bits and order 48, a factor `1-2^-100` yields a
mathematically subunit root, yet the recorded clearMantissaBits conversion
returns 1. This finding is retained. It is not a counterexample to the orbit
existence statements. The six adopted paths use normal successful returns;
their audit does not justify unobserved shrinking or refined branches.

For rows 4 and 14 the adopted J comes from the saved forward reduction.
Row 10 is the p320/r140 direct path, not an earlier lifted continuation.
C0 point evaluations may be reused between parameter radii and do not prove
uniqueness over the C0 filename's radius label. Row 13's original solver OS
exit was not observed; retained terminal/reduction evidence is the basis used.
A study-wide trial-checker TypeError remains a historical exception, not a
replacement for the selected checks.

## Partial evidence and history

Rows 11, 15 and 16 retain Arb-only existence certification. The second route
on rows 11/16 remains on hold. For row11, saved audit counts were 813/5919
C0, 388/5741 C1-prefix and 0/278 tail: 1201/11938 checked, 10737 unchecked.
For row16, 5955 C0 and 5997 C1 steps were checked, with 25 C1 tail steps
unconfirmed and an unresolved initial compressed-hash/terminal connection.
Old terminal PASS records for these rows are not adopted certificates here.
Row15's terminal failure and incomplete flow remain distinct from success.

The selected Arb row1 has a documented full same-algorithm rerun after a
rejection-guard revision, agreeing on saved mathematical fields. That was
not a replay of the original frozen CLI. No full replay of all selected
Arb rows is claimed by this observation or by the default packet check.

The prior mathematical reading recorded **545 PASS** and **124
PASS_SAVED_RECORD**, not 669 newly proved numerical validations. The reader
had existing method-audit experience and allowed exposure to a summary of a
previous report. It was not a blinded new session or human review. An earlier
inconsistent all-values-match reader conclusion was declined; it is not used
as evidence. This packet adapts already-shared results and is not a new
independent mathematical audit.

## Time units and catalogue limits

The article's 241/1963/2338/507 seconds for B rows 1/2/3/12 are retained.
Rows 1/12 use certificate elapsed_seconds; rows 2/3 round worker wall times
1963.022/2337.644, distinct from native times 1962.25/2336.77. C's timing
examples concern the old order-24 rows 4/13/14. They are not a performance
benchmark or an audited cloud invoice total.

Seven positive catalogue correspondences are numerical decimal-token matches.
Earlier labels for rows 4, 10 and 12 were not identified in the reviewed tables.
The default checks do not repeat the full old-catalogue search or establish
priority. These are already-published numerical candidates, not newly
discovered orbits. Local uniqueness does not establish geometric distinctness,
least period, stability or completeness.

Lean level L0 records the absence of a practical formalization route here
for the large interval integrations and Krawczyk argument; separate
implementations are used instead. It is not a formal theorem. There is no
human expert peer review. The version date is the production date of this
article, not evidence of a future public release, DOI assignment or hosted
CI result.
