# Interval certification of selected equal-mass three-body periodic brake orbits

[日本語](README.ja.md)

Thirteen parameter boxes around published numerical candidates admit a
nonconstant, collision-free periodic solution of the planar equal-mass
three-body equations. In ten boxes, two flow implementations certify the
same normalized solution. Readers can inspect the [article](paper/article.pdf),
check the finite certificate algebra with exact rational arithmetic, and
regenerate all seven paper tables from the included records.

Version **1.0.0** contains 23 saved certificates: Arb for rows
**1–6, 10–16**, handwritten MPFR for **1, 2, 3, 12**, and CAPD for
**4, 5, 6, 10, 13, 14**. Equal rational centers, nested boxes, and uniqueness
in the outer box identify the same zero in each of the ten pairs. Residual
agreement alone is not the argument. Rows **11, 15, 16** have Arb-only
certification; the second route for **11 and 16** remains on hold.

## Check the records

Python 3.11+ and Poppler (`pdftotext`, `pdfinfo`) are needed for the full
packet check on macOS or Linux. The mathematical checkers use only the
Python standard library. No network access or numerical integration occurs.
Run exactly the command used by CI and the clean-clone preflight:

```sh
python3 -B check_all.py
```

The individual commands are:

```sh
python3 -B verify.py
python3 -B regenerate.py --check
python3 -B negative_controls.py
python3 -B check_packet.py
```

Success is exit code 0 with `status: PASS`. `verify.py` checks 23 finite
terminal implications and ten same-zero comparisons; `negative_controls.py`
rejects 27 deliberately corrupted inputs. `regenerate.py --check` recomputes
numerical cells and checks them against the TeX and generated fragments.
Table 1's qualitative descriptions are manual text. Incomplete-run statuses
and the three unmatched catalogue entries are saved-record checks, not new
integrations or an exhaustive catalogue search. Optimized Python is refused.

The input file hashes and export transformations are in [provenance.json](provenance.json).
[FORMAT.md](FORMAT.md) explains interval encodings and the default checks;
[REPRODUCING.md](REPRODUCING.md) separates them from optional heavy work.

## Scope of the argument

The finite checker proves implications **conditional on the saved flow and
derivative enclosures**. It recomputes determinant, contraction, strict
Krawczyk inclusion, nondegeneracy and positive saved distance bounds. It does
not regenerate ODE/Taylor coefficients or verify the compiler, arithmetic
libraries or recorded executions. CAPD's six selected paths additionally have
saved-record and fixed-source audits; their summaries, arithmetic code and
large-evidence inventory are in [audit/](audit/README.md). Giant traces are
not required by the default checks and are not bundled.

The model and quarter-period symmetry argument are shared premises of the
implementations. Local uniqueness is in a normalized parameter box; thirteen
boxes are not a count of geometrically distinct orbits. Least period,
stability, completeness and priority are not established. Seven of the ten
paired rows numerically match earlier self-dual entries; earlier labels for
rows 4, 10, 12 remain unconfirmed in the reviewed tables. These candidates
were already published numerically.

Lean level **L0**: large interval integrations and Krawczyk tests have no
practical Lean formalization route in this work, so separate implementations
are used. L0 is not formal proof. There has been no human expert peer review.
The [verification record](VERIFICATION-RECORD.md) distinguishes the previous
545 computational PASS checks from 124 saved-record comparisons, and states
the prior reader's audit experience and permitted exposure. Building this
packet is integration of shared evidence, not a new independent review.

## Research and attribution

Tatsuhiko Takimoto is responsible for the work. AI systems assisted with
implementation, calculations, mathematical discussion and manuscript editing;
the article records that use. [claims.json](claims.json) states the exact
scope, and [CITATION.cff](CITATION.cff) provides citation metadata.

Hristov and collaborators supplied the numerical brake-orbit catalogue and
its coordinates. Earlier work by Montgomery, Kapela and Zgliczyński, and
Kapela and Simó supplies geometric and computer-assisted foundations.
Johansson's Arb, python-flint, MPFR/GMP, and CAPD make these calculations
possible. I thank the authors and developers for those foundations and hope
this work adds something to the field they have developed. This acknowledgment
does not imply their review or endorsement. Full references are in the article;
[SOURCES.md](SOURCES.md) identifies the retained numerical excerpts.

Code is Apache-2.0; the author's text and generated certificate data are
CC BY 4.0. Third-party numerical excerpts and library dependencies retain
their own attribution and terms; see [LICENSE](LICENSE).
