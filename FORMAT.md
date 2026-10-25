# Certificate and table formats

`data/index.json` lists 23 selected certificates by family (A = Arb, B =
handwritten MPFR, C = CAPD), original one-based catalogue row, exact decimal
center and power-of-two radius exponent. Each path names a JSON file under
`data/certificates/`. Original file hashes and any removed local metadata
fields are recorded in `provenance.json`. JSON reserialization is not claimed
to be byte-identical to the original; retained numerical payloads are unchanged.

A dyadic is `n * 2**e`. It is represented as `[n,e]` or `{"n":n,"e":e}`;
integer strings are converted to integers before use. Arb/MPFR balls use
`{"mid":dyadic,"rad":dyadic}` with a nonnegative radius. CAPD intervals are
two dyadic endpoints `[lower,upper]`. Arbitrary rational values use a pair
of integer strings `[numerator,denominator]`, in fields identified below.
The decoder does not infer encodings from decimal approximations.

For A and B, `krawczyk.C`, `.J` and `.F` contain the finite matrices and
residual enclosure; A's residual entries have an additional singleton-list
axis. `nondegeneracy_interval` is a ball. `distance_squared_lower` is a ball
whose lower endpoint is used. `steps` and `integrated_time` are rational
pairs; their sum is checked exactly. `box_radius_power_two` sets the radius.

For C, `terminal.mpfr.C`, `.J`, `.F0` and `.nondegeneracy` contain the terminal
intervals. `all_time_distance2_lower` is a rational pair;
`parameter_radius_exponent` sets the box. The normalized clock reaches 128
and the variable physical period is included. Center row 13 is additionally
bound to `data/catalogue/C13-row.txt`. Positive slab counts are retained
coverage data, not a new replay of each slab.

The checker recomputes `det(C) != 0`, `Z = I-CJ`, the infinity row-sum bound
`q < 1`, and `K_delta = -CF + Z[-r,r]^3` strictly inside `(-r,r)^3` using
Python Fraction interval endpoints. It checks that the nondegeneracy interval
excludes zero, that the saved distance-squared bound is positive, that the
period is positive throughout the box, and that the initial velocity cannot
vanish. The same-zero comparison checks exact centers and box containment;
outer-box uniqueness supplies the identification. The enclosures and their
connection to the differential equation remain premises of this finite check.

`data/settings/` contains selected fields of saved run configurations.
`data/status/` contains projections of retained unsuccessful/partial/pilot
records. Empty inventories are not proofs of nonexistence. Templates under
`data/table-templates/` preserve headings, captions and qualitative text;
`regenerate.py` creates complete table fragments under `generated/` and
checks exact equality with the seven blocks in `paper/article.tex`.

Table 4's linear distance lower bound is rounded down to six decimals by an
integer square-root calculation. Table 5 displays a square-root approximation
and a contraction approximation. Table 6 displays a squared-distance
approximation and the saved `q_upper`, not the complete inclusion ratio.
These display roundings are not the inequalities used by the checker.
Table 7 recomputes seven positive matches using the retained original decimal
tokens and the bound `5e-79`; the three nonmatches and catalogue totals are
saved-record evidence. Table 1 contains manually maintained prose about
proof layers; it is not mechanically proved by table regeneration.
