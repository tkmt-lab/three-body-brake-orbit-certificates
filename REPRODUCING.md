# Reproducing the checks

From a clean clone on macOS/Linux, with Python 3.11+ and Poppler installed:

```sh
python3 -B check_all.py
```

No pip dependency, solver, API, private cache or repository outside this clone
is used. The same command is in `.github/workflows/verify.yml`. Individual
commands and their scope are in README.md. The check includes manifest hashes,
privacy of text/PDF metadata and extracted text, exact finite algebra,
negative controls, and deterministic table regeneration. CPU/time limits are
stated by the runner and workflow; numerical work is sequential and single-thread.
Only the optional paths below require scientific libraries. CI pins the checkout
and Python setup actions to commit hashes and Python to 3.11.14. Poppler comes
from the Ubuntu 24.04 distribution package; its security-updated binary is
not pinned, and its role is PDF text inspection, not numerical algebra.

To regenerate table fragments after a deliberate source-data change:

```sh
python3 -B regenerate.py
python3 -B regenerate.py --check
```

Generation refuses disagreement with the TeX table blocks. A data change is
not validated by updating a hash manifest. Review changed numeric inputs and
the mathematics separately. `python3 -B check_packet.py --write-manifest`
updates file hashes only when explicitly requested; ordinary checking never
repairs mismatches. `.git`, bytecode, `.DS_Store` and the explicit `work/`
output directory are excluded. No other unlisted file is silently ignored.

## Optional heavy integration — not run in CI or for this package build

`optional/arb/` freezes the project-owned Taylor/Arb implementation and driver.
The dependency is python-flint 0.9.0 (observed FLINT 3.6.0); its full observed
environment is included. The driver is adapted at execution only to read the
20-row catalogue excerpt. The excerpt's whole-file hash differs from the old
full-table hash while the selected numerical tokens are identical.

For example, with the dependency already installed in the Python being used:

```sh
python3 -B optional/run.py --run-heavy --family A --row 1 --seconds 1800 --output work/arb-row1
```

The default mode refuses to run without `--run-heavy`. Choose a new output
directory: existing results are never overwritten. The process group is stopped
at the requested wall time or 2 GB sampled RSS; this is polling, not a hard
allocator limit. All output is saved. A time-limited partial run is not a
counterexample. Set thread environment variables to one as in CI.

`optional/mpfr/src/` freezes the separate MPFR endpoint implementation.
The observed libraries are MPFR 4.2.2 and GMP 6.3.0. Build on a machine with
those development libraries, supplying the appropriate include/library flags:

```sh
mkdir -p work/bin
c++ -std=c++17 -O2 -DFAST_INTERVAL optional/mpfr/src/certify.cpp -lmpfr -lgmp -o work/bin/mpfr-certify
python3 -B optional/run.py --run-heavy --family B --row 2 --seconds 3600 --binary work/bin/mpfr-certify --output work/mpfr-row2
```

The original Makefile records the observed Homebrew layout; the command above
uses compiler search paths. `optimization.json` describes baseline versus
sign-selected interval multiplication. Omit `-DFAST_INTERVAL` for a baseline
build; the saved comparison checked 30,000 multiplication pairs. This
package build has not rebuilt or rerun either binary, and source provenance
is not binary reproducibility.

`optional/capd/` freezes the project-owned CAPD drivers. Install/build the
external CAPD revision `37057e5899bccaa6cd5701057b7e6ea5e232ef8b` separately
with MPFR support; no CAPD library or modified third-party header is vendored.
The observed CAPD build used Debian g++ 12.2.0, MPFR 4.2.0 and GMP 6.2.1.
Use its `capd-config --cflags --libs` to compile the chosen driver:

```sh
c++ optional/capd/batch1_solver.cpp -o work/bin/capd-batch1 $(capd-config --cflags --libs)
python3 -B optional/run.py --run-heavy --family C --row 4 --phase c1 --seconds 7200 --binary work/bin/capd-batch1 --output work/capd-row4-c1
```

Use batch1_solver for rows4/14, batch2_solver for5/6, precision10_solver for10,
and row13_radius90 for13. Run C0 and C1 separately with their saved settings.
Row13's historical C0 was computed by row13_solver; its radius-independent
point calculation was reused for the radius90 certificate. Reintegrating one
phase does not by itself reproduce the adopted terminal certificate. The
terminal drivers, saved certificate fields and paper describe the reduction;
producing and auditing an entire replacement C certificate requires both
phase histories and their correct connection. Do not treat a newly printed
terminal PASS as replacement for that evidence.

Full archived CAPD trace audits additionally require large streams and their
receipt/cache tree, which are not bundled; see audit/README.md and its size/hash
inventory. No automatic fetch URL is provided. Finite verification, all seven
tables and negative controls are fully runnable without these optional inputs.
