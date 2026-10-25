# Saved CAPD audits

For rows 4, 5, 6, 10, 13, 14, each directory contains the available original
coverage/terminal/provenance summaries with local path metadata removed,
and the integer/dyadic arithmetic and audit programs. These are records of
previous source and saved-trace checks, not checks rerun by check_all.py.
The adopted normal-call scope and the general shrinking-rounding counterexample
are explained in VERIFICATION-RECORD.md and the article.

The archived audit entrypoints reference an auxiliary receipt/cache tree and
large compressed streams. They are provided for inspecting the mathematical
checks, not as standalone commands supplied with all inputs. Local hard-coded
RAW paths have been relocated to `external-evidence` next to each program.
Before a full replay, the receipt tree and streams must be transferred with
hash verification and their relative paths bound explicitly. Do not replace
missing records with fabricated summaries or silently read a research checkout.

`large-evidence-inventory.json` lists saved size/hash bindings across compressed
and uncompressed layers; its sum is not a download size because layers can
represent the same data. The giant streams are outside this Git package.
No download service, public asset or URL is asserted. Arranging hash-checked
archival provision is a remaining task for anyone requesting a full-trace
replay; the finite default checks and all seven table checks require none of
these streams. No full trace audit, source-method audit or heavy integration
was rerun while assembling version 1.0.0.
