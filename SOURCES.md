# Sources and numerical excerpts

The candidate Euler coordinates are the first 20 rows of
`euler_100_4860.txt` from Hristov et al.'s 2025 catalogue,
<https://db2.fmi.uni-sofia.bg/3bodysym4860/>. The corresponding free-coordinate
records come from `free_100_4860.txt` at the same source. The article's
reference is <https://arxiv.org/abs/2503.00432v1>. Whole-source SHA256 values
are, respectively, `e93d789a7c1dff2931ff2772754d4072d464704ec35cb9f466274bfc5d2977d3`
and `74a23a8392a4ec75b3e855973a6d4cfcacc90de84784688f516a366d6185dccb`.
The packet retains the original decimal tokens needed for its certificates.

Earlier matching records are from <https://db2.fmi.uni-sofia.bg/3bodyfree/>,
`sol_80.txt` and `sd_80.txt`. The saved comparison records identify their
original row numbers and tokens; positive matches are checked afresh against
those tokens. Reported full-table sizes (24582 and 236) and three unconfirmed
labels are saved-record information, not a new full-table search. Catalogue
authors retain the credit for the numerical discovery and original data.

`provenance.json` records original size/hash, exported size/hash and fields
removed for local-path privacy. The release production record retains the
internal source-to-export mapping; no private path is needed to run the
public checks. Sanitized JSON is a derivative, not a byte-identical original.
Hash integrity binds the distributed data, and is not itself a mathematical
proof of an enclosure.

The fixed CAPD revision for the recorded audits is
`37057e5899bccaa6cd5701057b7e6ea5e232ef8b`; the observed environment and
source fingerprints are under optional/capd/ and audit/. The Arb environment
records python-flint 0.9.0 / FLINT 3.6.0. The handwritten MPFR records identify
MPFR 4.2.2 / GMP 6.3.0. Code dependencies are external and not fetched by
packet checking. Consult the article bibliography for the software papers
and the earlier geometric and rigorous numerical work.
