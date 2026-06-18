# Third-Party Tool License Audit

Date: 2026-06-17
Scope: every external tool currently checked into `tools/per-tool/` or invoked from
`tools/install/`. The audit decides whether each tool can be **vendored**
(copied into this repository's git history and redistributed under the project
license) or must be installed by the user from upstream.

## Summary Matrix

| Tool | Upstream License | SPDX | Vendorable? | Source of truth |
|---|---|---|---|---|
| SEG | NCBI FTP, no LICENSE file (US-gov public-domain convention) | n/a | ⚠️ Gray area — practically OK, strictly unclear | `tools/per-tool/SEG/README` (NCBI), `ftp://ftp.ncbi.nih.gov/pub/seg/seg/` |
| PScore | eLife article supplement, CC-BY 4.0 | `CC-BY-4.0` | ✅ Yes, with attribution | https://elifesciences.org/articles/31486, archive `elife-31486-code2-v2.tgz` |
| PLAAC | MIT (Whitehead Institute, BBRI, UMass) | `MIT` | ✅ Yes | https://github.com/whitehead/plaac, `tools/per-tool/PLAAC/plaac-master/LICENSE.TXT` |
| catGRANULE v2 | MIT (Jonathan Fiorentino, 2025) | `MIT` | ✅ Yes | https://github.com/tartaglialabIIT/catGRANULE2.0 (LICENSE file) |
| IUPred3 | ELTE academic license — non-commercial, **no redistribution** | non-OSI | ❌ No | `tools/per-tool/IUPred3/iupred3/LICENSE`, registration at https://iupred3.elte.hu/ |
| ESpritz | Tosatto lab academic license — non-commercial, **no redistribution** | non-OSI | ❌ No | `tools/per-tool/ESpritz/espritz/LICENSE` |
| DeepCoil | **No LICENSE declared** upstream (defaults to all rights reserved) | none | ❌ No | https://github.com/labstructbioinf/DeepCoil (no LICENSE, no `pyproject.license`, no README notice) |
| AlphaFold / LocalColabFold | LocalColabFold MIT; AlphaFold weights subject to DeepMind terms | mixed | ⚠️ Wrapper yes, weights no — moot (DBs are tens of GB and not vendorable in any case) | https://github.com/YoshitakaMo/localcolabfold |
| CD-HIT | GPL v2 | `GPL-2.0` | ✅ Yes (but GPL is copyleft) | `tools/per-tool/cdhit/license.txt` (the entry has since been removed), https://github.com/weizhongli/cdhit. **Currently unused in this codebase.** |

Legend: ✅ vendorable, ⚠️ gray area, ❌ must be downloaded by the end user.

## Per-Tool Notes

### SEG (NCBI)

NCBI publishes the C source on its FTP server without a LICENSE file. The
bundled `README` and `seg.doc` carry only attribution requirements ("Wootton &
Federhen 1993"). U.S. federal works produced by NCBI staff are conventionally
treated as public domain inside the United States, but the SEG distribution
predates current NCBI source-policy headers and contains no explicit grant.

Pragmatic decision: SEG is vendored across many bioinformatics packages
(`blast+`, `repeatmasker`, distro packages) on the public-domain reading. We
adopt the same posture but record this caveat. If the project ever needs a
strict-audit posture, replace the vendor with a download-on-install script.

### PScore (Vernon et al., eLife 2018)

eLife distributes the article and its supplementary files under
**CC-BY 4.0**. "Source code 2" (`elife-31486-code2-v2.tgz`) is part of those
supplementary materials and inherits CC-BY. Redistribution is allowed with
attribution. We must keep the citation:

> Vernon, R.M., Chong, P.A., Tsang, B., Kim, T.H., Bah, A., Farber, P.,
> Lin, H., Forman-Kay, J.D. (2018). Pi-Pi contacts are an overlooked
> protein feature relevant to phase separation. *eLife* 7:e31486.

### PLAAC (Lancaster, Nowick, Krieger, Lindquist 2014)

Standard MIT license. We may vendor `cli/` and `web/bin/plaac.jar`. We must
keep the `LICENSE.TXT` next to the code.

### catGRANULE v2 (Monti, Fiorentino et al., 2024)

MIT (GitHub upstream `tartaglialabIIT/catGRANULE2.0/LICENSE`, Copyright (c)
2025 Jonathan Fiorentino). Even though the release archive
(`catGRANULE2.0-1.0.0.tar.gz`) does **not** ship a LICENSE file inside the
tarball, the upstream repository does, and the tarball is a snapshot of that
repository. Vendoring is permitted.

This project currently uses the **paper-formula v1** scorer (`catgranule_v1`),
not v2. The v2 tree (~25 MB after extraction) is documented but not on the
runtime path, so a vendor would be carrying weight we do not exercise. The
practical decision is to keep v2 install-on-demand.

### IUPred3 (Erdős, Pajkos, Dosztányi 2021)

ELTE academic license. The LICENSE explicitly states:

> 6. I will not redistribute the software to others. Academic Licensee shall
>    suggest to other interested partners to contact ELTE directly. The
>    program shall not be made available to users other than the Academic
>    Licensee.

This rules out vendoring. Users must accept the click-through at
https://iupred3.elte.hu/ and download `iupred3.tar.gz` themselves. We may
document the install path and provide a wrapper, which we already do at
`tools/wrappers/run_iupred3.sh`.

### ESpritz (Walsh, Martin, Di Domenico, Tosatto 2012)

Tosatto-lab academic license (Padova). LICENSE clause 5:

> Copying Restrictions. You will not sell or otherwise transfer these
> programs or derivatives to any other party, whether with or without
> consideration, for any purpose.

No vendoring. Distribution is via direct request / private link. Same posture
as IUPred3: document, do not ship.

### DeepCoil (Ludwiczak et al. 2019)

The upstream GitHub repo `labstructbioinf/DeepCoil` (`master` branch as of
2026-06-17) has:

- no `LICENSE` / `LICENSE.txt` / `COPYING` file
- no `license` field in `pyproject.toml`
- no license notice in `README.md`

Under default copyright law, this means **all rights reserved**. We may run
the public PyPI release locally, but we may not redistribute the source.
DeepCoil therefore stays out of git: keep `tools/per-tool/DeepCoil/README.md` plus
`tools/install/install_deepcoil_env.sh`, which `pip install`s the upstream
release into an isolated conda env.

### AlphaFold / LocalColabFold

- LocalColabFold wrapper: MIT (Yoshitaka Moriwaki).
- AlphaFold model parameters: CC-BY 4.0 for the weights, but the original
  AlphaFold codebase is Apache-2.0 with a separate weights notice.
- MMseqs2 ColabFold databases: tens of GB, not vendorable on size alone.

The AlphaFold pLDDT feature was evaluated and found not to improve AUC.
Production drops it; the install/storage scaffolding for that experiment
is not part of this public release.

### CD-HIT (Li & Godzik 2006; Fu et al. 2012)

GPL v2. Vendoring is allowed but GPL is copyleft — derived works would be
forced under GPL too, and our project license has not been finalized.
However:

- **Nothing in `src/`, `scripts/`, or `tests/` invokes `cd-hit`.**
- the project originally listed "CD-HIT or equivalent sequence-identity filtering"
  as an unfinished future work item.
- The previous `Tools/cdhit` entry in this repo was a dangling git submodule
  (`160000` gitlink, no `.gitmodules`, modified content). It exists only
  because someone cloned the upstream repo on 2026-05-07 to confirm the
  binary built.

Decision: remove the gitlink entry. If sequence-identity filtering is ever
implemented, the install path can be added back through
`tools/install/install_cdhit.sh`.

## Vendor Policy

Based on the matrix above, the redistribution-safe tools to consider for
inclusion in git are:

1. **PLAAC** (MIT, ~5 MB extracted)
2. **PScore** (CC-BY 4.0, ~50 MB extracted including DBS/; the trained
   reference DB is the bulky part)
3. **SEG** (gray area; ~660 KB of K&R C source)

Tools that **must not** enter git history:

- IUPred3, ESpritz, DeepCoil, AlphaFold weights/DBs, catGRANULE v2 release
  tarball (kept out by size + because we use the v1 reconstruction).

Tools that are removed entirely:

- CD-HIT (unused, dangling submodule).

## Attribution Obligations Carried Forward

If any of the vendorable tools are bundled, the following must remain in the
repository next to the corresponding source tree:

- PLAAC: `LICENSE.TXT` (MIT) and the original copyright header.
- PScore: a `LICENSE` or `NOTICE` carrying the CC-BY 4.0 grant and the
  Vernon et al. citation.
- SEG: the upstream `README` and `seg.doc` (attribution to Wootton &
  Federhen).
- CD-HIT (if ever re-vendored): the GPL v2 `license.txt` plus a notice that
  the combined work falls under GPL.

This file itself is the audit record.
