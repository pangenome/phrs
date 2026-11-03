# Tan et al. 2026 — methods extract (ABSTRACT-LEVEL ONLY; full text still unread)

Tan K-T, Ong RJX, Kee BBR, Yap RKH, Ng AJT, Rajapaksha GK, Garrity-Janger M, Lin Q, Chong JR,
Kyi CT, HPRC, Meyerson M, Li H. "Widespread structural variations at human chromosome ends."
bioRxiv 2026. doi:10.64898/2026.09.03.748401. **PMC13684503, PMID 42922562.** Posted 2026-09-06,
PMC release 2026-10-22, CC BY-NC 4.0. Data: Zenodo 10.5281/zenodo.21807156.

**Status: the full text could not be obtained from this host.** Everything below is supported only by
the abstract and by metadata. See "Retrieval record" at the end - it is complete, so this does not
need to be attempted again from the same machine.

## 1. Dataset (abstract-supported)

"a collection of complete subtelomeric sequences spanning all non-acrocentric chromosome arms, derived
from **860 haploid assemblies** across **six ancestry groups**."

*Not available:* per-population counts, number of individuals, how many ends were retained per
assembly, and the operational definition of "complete" (the input assembly count must not be read as a
count of recovered complete ends).

## 2. Structural-variation method — UNVERIFIED

No aligner, parameters, identity or length thresholds, haplotype-classification rule, or
boundary-placement procedure appear in the abstract. Their phrase "mosaics of blocks shared between
chromosome arms" is a description, not a method.

## 3. Dating / age estimation — PRIORITY ITEM, qualitative only

Two qualitative statements, both in the abstract:

- diversity "is associated with chromosome arm rather than ancestry, suggesting that cross-arm
  paralogy block duplications **predate human population divergence**";
- the blocks "are **broadly conserved across the great apes**".

**No mutation rate, molecular clock, coalescent analysis, calibrated duplication age, or
divergence-time estimate appears in the abstract.** Whether such an analysis exists in the full text
is unknown. Consequently we cannot yet compare their dating approach with our `d = 2*mu*T` bound
(5-12 Myr upper bound, 2-5 Myr at mu = 3e-8), and we must not attribute a rate to them.

## 4. Concerted evolution / homogenisation — partially supported

Abstract states: "most non-acrocentric autosomal arms (**54%, 21 of 39**) carry **multiple haplotypes
differing by up to 100-200 kb**"; divergent haplotypes "differ in gene content, driving gene
copy-number variation between individuals among olfactory receptors and other genes"; and the study
"revealed **rare subtelomeric recombination**".

*Not available:* gene conversion, homogenisation, concerted evolution, or conversion-tract length are
not mentioned in the abstract. "Rare recombination" must not be rewritten as gene conversion.

## 5. Telomere length

Abstract claims the dataset "enables the accurate measurement of telomere length at individual
chromosome ends from long-read data". Algorithm, anchoring filters, and validation are unavailable.

## 6. Citation of Guarracino et al. 2026 — UNVERIFIED

Whether they cite "Concerted evolution and unorthodox recombination of human subtelomeres"
(doi:10.64898/2026.07.10.737660) cannot be determined without the reference list. The abstract
mentions no critique. This is absence of evidence, not evidence of absence.

## 7. Numbers available for comparison

| quantity | their abstract | comparability note |
|---|---|---|
| input assemblies | 860 haploid | not recovered ends; our set is 465 |
| ancestry groups | six | we do not stratify by ancestry |
| arms with multiple haplotypes | 21/39 (54%) non-acrocentric autosomal | our parallel figure is 41 signal-bearing arms of 48 |
| within-arm haplotype differences | up to 100-200 kb | structural length, **not** nucleotide divergence |
| divergence | none quoted | cannot be compared with our d_sub |
| OR copy-number variation | qualitative | our enrichment: 39 families, 187 copies, 8 families hold 142 |
| recombination | "rare" | ours: 2.3% of arms carry a strong candidate |

## Retrieval record (do not repeat from this host)

- `biorxiv.org/.../2026.09.03.748401v1.full.pdf` - HTTP 403 then **429**, with and without a
  browser user agent, both the `v1.full.pdf` and `biorxiv/early/...` routes.
- `api.biorxiv.org/details/biorxiv/<doi>` - metadata only, no abstract body.
- Europe PMC `fullTextXML` for PMC13684503 - returns front matter plus an explicit
  **"Full Text Availability"** section: *"The license terms selected by the author(s) for this
  preprint version do not permit archiving in PMC."* So PMC will never serve it.
- OpenAlex W7211892909 - confirms OA (`oa_status: green`), `pdf_url` = the bioRxiv full PDF above,
  i.e. the only open copy is the one this host is rate-limited from.
- Semantic Scholar API - 429.
- sciety / lifescience.net / biorxiv.science - abstracts only; `.full`/`.full.pdf` routes dead.
- A delegated Codex CLI run also failed to retrieve it, and additionally **failed to write files**
  (`bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted`) - use the `sxgc.*` user lanes
  (danger-full-access) rather than the builtin `codex-exec*` agents on this host.

**Ask:** download the PDF in a browser (institutional access or the bioRxiv page) and drop it at
`~/Downloads/papers/Tan2026_chromosome_ends.pdf`, then sections 2, 3 and 6 can be completed and the
dating comparison in `lit_Tan2026_chromosome_ends.md` finished.
