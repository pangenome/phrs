# M5 — Matched SD specificity control (Referee #1 major point 5; Referee #3 major point 3)

Status: IN PROGRESS (skeleton created; appending after each major step)

Question: Is the sequence-to-contact relationship a property of SUBTELOMERES, or of
high-identity segmental duplications (SDs) generally?

## Log (append-only)

- [step 0] Skeleton notes file created.

### Inventory so far (step 1, reconnaissance)
- `PEER_REVIEW.md` R1 major 5 (verbatim ask read); tracker row `R1-M5 / R3-M3` = "matched
  non-subtelomeric duplication control (matched identity/length/copy) · L0".
- Baseline number to beat: CHM13 Hi-C arm-level Mantel rho = 0.656; HG002 Hi-C 0.657
  (`end-to-end-report/report/05_hic_validation.md`).
- Hi-C machinery: `.mcool` per sample present under
  `/moosefs/guarracino/HPRCv2/PHR_III/HiC/<sample>/hicpro_output/cool/<sample>.mcool`
  (CHM13, HG002, HG00658, HG02148, HG02559, NA19036); HG002 Pore-C and CiFi mcools under
  `porec_output/porec.mcool` (1.17 GB) and `cifi_output/cifi.mcool` (141 MB).
- Genome-wide CHM13 window table on disk: `data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz`
  (31,202 windows; per-window cross-chromosome haplotype counts + avg identity). This is the
  only genome-wide, chromosome-resolved sharing substrate on disk.
- Genome-wide CHM13-vs-haplotype PAFs exist (`/moosefs/guarracino/HPRCv2/aln-ba8599f.../`,
  465 files) but target contigs are accessions (e.g. `HG00097#1#CM094072.1`) with no
  contig->chromosome map on disk for the 465 haplotypes => not directly usable.
- No segmental-duplication annotation (UCSC hprc2annot / Jeong et al. SD track) is present
  anywhere under /moosefs/guarracino/HPRCv2 (searched `*segdup*`, `*sd.bed*`, `*chm13.dv*`:
  zero hits). SD control blocks must therefore be derived de novo from alignment data.
- tools: `wfmash` present (/home/erikg/bin/wfmash); `odgi`, `vg`, `minimap2`, `mashmap`
  NOT available; `samtools` + `bedtools` present; python lacks h5py/cooler globally,
  scratch venv /tmp/m5venv has h5py 3.16.0 (mcool is HDF5-readable).

### Machinery gate (step 2): exact reproduction of the published contact extraction
- Contact-contract: `nan-sum of `clr.matrix(balance=True).fetch(A,B)` over the bin-pair block,
  divided by the number of non-NaN bin pairs (`hic_contact_norm`), cooler 0.10.4 = the version
  that wrote the mcools (mcool header `generated-by: cooler-0.10.4`).
- GATE G2: recomputed `hic_contact_norm` from `CHM13.mcool::resolutions/50000` for the first 60
  rows of `data/human_CHM13_hic_50000bp_seqlevel.tsv` -> **60/60 bit-exact, max abs diff 0.0**.
- GATE G1: Spearman(jaccard, hic_contact_norm) over all 652 published CHM13 PHR pairs =
  **0.71625** (published 0.716) — reproduced from the published table.
- GATE G3 (documented non-reproduction): the report's arm-level Mantel rho = 0.656 (n=38 arms)
  is NOT reproduced by aggregating the published pair table to arm pairs
  (+0.7499 over the 652 arms pairs at 41 arms, similarity = 1 - dist from
  `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`; variant aggregations 0.725-0.750). The paper's
  arm contact matrix is built inside `analyze_hic_communities.py` (matrix-pixel sums), not from
  the exported pair table. For M5 I therefore use the **pair-level Spearman (exactly
  reproducible, = 0.716)** as the primary comparand for every class, and report an arm-level
  Mantel computed identically for all classes as a secondary statistic.

### Class construction (step 3) — from the only genome-wide sharing substrate on disk
Source: `data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz`, 31,202 rows (CHM13 100 kb windows
aligned against 25272 haplotype contigs at >=98% identity, xm5). 625 chrY rows dropped (the Hi-C
coordinate space is the masked-Y CHM13 fasta, which has no chrY; the UCE mcools have no chrY
either). 17 rows dropped as shorter than 100 kb (chrM). Terminal windows of the 100 kb grid are
anchored at the chromosome end so every unit is exactly 100 kb (length-matched by construction).

Class definitions (all from that one table, so identity/copy-number are estimated identically):
- `T`     : inside a `data/chm13.phrs.bed` interval with >=1 cross-chromosome partner carrying
            >=2 haplotype copies  -> **89 windows / 23 chromosomes**
- `P`     : >1 Mb from the telomere, <=1 Mb from a `data/chm13-annotations.bed` centromere,
            cross-chromosome duplicated -> **470**
- `I`     : >1 Mb from telomere and centromere, cross-chromosome duplicated -> **804**
- `Utelo` : in the 1 Mb telomere band, not a PHR, no cross-chromosome duplicate -> **357**
- `Ucen`  : <=1 Mb of a centromere, no cross-chromosome duplicate -> **982** (position-matched
            unique control: same geography as `P`, no duplication)
- `Ugen`  : unique, >1 Mb from telomere and centromere -> **27,798**
- excluded-but-retained-for-inspection: `telo_band_shared_not_PHR` (56), `T_nosharing` (4).

Window set analysed: T 89, P 470, I 470 (of 804), Utelo 357, Ucen 470 (of 982), Ugen 470 (of
27,798) = 2326 windows x 100 kb = 233 Mb, sampled round-robin over chromosomes (seed 20260518).
**Matching imbalance to be reported, not hidden**: T windows are far more multi-copy than the
SD controls — median distinct partner chromosomes 14 (T) vs 2 (P) vs 5 (I); median haplotype
support 2014 (T) vs 458 (P) vs 340 (I). Copy-number matching is therefore restrictive and is
reported with its balance table.

### Alignment regimes for the sequence axis (step 4)
- Regime A `p95`: `wfmash -p 95 -L -f -M` — the identity regime of the PHR call's all-vs-all,
  with **all** mappings retained. This matters: default wfmash filtering keeps one best mapping
  per segment, which would penalise exactly the high-copy windows. Pilot (30 windows):
  default 36 records/142 kb, `-M` 128/318 kb, `-f` 229/505 kb, `-f -M` 351/732 kb.
- Regime B `p90_pggb`: `wfmash -p 90 -s 1000 -l 3000 -n 10 -k 15 -L` — the exact wfmash
  invocation pggb uses to build the graph the published odgi Jaccard is computed on (filtered,
  <=10 mappings/segment). Second alignment regime for the sequence axis.

### Pipeline gotchas found (step 5) — recorded so they are not repeated
- `wfmash` groups mappings by the last `#` of the sequence name (`-Y '#'` default) and then skips
  self-group mappings. PanSN-style window names put *every* window in one group and silently
  suppress all paralogous mappings (measured: 0 mappings for a set of 5 sequences including
  exact duplicates). Window FASTA names are therefore `chr:start-end` (no `#`).
- `wfmash` default filtering keeps ~1 mapping per query segment (pilot: 36 records vs 351 with
  `-f -M`), which penalises high-copy windows specifically. Production sequence regimes use
  `-f -M`.
- `wfmash` writes its temporary mapping file into the *current directory* unless `-B` is given;
  all invocations pass `-B <workdir>` so nothing lands in the repo.
- PAFs are written to `.part` and renamed, so an interrupted run can never masquerade as a
  complete alignment (an early jaccard pass on a partially written PAF produced 9,539 pairs
  instead of the final count — caught and discarded).
- The independent per-base recomputation in the oracle initially used a *global* coverage array
  instead of a per-pair one and falsely flagged the estimator; fixed, and all 4 semantic checks
  plus all per-base equality checks now pass.
