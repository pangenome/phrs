# Exchange regimes: ratios of genealogical depths, not rates of events

Supersedes the earlier version of this note, which expressed mixing as "one exchange per N generations".
That unit does not exist for these data and produced numbers that contradicted each other. This note
gives the formulation the data actually support.

## Why "generations per exchange" is the wrong unit

In a region of **complete mixing** (the pseudoautosomal regions, with an obligate crossover in PAR1, and the
acrocentric short arms) the coalescence time of two copies there is **62,500 generations**,
i.e. $\sim 2N_e$. Crossovers do not shorten coalescence to one generation: they place both loci in
**one pool**, and a single pool's genealogy is as deep as the population's. So there is no interval
between events to measure, and any statement of the form "one exchange every N generations" implies a
clock that these data do not contain.

## What the data give: a dimensionless ratio

$$R = \frac{T_{\text{coalescence}}(\text{non-homologous pair})}{T_{\text{coalescence}}(\text{alleles in one pool})}$$

$R$ answers "how many times deeper is this pair's common ancestor than a single pool's". $R = 1$ is
complete mixing; larger $R$ means a barrier to mixing.

| regime | $R$ | $T$ (generations) | $T$ (Myr) | mixing, $1/(R-1)$ |
|---|---|---|---|---|
| complete mixing (acrocentric short arms, PAR1, PAR2) | **1.00** | 62,500 | 1.25 | **∞: one pool** |
| **telomeric end, 6-10 kb** | **1.82** | 113,750 | 2.28 | **1.22** |
| mid compartment, 20-100 kb | 3.44 | 215,000 | 4.30 | 0.41 |
| average block (length-weighted, all records) | 5.93 | 370,625 | 7.41 | 0.20 |
| far interior, 200-500 kb | 6.63 | 414,375 | 8.29 | 0.18 |
| weakest blocks (ratio 44) | 44.0 | 2,750,000 | 55.0 | 0.023 |

Data table: `submission/supporting_material/ED_exchange_regimes_table.tsv`; paste-ready LaTeX:
`response/patches/ED_table_exchange_regimes.tex`.

## The model-consistent rate column

Under coalescence with migration between two compartments, $T = 2N_e + 1/(2m)$, so

$$4N_e m = \frac{1}{R-1}$$

which is the last column: mixing **relative to complete mixing**, where complete mixing is infinite
(no barrier) and a fully isolated pair would be zero. This is the quantity worth quoting; it needs no
per-generation event count.

## Why $R$ is the robust quantity and the times are not

$R$ is a ratio of two divergences measured identically, in the same GC-rich compartment, so the
substitution rate **cancels**. That matters for two reasons we already know about:

1. **GC-biased gene conversion.** This compartment is 62% GC at the tip; gBGC raises the substitution
   rate toward GC, so $\mu_{\text{eff}} > \mu$ and every absolute $T = d/(2\mu)$ is an **overestimate**.
2. **Elevated paternal-germline mutation.** Same direction, same overestimate.

Both inflate numerator and denominator together, so $R$ survives them while the absolute times do not.
**Report $R$ and $1/(R-1)$; treat the generations column as an upper bound carrying model dependence.**

## What this corrects in earlier statements

| earlier statement | status |
|---|---|
| "one exchange per ~740,000 generations per locus pair" | **withdrawn** — conflates a per-pair rate with the ratio anchored on complete mixing, and implies an event clock that does not exist |
| "one inter-arm exchange somewhere every ~920 generations" | **withdrawn** — a system-wide figure computed from the average block and presented alongside per-pair figures as if comparable |
| "17% of the obligate rate" for the average block (now: of complete mixing) | **retained in substance**, as $1/R = 1/5.93 = 0.17$ |
| "the most cohesive blocks approach the obligate rate" | **retained**: $R = 1.82 \Rightarrow$ mixing 1.22 relative to complete mixing, i.e. near-pool |
| "the ratio is not linear in rate" | **retained and sharpened**: the mapping is $1/(R-1)$, hyperbolic, and undefined-favouring at the complete-mixing end rather than a multiple |

## Statistical caution on the gradient

The gradient must be quoted from one statistic throughout. The per-arm median ratio (one vote per
chromosome arm) is **1.82 at 6-10 kb, 3.44 at 20-100 kb, 6.63 at 200-500 kb**; the length-weighted
cohort figure of 5.93 is a different quantity and must not be placed on the same axis. Mixing the two
(an earlier error) produced an apparent contradiction between the tip ratio and the average rate.
