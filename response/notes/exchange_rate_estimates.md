# What exchange rate and transfer amount do these divergences imply?

Derived from the cohort numbers: allelic d = 0.150%, paralog d = 0.890% (ratio 5.93x), mu = 1.2e-8,
20 yr/generation, T = d/(2*mu).

## The bounds are a FLOOR and a CEILING, not a rate range

**Correction to a natural first guess:** PAR1/PAR2 and the acrocentric p arms are not an "upper bound
that drives divergence to zero". They recombine obligately - every individual, every generation, the
entire substrate - and their divergence still sits at **0.12-0.15%**, i.e. exactly the allelic level.
The reason is that obligate recombination makes two loci behave as alleles *within* an individual, but
the genealogy of those alleles in the population still coalesces at ~2*N_e. So:

- **Floor = the allelic level (0.15%) = a single freely recombining pool.** Cannot go below this by
  any amount of exchange.
- **Ceiling = no exchange since duplication** (percent-level, and unbounded in principle).
- Our blocks at 0.890% sit **5.9x above the floor**.

So the meaningful quantity is not "rate" on an absolute scale but **exchange intensity relative to
obligate**, given by T_allelic / T_paralog:

| block class | ratio to allelic | intensity vs obligate |
|---|---|---|
| acrocentric p, PAR1 | 1.00 | **100%** (obligate; at the floor) |
| PAR2 | 1.25 | 80% |
| near-allelic (e.g. 19p 1.33, 8p 1.23) | 1.2-1.5 | **67-83%** |
| cohesive tier (C3/C11/C1) | 3.3-3.5 | 29% |
| **average block** | **5.9** | **17%** |
| weak (C6, C2) | 20-44 | **2-5%** |

## Implied absolute rate (with assumptions stated)

Inverting the coalescent, if exchange dominates coalescence between a pair of copies:

    T_allelic = 62,500 generations (1.25 Myr)      [the floor]
    T_paralog = 370,833 generations (7.42 Myr)     [the average block]

    exchange rate per locus pair  c = 1/(2*T_paralog) = 1.35e-6 per generation
                                  = 1 event per ~742,000 generations (14.8 Myr) per pair

Over the **805 heterologous arm pairs**:

    805 x 1.35e-6 = 1.1e-3 events per generation genome-wide
    = one inter-arm exchange somewhere in the subtelomere every ~920 generations (~18,000 years)

**Amount transferred:** the tract size is independently measured at **~50 kb** - the longest shared
block is 49.5-54.0 kb in *every* heterologous pair examined, and the 9q/3q module itself is 54.1 kb.
So the flux is ~54 bp per generation summed over the whole system, and a given 500 kb window is fully
replaced once per ~1.7e8 generations. Individually: **short tracts, transferred rarely per locus,
often genome-wide** - which is exactly the qualitative picture that motivated this analysis.

## Independent cross-check from population frequency

Discrete transfer modules are observed at population frequency ~2% (our 2.3% of arms; Tan et al. report
2.5% and 0.58% at 4q/10q). Under a neutral mutation-drift equilibrium, p = 4*N_e*r, so with N_e = 3e4:

    r = p/(4 N_e) = 1.7e-7 per generation per locus  (1 per ~6e6 generations)

versus 1.35e-6 from the divergence inversion. **Agreement within ~8x** - i.e. the same order of
magnitude from two completely independent routes (sequence divergence vs population frequency). That
is reassuring given how crude both are.

## Assumptions and their failure modes

1. Exchange dominates coalescence for these copies; if drift matters, c is overestimated.
2. Tracts (50 kb) are shorter than the window (500 kb), so a tract does not reset the whole locus.
3. mu = 1.2e-8 and constant; if mu_subtelo is 2-4x higher, all times shrink and c rises proportionally.
4. `2*N_e = 62,500` generations is inferred from the allelic level itself, so the floor is self-consistent
   but not independently measured here.
5. The population-frequency check assumes neutrality and no ascertainment bias toward common modules.

## Why this matters for the argument

It converts the divergence result into a statement about mechanism that can be tested elsewhere:
**the average subtelomeric block experiences inter-arm exchange at about one-sixth of the obligate
rate, with a measured tract size of ~50 kb, and a handful of blocks (ratio ~1.2) approach the obligate
regime.** That is the sentence the paper can defend; the absolute per-generation rate should be quoted
only with the assumptions attached.

## Is the ratio a multiple? No - it is a reciprocal, and the linear model fails measurably

The natural guess is that R (paralog/allelic) maps to recombination rate by some multiple, or perhaps
logarithmically. **Neither. The relation is R = T_paralog / T_allelic, and rates go as 1/R.**

Why, and why the linear reading is ruled out by our own data:

| R | T_paralog (generations) | c = 1/(2T) per generation | intensity vs full pooling (1/R) |
|---|---|---|---|
| 1.00 | 62,500 | 8.0e-6 | 100% |
| 1.20 | 75,000 | 6.7e-6 | 83% |
| 2.00 | 125,000 | 4.0e-6 | 50% |
| 3.50 | 218,750 | 2.3e-6 | 29% |
| **5.93** | **370,625** | **1.35e-6** | **17%** |
| 12.0 | 750,000 | 6.7e-7 | 8% |
| 44.0 | 2,750,000 | 1.8e-7 | 2% |

**The escape model contradicts itself at the obligate extreme.** If exchange reset a pair at rate c per
generation, PAR1/PAR2 - which cross over essentially every male meiosis, c ~ 1 - would give
T = 1/(2c) = 0.5 generations and divergence ~0.000001%. Observed PAR2 (chrXq-chrYq) is **0.124%**, which
is T = 51,700 generations, i.e. **exactly the allelic floor**, not zero. So exchange does not drive
divergence to zero; it drives it to the depth of a single population gene pool, where it stops.

**Consequences:**

1. **R is a ratio of genealogical depths, not a rate.** "5.93x" means the typical block's genealogy is
   ~6x deeper than that of two alleles in one pool. That is the honest statement.
2. **The safe rate-like quantity is 1/R** - the exchange intensity as a fraction of what full pooling
   would give. R = 5.93 -> 17%; R = 1.2 -> 83%; R = 44 -> 2%. This needs no model of how an event
   resets a tract, because it is defined relative to the obligate case we can measure.
3. **Equal steps in 1/R are equal steps in intensity; equal steps in R are not.** The mapping is
   hyperbolic. Reporting "N times more recombination" from a ratio of divergences is therefore wrong
   unless N is read as a change in 1/R.
4. **Absolute per-generation rates remain model-dependent.** c = 1/(2T) gives 1.35e-6 per generation for
   the average block, but since that same model fails by six orders of magnitude at the PAR end, the
   number should be quoted only with the caveat attached - unlike 1/R, which is model-free.
5. **The 8x agreement between c = 1/(2T) and the population-frequency estimate (1.7e-7) should therefore
   be treated as suggestive, not as validation**, because both are built on the same kind of model.

**Bottom line for the manuscript:** present R and 1/R. Do not present a per-generation recombination rate
as if the ratio gave one by multiplication.
