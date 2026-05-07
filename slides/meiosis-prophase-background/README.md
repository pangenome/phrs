# Meiotic Telomere Biology Review Deck

Standalone background deck for meiosis prophase I telomere/chromosome organization.
This material is intentionally separate from the current draft presentation so images
can be inspected and manually reused later.

## Outputs

- `meiotic_telomere_biology_review.typ` - reproducible Typst deck source.
- `meiotic_telomere_biology_review.pdf` - rendered review deck.
- `assets/` - downloaded figure image assets.
- `image_manifest.tsv` - per-image source URL, article, figure/panel, license,
  attribution, local file, and notes.

Build command:

```bash
typst compile meiotic_telomere_biology_review.typ meiotic_telomere_biology_review.pdf
```

Run the command from `slides/meiosis-prophase-background/`.

## Deck Coverage

The rendered deck covers the requested biology:

- Overview of prophase I stages using the Frontiers/PMC Figure 1 staging schematic.
- Zygotene human oocyte and spermatocyte examples from Cheng et al. 2009 and
  Gruhn et al. 2013.
- Zygotene bouquet/telomere clustering and telomere-adjacent synapsis initiation
  using the Blokhina et al. 2019 PLOS Genetics zebrafish figure.
- Pachytene full synapsis in human oocyte/spermatocyte side-by-side from Gruhn
  et al. 2013.
- Diplotene/desynapsis/chiasmata context using the open-access Bisig et al. 2012
  mouse figure as a substitute for paywalled human oocyte material.
- Telomere-NE tether mechanism using Dunce et al. 2018 TERB1/TERB2/MAJIN figures
  and Chen et al. 2021 SUN1-SPDYA Figure 5.
- Final citation/license slide with open-license and skipped/source-only notes.

## Downloaded Assets

Downloaded and used in the deck under open licenses:

- `assets/sanchez_saez_2021_frontiers_fig1.png` - Wang and Pepling 2021 Frontiers
  Figure 1. The brief referred to Sanchez-Saez/Bolcun-Filas, but the supplied
  URL resolves to Wang and Pepling 2021; the deck and manifest record the
  resolved source.
- `assets/cheng_2009_plosgenet_fig1_large.png` - Cheng et al. 2009 PLOS Genetics
  Figure 1, human fetal oocyte stages.
- `assets/gruhn_2013_plosone_fig7_large.png` - Gruhn et al. 2013 PLOS ONE Figure 7,
  human zygotene spermatocyte/oocyte examples.
- `assets/gruhn_2013_plosone_fig1_large.png` - Gruhn et al. 2013 PLOS ONE Figure 1,
  human pachytene oocyte and spermatocyte.
- `assets/blokhina_2019_plosgenet_fig1_large.png` - Blokhina et al. 2019 PLOS
  Genetics Figure 1, zebrafish bouquet progression.
- `assets/bisig_2012_plosgenet_fig1_large.png` - Bisig et al. 2012 PLOS Genetics
  Figure 1, pachytene-to-diplotene SC/centromere progression.
- `assets/dunce_2018_natcomm_fig1.png` - Dunce et al. 2018 Nature Communications
  Figure 1, TERB1/TERB2/MAJIN tether schematic.
- `assets/dunce_2018_natcomm_fig8.png` - Dunce et al. 2018 Nature Communications
  Figure 8, SIM telomere attachment plate imaging.
- `assets/mu_2021_natcomm_fig5.png` - Chen et al. 2021 Nature Communications
  Figure 5, SUN1-SPDYA ring architecture. The brief listed this as Mu et al.;
  the Nature article page citation is Chen et al. 2021.

Downloaded with non-CC/unclear rights status:

- `assets/liebe_2004_mbc_fig1_pmc_review_use.jpg` - Liebe et al. 2004 MBoC Figure 1.
  This is a useful classic mouse 3D telomere FISH reference for zygotene bouquet
  versus pachytene telomere positioning. It is free in PMC, but automated review
  did not find an explicit CC reuse license. The manifest records that status;
  the slide itself cites the paper normally and keeps the biology prominent.

## Skipped Sources

Skipped rather than silently using paywalled or unclear-license content:

- Mytlis et al. bioRxiv/Science bouquet figures: the bioRxiv page returned HTTP
  403 during automated retrieval in this environment, and the Science article is
  paywalled. The Blokhina et al. PLOS Genetics figure was used instead as the
  open CC-BY bouquet visual.
- Lenzi et al. 2005 AJHG human oocyte diplotene/chiasmata figures: paywalled /
  unclear reusable image access. Bisig et al. 2012 PLOS Genetics was used as the
  open-access diplotene/desynapsis substitute.
- Shibuya et al. 2015 Cell MAJIN/shelterin figures: paywalled. Dunce et al. 2018
  and Chen et al. 2021 Nature Communications figures were used for mechanism.
- Morimoto et al. 2012 JCB KASH5/dynein figures: not used because image reuse
  rights were not verified for this deck. The mechanism is noted textually.
- Link et al. 2014 PLOS Genetics SUN1/LINC figures: open-access and useful, but
  not downloaded for this compact deck because Dunce et al. 2018 and Chen et al.
  2021 cover the mechanism more directly.
- Kouznetsova et al. 2005 JCS diplotene/oocyte figures: not downloaded because
  Bisig et al. 2012 provides the compact open-access diplotene/SC persistence
  visual used here.
- Scherthan et al. 2000 MBoC human/mouse pachytene NE source: skipped because
  the PMC URL supplied in the brief did not resolve to the expected telomere
  article during automated retrieval.
- Enguita-Marruedo et al. 2018 Chromosoma: the PMC ID in the brief did not
  resolve to the expected article during automated retrieval, so no figure was
  used from that source.
- Fan et al. 2021 PLOS Genetics timing schematic: optional human-specific timing
  source. Not downloaded because the overview slide already covers staging/timing
  and the deck prioritized telomere/FISH/prophase images.

## Notes For Reuse

- Treat the PLOS, Frontiers, and Nature Communications figures as reusable with
  attribution under their recorded CC licenses.
- Treat `liebe_2004_mbc_fig1_pmc_review_use.jpg` as useful internal review material,
  not as a cleared external reuse asset.
- The deck intentionally uses full downloaded figures rather than destructive
  crops. Panel callouts in slide text identify the requested panels, while the
  original figure context remains intact for review.
