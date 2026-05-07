// Meiosis prophase I telomere biology review deck
// Build: typst compile meiotic_telomere_biology_review.typ meiotic_telomere_biology_review.pdf

#set page(
  width: 13.33in,
  height: 7.5in,
  margin: (x: 0.38in, y: 0.30in),
)
#set text(size: 11pt, lang: "en")
#set par(justify: false, leading: 0.62em)

#let navy = rgb("#17365d")
#let teal = rgb("#0f766e")
#let amber = rgb("#b45309")
#let red = rgb("#b91c1c")
#let pale = rgb("#f7fafc")
#let line = rgb("#cbd5e1")
#let muted = rgb("#64748b")

#let header(num, stage) = block(
  width: 100%,
  fill: rgb("#e8f0fa"),
  stroke: 0.6pt + rgb("#c7d8ee"),
  radius: 2pt,
  inset: (x: 0.09in, y: 0.045in),
)[#text(size: 8pt, fill: navy)[*#num* #h(0.18in) #stage]]

#let point_box(body, fill: rgb("#eef9f6"), stroke: teal) = block(
  width: 100%,
  fill: fill,
  stroke: (left: 3pt + stroke),
  radius: 3pt,
  inset: (x: 0.10in, y: 0.08in),
)[#text(size: 12pt, weight: "bold", fill: stroke)[#body]]

#let source_line(body) = text(size: 7pt, fill: muted)[#body]

#let image_slide(
  num,
  stage,
  title,
  img,
  point,
  bullets,
  source,
  fig_width: 64%,
  img_height: 5.82in,
  rights: none,
) = {
  header(num, stage)
  v(0.05in)
  grid(
    columns: (fig_width, 100% - fig_width),
    gutter: 0.24in,
    [
      #if rights != none {
        block(
          width: 100%,
          fill: rgb("#fff7ed"),
          stroke: 0.7pt + rgb("#fed7aa"),
          radius: 2pt,
          inset: (x: 0.07in, y: 0.045in),
        )[#text(size: 8pt, weight: "bold", fill: amber)[#rights]]
        v(0.035in)
      }
      #image(img, fit: "contain", width: 100%, height: img_height)
      #v(0.035in)
      #source_line(source)
    ],
    [
      #text(size: 24pt, weight: "bold", fill: navy, hyphenate: false)[#title]
      #v(0.13in)
      #point_box(point)
      #v(0.13in)
      #text(size: 10.6pt)[#bullets]
    ],
  )
}

#let text_slide(num, stage, title, body) = {
  header(num, stage)
  v(0.13in)
  text(size: 26pt, weight: "bold", fill: navy)[#title]
  v(0.16in)
  text(size: 10.2pt)[#body]
}

#image_slide(
  "01",
  "Overview",
  "Prophase I staging map",
  "assets/sanchez_saez_2021_frontiers_fig1.png",
  [Use this as the orientation slide: leptotene, zygotene, pachytene, and diplotene are defined by SC assembly/disassembly and oocyte timing.],
  [
    - Zygotene: homologs begin synapsis while chromosome ends cluster in the bouquet.
    - Pachytene: synapsis is complete; telomeres are redistributed but remain NE-attached.
    - Diplotene: SC disassembles; chiasmata remain as physical homolog links.
    - Female meiosis reaches diplotene during fetal life and arrests there; male meiosis is continuous after puberty.
  ],
  [Wang and Pepling 2021, Front. Cell Dev. Biol. Fig. 1. CC BY. Downloaded from Frontiers image provider.],
  fig_width: 61%,
)

#pagebreak()

#image_slide(
  "02",
  "Zygotene / Human Oocyte",
  "Human fetal oocyte zygotene",
  "assets/cheng_2009_plosgenet_fig1_large.png",
  [Panel B is the zygotene human fetal oocyte: axial elements are organized, but synapsis is not yet complete.],
  [
    - SYCP3 marks chromosome axes in red; MLH1 is green; CREST/centromeres are blue.
    - Panel C on the same figure provides the pachytene comparison: longer, cleaner axes and mature recombination foci.
    - This is a useful human-oocyte anchor before moving to clearer telomere-specific model-organism images.
  ],
  [Cheng et al. 2009, PLOS Genetics Fig. 1. CC BY. Direct large PLOS image.],
  fig_width: 70%,
  img_height: 5.35in,
)

#pagebreak()

#image_slide(
  "03",
  "Zygotene / Human Spermatocyte",
  "Human zygotene synapsis begins at chromosome ends",
  "assets/gruhn_2013_plosone_fig7_large.png",
  [Human zygotene spreads show partial synapsis and terminally initiated pairing; panel A is the spermatocyte anchor.],
  [
    - Panel A: male zygotene with SYCP3 in red and CREST in blue.
    - Panels B/C: oocyte examples with SYCP1 added, showing the transition from unsynapsed axes to synapsed regions.
    - The figure is useful because it pairs microscopy with schematic insets that make partial synapsis readable.
  ],
  [Gruhn et al. 2013, PLOS ONE Fig. 7. CC BY 4.0. Direct large PLOS image.],
  fig_width: 62%,
)

#pagebreak()

#image_slide(
  "04",
  "Zygotene Bouquet",
  "Telomere bouquet as a pairing hub",
  "assets/blokhina_2019_plosgenet_fig1_large.png",
  [The zebrafish bouquet is the clearest CC-BY visual of telomere clustering and telomere-adjacent synapsis initiation.],
  [
    - Columns move from leptotene through zygotene into pachytene.
    - Magenta telomere signal is concentrated with short Sycp3 axes early, then spreads as synapsis progresses.
    - Rows E-G zoom the telomere-adjacent zone where synapsis starts and telomere associations occur.
  ],
  [Blokhina et al. 2019, PLOS Genetics Fig. 1. CC BY 4.0. Direct large PLOS image.],
  fig_width: 68%,
)

#pagebreak()

#image_slide(
  "05",
  "Mouse 3D Telomere FISH",
  "Classic 3D telomere FISH: bouquet to pachytene context",
  "assets/liebe_2004_mbc_fig1_pmc_review_use.jpg",
  [Structurally preserved spermatocyte nuclei show telomere redistribution across the zygotene-to-pachytene transition.],
  [
    - The brief flags panel A(iii) for wild-type zygotene bouquet topology and A(iv) for the pachytene comparison.
    - This figure is lower-resolution than the PLOS/Nature assets, but it is one of the clearest classic telomere-FISH references for 3D-preserved nuclei.
    - Keep it available as review context when deciding whether the final talk needs a direct bouquet-to-pachytene telomere positioning panel.
  ],
  [Liebe et al. 2004, Mol. Biol. Cell Fig. 1. PMC/free-to-read source; details in manifest.],
  fig_width: 50%,
  img_height: 4.85in,
)

#pagebreak()

#image_slide(
  "06",
  "Pachytene / Human",
  "Human pachytene: full synapsis in oocyte and spermatocyte",
  "assets/gruhn_2013_plosone_fig1_large.png",
  [Pachytene is the fully synapsed state: continuous SCs, mature recombination foci, and paired homologs across both sexes.],
  [
    - Panel A: human fetal oocyte pachytene.
    - Panel B: human adult spermatocyte pachytene.
    - This is the strongest human side-by-side image for the pachytene cellular background slide.
  ],
  [Gruhn et al. 2013, PLOS ONE Fig. 1. CC BY 4.0. Direct large PLOS image.],
  fig_width: 62%,
)

#pagebreak()

#image_slide(
  "07",
  "Diplotene",
  "Diplotene: SC disassembly, centromere persistence, chiasmata context",
  "assets/bisig_2012_plosgenet_fig1_large.png",
  [As prophase exits pachytene, the SC disassembles along arms; homologs remain connected at chiasmata and centromeric SC components can persist.],
  [
    - The figure follows pachytene-to-diplotene progression in mouse spermatocytes.
    - It is the open-access substitute for paywalled human diplotene oocyte figures.
    - Use this slide to flag the biological transition: telomere-led prophase architecture is resolving as homologs prepare for division.
  ],
  [Bisig et al. 2012, PLOS Genetics Fig. 1. CC BY. Direct large PLOS image.],
  fig_width: 64%,
  img_height: 5.72in,
)

#pagebreak()

#image_slide(
  "08",
  "Mechanism / TERB-MAJIN",
  "TERB1-TERB2-MAJIN attaches meiotic telomeres to the NE",
  "assets/dunce_2018_natcomm_fig1.png",
  [The core tether bridges telomeric shelterin to the inner nuclear membrane through MAJIN, TERB2, and TERB1.],
  [
    - This mechanism explains why telomeres stay NE-associated after the bouquet stage.
    - Use the schematic to connect cell-scale bouquet images to molecular attachment machinery.
    - Relevant pathway: telomere shelterin -> TERB1/TERB2/MAJIN -> SUN1/LINC/KASH5 -> cytoskeletal force.
  ],
  [Dunce et al. 2018, Nature Communications Fig. 1. CC BY 4.0. Springer Nature image.],
  fig_width: 66%,
  img_height: 4.75in,
)

#pagebreak()

#image_slide(
  "09",
  "Mechanism / Attachment Plates",
  "Pachytene telomere attachment plates remain active",
  "assets/dunce_2018_natcomm_fig8.png",
  [SIM imaging links telomeric DNA, TRF1, TERB2, and MAJIN at pachytene NE attachment plates.],
  [
    - This is the direct visual support for "pachytene telomeres are dispersed but still anchored".
    - It complements the human pachytene spread: the human image shows full synapsis, this mouse SIM image shows the attachment machinery.
    - The attachment is functional, not merely a leftover from zygotene bouquet clustering.
  ],
  [Dunce et al. 2018, Nature Communications Fig. 8. CC BY 4.0. Springer Nature image.],
  fig_width: 59%,
)

#pagebreak()

#image_slide(
  "10",
  "Mechanism / SUN1",
  "SUN1 ring architecture at meiotic telomeres",
  "assets/mu_2021_natcomm_fig5.png",
  [SUN1-associated architecture provides the LINC-side context for force transmission and telomere-led chromosome movement.],
  [
    - Figure 5 is useful if the mechanism slide needs a SUN1/LINC-specific visual.
    - It should be paired conceptually with KASH5/dynein, even though KASH5 is not the main image source here.
    - Keep as optional mechanism material for the final talk, not as required main-flow biology.
  ],
  [Chen et al. 2021, Nature Communications Fig. 5. CC BY 4.0. Springer Nature image.],
  fig_width: 64%,
  img_height: 5.65in,
)

#pagebreak()

#text_slide(
  "11",
  "Sources / Citations",
  "Sources and citations",
  [
    #set par(leading: 0.58em)
    #text(size: 9.2pt)[
      *Figures shown in this review deck*

      Wang X, Pepling ME. 2021. Regulation of Meiotic Prophase One in Mammalian Oocytes. Front. Cell Dev. Biol. 9:667306. DOI 10.3389/fcell.2021.667306. CC BY.

      Cheng EY et al. 2009. Meiotic Recombination in Human Oocytes. PLOS Genet. 5(9):e1000661. DOI 10.1371/journal.pgen.1000661. CC BY.

      Gruhn JR et al. 2013. Cytological Studies of Human Meiosis: Sex-Specific Differences in Recombination Originate at, or Prior to, Establishment of Double-Strand Breaks. PLOS ONE 8(12):e85075. DOI 10.1371/journal.pone.0085075. CC BY 4.0.

      Blokhina YP et al. 2019. The telomere bouquet is a hub where meiotic double-strand breaks, synapsis, and stable homolog juxtaposition are coordinated in the zebrafish, Danio rerio. PLOS Genet. 15(1):e1007730. DOI 10.1371/journal.pgen.1007730. CC BY 4.0.

      Bisig CG et al. 2012. Synaptonemal Complex Components Persist at Centromeres and Are Required for Homologous Centromere Pairing in Mouse Spermatocytes. PLOS Genet. 8(6):e1002701. DOI 10.1371/journal.pgen.1002701. CC BY.

      Dunce JM et al. 2018. Structural basis of meiotic telomere attachment to the nuclear envelope by MAJIN-TERB2-TERB1. Nat. Commun. 9:5355. DOI 10.1038/s41467-018-07794-7. CC BY 4.0.

      Chen Y et al. 2021. The SUN1-SPDYA interaction plays an essential role in meiosis prophase I. Nat. Commun. 12:3176. DOI 10.1038/s41467-021-23550-w. CC BY 4.0.

      #v(0.13in)
      Liebe B et al. 2004. Telomere attachment, meiotic chromosome condensation, pairing, and bouquet stage duration are modified in spermatocytes lacking axial elements. Mol. Biol. Cell 15(2):827-837. DOI 10.1091/mbc.e03-07-0524. PMC/free-to-read source; exact rights note in manifest.

      #v(0.13in)
      *Additional source notes*

      Mytlis et al. bioRxiv/Science, Lenzi et al. AJHG, Shibuya et al. Cell, and Morimoto et al. JCB were considered from the brief but not downloaded for this compact deck. Full direct URLs, figure/panel notes, rights labels, attribution strings, and skipped-source notes are in image_manifest.tsv.
    ]
  ],
)
