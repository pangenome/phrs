// Slide-ready CHM13/hs1 UCSC browser examples for review zoom v8.
// Generated for task review-zoom-v8-chm13-ucsc-example-selection.

#set page(width: 13.333in, height: 7.5in, margin: 0.36in)
#set text(font: "DejaVu Sans", fill: rgb("#111827"))
#set par(justify: false, leading: 0.62em)

#let muted = rgb("#4b5563")
#let rule = rgb("#d1d5db")
#let accent = rgb("#1f4f82")

#let example(num, title, reason, panel, source) = [
  #block(width: 100%)[
    #text(size: 9pt, fill: accent, weight: "bold")[CHM13 UCSC example #num]
    #v(0.04in)
    #text(size: 23pt, weight: "bold")[#title]
    #v(0.08in)
    #text(size: 12.5pt, fill: muted)[#reason]
  ]
  #v(0.12in)
  #line(length: 100%, stroke: 0.55pt + rule)
  #v(0.14in)
  #box(width: 100%, height: 3.45in)[
    #align(center + horizon)[
      #image(panel, width: 100%)
    ]
  ]
  #v(0.10in)
  #line(length: 100%, stroke: 0.55pt + rule)
  #v(0.06in)
  #text(size: 6.8pt, fill: muted)[#source]
]

#example(
  "01",
  [chr4q - C1 D4Z4/DUX4L],
  [Best single real-UCSC view of the D4Z4/DUX4L system: repeated DUX4L/MIR8078 labels sit inside the PHR interval.],
  "source_panels/01_chr4q_c1_d4z4_dux4l_ucsc_panel.png",
  [UCSC hs1/CHM13; PHR BED chr4:193392741-193572740; browser chr4:193,304,946-193,574,945; manifest row 28; panel 01_chr4q_c1_d4z4_dux4l_ucsc_panel.png],
)

#pagebreak()

#example(
  "02",
  [chr3q - C3 OR4F/f7501 cluster],
  [C3 is a large interpretable community; this panel shows OR4F5 plus the broader duplicated gene context in one clear browser window.],
  "source_panels/02_chr3q_c3_or4f_f7501_ucsc_panel.png",
  [UCSC hs1/CHM13; PHR BED chr3:200846363-201101362; browser chr3:200,723,449-201,105,948; manifest row 26; panel 02_chr3q_c3_or4f_f7501_ucsc_panel.png],
)

#pagebreak()

#example(
  "03",
  [chr15q - C8 OR4F endpoint],
  [Clean OR4F-rich endpoint with multiple OR4F labels and WASH/DDX11L context; useful contrast to the multi-arm C3 example.],
  "source_panels/03_chr15q_c8_or4f_endpoint_ucsc_panel.png",
  [UCSC hs1/CHM13; PHR BED chr15:99625359-99750358; browser chr15:99,565,696-99,753,195; manifest row 11; panel 03_chr15q_c8_or4f_endpoint_ucsc_panel.png],
)

#pagebreak()

#example(
  "04",
  [chr18p - C2 repeat/gene context],
  [C2/chr18p is report-backed as TAR1-rich; the real UCSC view visibly shows repeated duplicon blocks, IL9RP4, and terminal PHR structure.],
  "source_panels/04_chr18p_c2_tar1_repeat_gene_ucsc_panel.png",
  [UCSC hs1/CHM13; PHR BED chr18:2017-267017; browser chr18:1-397,502; manifest row 16; panel 04_chr18p_c2_tar1_repeat_gene_ucsc_panel.png],
)

#pagebreak()

#example(
  "05",
  [chr21p - C7 acrocentric p-arm],
  [C7 acrocentric p-arm example: long PHR, rDNA-adjacent duplicated structure, and MTCO/SNX18 pseudogene labels at slide scale.],
  "source_panels/05_chr21p_c7_acrocentric_p_ucsc_panel.png",
  [UCSC hs1/CHM13; PHR BED chr21:2505-487505; browser chr21:1-727,502; manifest row 22; panel 05_chr21p_c7_acrocentric_p_ucsc_panel.png],
)

#pagebreak()

#example(
  "06",
  [chrXp - C15 PAR1/SHOX],
  [Shows PAR1 coding genes including SHOX while making the key point that the PHR begins internal to the telomeric end, not at the telomere itself.],
  "source_panels/06_chrXp_c15_par1_shox_ucsc_panel.png",
  [UCSC hs1/CHM13; PHR BED chrX:1579-501579; browser chrX:1-750,002; manifest row 36; panel 06_chrXp_c15_par1_shox_ucsc_panel.png],
)
