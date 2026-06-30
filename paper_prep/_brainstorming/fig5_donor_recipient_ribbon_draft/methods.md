# Fig. 5 Methods Note

The Fig. 5 donor-recipient ribbons use the PAN027 paternal haplotype as query
against a joint two-haplotype PAN011 father target. They are derived from the
existing SweepGA/FastGA frequency-32, many:many, no-scaffold raw PAF, filtered
with sweepga 0.1.1 to 1:1, 4:4, and 10:10 mapping classes. The manuscript figure
uses the 10:10 class. Calls are summarized in 2 kb query-space windows after
excluding centromeric windows, windows without alignment support, and windows
with excessive inter-chromosomal PAF depth (maximum depth 100).

For each retained window, `impg similarity` 0.4.1 was run with POA, no merge,
many:many mapping, and scaffold jump 0. A window was classified by the best
inter-chromosomal PAN011 father match when that match beat the same-chromosome
or homologous father match. Adjacent 2 kb windows were merged only for display
when donor chromosome, donor haplotype, donor sequence, and nearby recipient and
donor coordinates were consistent. PHR and arm-community labels were assigned
post hoc and were not used to call the windows.
