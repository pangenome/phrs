# Modifiche proposte a submission/paper.tex dopo l'esclusione delle 28 estremità

Stato al 2026-11-25 21:30 MST. paper.tex non è stato modificato. Numeri delle righe riferiti a paper.tex di oggi.

Fonte dei numeri: chiamata nuova `w500clean.f773342` (stesso programma del paper con la correzione del conteggio dei contig con "chr" nel nome; 28 estremità tolte come bersagli e come partner), `/moosefs/guarracino/HPRCv2/PHR_III/flank5mb/check/out/`. Le 28 estremità sono in `remove28.txt` (22 giunzioni sbagliate e 3 sbagliate per meno di 10 kb da `scripts/phr-detection/flank_joins_read_check.tsv`, più HG03521#2 3q, HG02922#1 Xq e NA18982#1 18q).

## Numeri già ricalcolati

| riga | ora | dopo | come |
| --- | --- | --- | --- |
| 163 | 15,668 PHRs | 15,646 | w500clean.calls.tsv |
| 173 | median 105 kb, mean 144 kb | uguale (105 kb; media 143.8 → 143.7 kb) | w500clean.calls.tsv |
| 841, 845, 856 | 18,826 flanks | 18,799 | 18,827 meno le 28 |
| 856 | 15,668 PHRs (83.2%) | 15,646 (83.2%) | |
| 857-858 | stop after 1, 2, 3 failed windows: 15,624, 15,664, 15,668 | 15,605, 15,643, 15,646 | `scripts/phr-detection/phr_sensitivity.py`: sui dati del paper riproduce esattamente 15,624 / 15,664 / 15,668, 15,668 e 40 / 38 / 37 bracci con 14,369 / 13,422 / 10,030 |
| 858-859 | 25-400 kb gives 15,668 | 15,646 | come sopra |
| 859-860 | 3, 6, 15 chromosomes: 40, 38, 37 arms (14,369, 13,422, 10,030 PHRs) | 40, 38, 37 arms (14,348, 13,358, 10,015 PHRs) | come sopra |
| 859 | 41 arms | 41 | |

Le 5 PHR che cambiano lunghezza per la correzione del programma (+5 kb ciascuna): HG00673#2, HG04199#1, HG04228#2, NA20799#1 chr9q e NA20827#2 chr14p. Le altre 15,641 PHR sono identiche al paper (nome e sequenza).

## Testo da riscrivere

- **Righe 837-841 (esclusione delle estremità).** Oggi descrive solo NA18982 18q, esclusa "prima dell'allineamento" (in realtà la pipeline l'aveva tolta dopo la chiamata). Proposta: descrivere la ricerca di buchi di N e di array telomerici interni in tutte le 18,827 estremità (112 segnalate), il controllo con le read ONT e HiFi di HPRC (read che attraversano la giunzione; read che finiscono nel telomero alla giunzione; lato verso il telomero senza read), e l'esclusione di 28 estremità come bersagli e come partner prima della chiamata delle PHR.
- **Riga 845.** L'allineamento tutti contro tutti è stato fatto su 18,827 estremità; le righe che coinvolgono le 28 escluse sono state tolte prima della chiamata. Il testo dovrebbe dirlo così.
- **Righe 951-953 (Leiden dei bracci).** "we selected the maximum-silhouette partition, 15 communities at resolution 1.16": la divisione usata nel paper ha silhouette 0.332, mentre la divisione più frequente alla stessa risoluzione ha 0.347 (studio con 1,000 semi sulla matrice del paper). Da decidere con i risultati sul grafo nuovo.

## In attesa del grafo pggb nuovo (job 2881196)

- Righe 182, 942, 944, 1033: "15,668" diventa 15,646; nodi e archi del grafo (477,726 e 786,318) da sostituire con quelli del grafo nuovo.
- Communities dei bracci (15 Leiden, 14 UPGMA) e delle sequenze: da ricalcolare sul grafo nuovo; lo studio con 1,000 semi darà tutte le divisioni di Leiden con frequenza e silhouette.

## Da chiarire

- **Riga 176, "3.51 Mb".** Riprodotto: è la somma delle lunghezze delle PHR di CHM13 fuori dagli acrocentrici p e da X/Y (3.510 Mb; la somma sui bracci della lunghezza media per braccio, contando 0 le estremità senza PHR, dà 3.511 Mb). Con o senza le 28 estremità resta 3.51 Mb (CHM13 non è tra le 28). Il BED data/chm13.phrs.no_acro.bed è un'altra cosa (4.075 Mb).
- **Righe 397, 624, 935 ("18,827 windows").** Le analisi di divergenza e dei telomeri includevano le 28 estremità. Vanno rifatte senza, oppure il testo deve dire che le includono.

## Altre correzioni trovate durante la revisione

- **Righe 1086-1095 (Methods, read multiple; erano 1062-1069 prima degli ultimi commit di paper.tex).** Il testo dice che per Hi-C, Pore-C e CiFi teniamo le read che mappano in più punti e le mettiamo a caso su una delle estremità. Vale solo per Hi-C (HiC-Pro con `MIN_MAPQ = 0` e `RM_MULTI = 0`). Per Pore-C e CiFi, `pairtools parse2` 1.1.3 ha `--min-mapq` 1 per default (non passato nei comandi), mette `!` come cromosoma ai lati con MAPQ 0, e `cooler cload` scarta queste coppie. Conteggi: Pore-C UU 303,703,855 (141,267,621 + 162,436,234) contro una somma di 286,293,292 contatti nel `.mcool` a 5 kb; CiFi UU 39,983,047 contro 35,560,326; le coppie con lati multipli sono 1,359,086,401 MM e circa 360 milioni MU e UM (Pore-C, flow cell 1). Proposta: dire che Pore-C e CiFi usano solo coppie di monomeri con MAPQ >= 1 su entrambi i lati, e che solo Hi-C tiene le read multiple. Questo risponde anche a R1-M4 (contatti ancorati a un monomero mappabile in modo unico) per Pore-C e CiFi, e la Fig. 4A (HG002 Pore-C) usa già solo contatti unici.

## Somiglianza e contatto 3D: modifiche dalle prove per R1-M4, R1-M5, R3-M3 (2026-11-27)

Numeri e prove in `paper_prep/manuscript_revision/similarity_contact_evidence.md` (sezioni indicate). Righe del paper.tex attuale.

- **Righe 344-353 (Results) e 1130-1133 (Methods): togliere l'argomento del "limite inferiore".** Il testo dice che le read multiple, messe a caso, attenuano i contatti in modo simmetrico, spingono rho e B/W verso il nullo, e che un artefatto metterebbe le PHR sotto i fianchi. Le simulazioni nulle dicono il contrario: senza nessun contatto vero tra estremità, i soli errori di mappatura danno un Mantel di 0,755-0,815 in Hi-C e 0,643 in Pore-C, perché una read finisce sulle PHR più simili (sezioni 6-7). Anche l'ordinamento PHR 0,657 sopra fianchi 0,520 non esclude più l'artefatto.
- **Dopo la riga 343 (Results): nuovo paragrafo.** Simulazioni nulle a livello di read con profili di errore presi dalle read vere; calibrazione sulla profondità (gli errori spiegano fino al 56% dei contatti tra PHR in Hi-C HG002, fino all'86% in CHM13, dal 9,0% all'82,0% in Pore-C); Mantel del paper sulle matrici senza i contatti falsi attesi (0,40-0,66 in Hi-C, 0,47-0,48 in Pore-C, p < 0,0001) (sezioni 6-8).
- **Results: una frase sul controllo con duplicazioni segmentali appaiate.** In Pore-C e CiFi le PHR stanno sopra i controlli (0,420 contro 0,202; 0,191 contro 0,071), in Hi-C no (0,488 contro 0,492), per gli errori di mappatura (sezione 4).
- **Results o Extended Data: Hi-C con sole read uniche.** Mantel 0,487 in HG002 (MAPQ >= 2) e 0,503 in CHM13 con il bilanciamento che non scarta i bin delle PHR (`cooler balance --mad-max 0`); con il bilanciamento di default CHM13 dà -0,023 perché scarta 19 delle 38 PHR (sezione 3).
- **Methods: nuove sottosezioni.** Simulazioni nulle (Hi-C con HiC-Pro, Pore-C e CiFi con la catena del paper; profili di errore campionati lungo tutti i file), calibrazione, sottrazione degli artefatti alla risoluzione di ogni statistica; controllo con duplicazioni (SEDEF, appaiamento su identità, lunghezza e copie, nuovo grafo pggb, rho per coppia con bootstrap); bilanciamento delle matrici con sole read uniche.
- **Extended Data: tabelle** delle sezioni 4, 7 e 8 del file di prove.
- **Discussion: una frase.** Il contatto Hi-C tra PHR è gonfiato dagli errori di mappatura; la conclusione poggia sulle statistiche corrette e sui saggi con sole read uniche.
- **Fig. 4A: dati invariati.** Nella legenda o nel testo si può dire quale quota dei contatti Pore-C tra PHR può venire dagli errori (9,0-82,0%).


## Somiglianza e contatto 3D: testo proposto (bozza del 2026-11-27, 21:18 MST)

Testo in inglese pronto per paper.tex. Ogni numero viene da `paper_prep/manuscript_revision/similarity_contact_numbers.md` (stampato da `scripts/hic/similarity_contact_numbers.py`). I valori di Pore-C con le read multiple (job 2890449) arriveranno sabato 28/11 e andranno aggiunti alla frase sulla sottrazione. Righe del paper.tex attuale.

### Results: sostituisce le righe 344-353 (da "Both statistics are biased" a "which is not observed.")

```latex
Misplaced reads could create this pattern on their own, because a read from one PHR is most likely to be placed on the PHRs most similar to it.
To measure this, we simulated reads whose two ends come from the same chromosome end, so that every contact between two ends is a mapping error, and we processed them with the same pipelines and with error profiles sampled from the real reads (Methods).
In HG002 Hi-C, the rate at which reads from one PHR were placed on another followed the sequence similarity of the two PHRs (Spearman $\rho=0.774$ over 2,830 pairs).
The simulated matrices alone gave an arm-level Mantel $\rho$ of 0.755 to 0.815 in Hi-C and 0.643 in Pore-C at 50\,kb, as high as the observed values, so a rank statistic cannot separate placement errors from proximity and the flank comparison above does not exclude them.
Scaled to the depth of each dataset, placement errors explain up to 55.7\% of the contacts between PHRs in HG002 Hi-C, up to 85.9\% in CHM13 Hi-C, 9.0\% to 82.0\% in Pore-C and 2.8\% to 7.9\% in CiFi.
After we subtracted the upper bound of the expected errors from each matrix, the arm-level Mantel $\rho$ stayed significant in every assay ($p<10^{-4}$, 10,000 permutations): 0.404 at 50\,kb and 0.574 at 20\,kb in HG002 Hi-C, 0.468 and 0.661 in CHM13 Hi-C, 0.479 and 0.473 in Pore-C, and 0.427 and 0.502 in CiFi with multi-mapping reads kept (Extended Data Table~\ref{tab:ed_simcorr}).
Hi-C restricted to uniquely mapped reads (MAPQ $\geq2$ on both reads) gave $\rho=0.487$ in HG002 and 0.503 in CHM13 at 50\,kb ($p<10^{-4}$).
```

### Results: subito dopo, controllo con duplicazioni segmentali appaiate

```latex
To ask whether subtelomeric PHRs differ from other interchromosomal duplications, we compared each PHR pair with up to three segmental duplications at least 1\,Mb from any chromosome end, matched on identity, length and copy number (Methods).
In HG002 Pore-C, similarity tracked contact more strongly for PHR pairs than for the matched duplications (Spearman $\rho=0.420$ against 0.202, difference 0.218, 95\% bootstrap interval 0.119 to 0.311), and likewise in CiFi (0.191 against 0.071, difference 0.120, 0.032 to 0.213).
In HG002 Hi-C with all reads, PHRs and matched duplications did not differ (0.488 against 0.492, difference $-0.004$, $-0.119$ to 0.118), as expected if placement errors inflate contact between any two similar loci in short reads.
With uniquely mapped Hi-C reads, PHRs lay above the matched duplications in HG002 (0.362 against 0.232 at 20\,kb, difference 0.130, 0.003 to 0.247).
```

### Methods: sostituisce le righe 1086-1095 (da "We reprocessed Hi-C" a "not reads confidently assigned to a single end.")

```latex
We reprocessed Hi-C, Pore-C and CiFi from raw reads against the corresponding T2T reference (CHM13v2.0 or the matched HPRC v2 haplotype).
We processed Pore-C (NlaIII) and CiFi (DpnII) with pore-c-py 2.1.5, minimap2 2.30 and pairtools 1.1.3, and Hi-C with HiC-Pro 3.0.0 and Bowtie2 2.5.3, binning contacts with cooler 0.10.4 and samtools 1.23 under Python 3.12.12 (numpy 1.26.4, pandas 3.0.0, scipy 1.17.0).
For Hi-C we kept multi-mapping reads (HiC-Pro \texttt{MIN\_MAPQ} 0, \texttt{RM\_MULTI} 0, Bowtie2 local mode), so a read that maps equally well to several paralogous ends is placed at one of them at random, and contact inside PHRs is an aggregate over these near-identical ends.
For Pore-C and CiFi, \texttt{pairtools parse2} kept only uniquely mapped monomers (its default \texttt{--min-mapq 1}).
```

### Methods: togliere le righe 1130-1133

Da "Because uniform placement of a multi-mapping read" a "is a lower bound." senza sostituzione.

### Methods: nuova sottosezione dopo "Hi-C, Pore-C and CiFi pipeline"

```latex
\subsection*{Read-level simulation of placement errors}

For each chromosome end we defined an end domain (the PHR window plus its 100\,kb centromere-ward flank) and simulated reads with no contact between ends.
For Hi-C, the first read of each pair started uniformly in the end domain and the second lay at a distance $s$ with $P(s)\propto1/s$ on [1\,kb, 500\,kb], towards the telomere or the centromere at random, at 40 read pairs per kb of end domain (863,400 pairs over 77 HG002 ends, 425,600 over 38 CHM13 ends).
For Pore-C and CiFi, each read held $K\geq2$ monomers ($K-2$ Poisson with mean 2, at most 20), each the NlaIII or DpnII restriction fragment at a position drawn in the same way, joined by the restriction site, at 20 reads per kb of end domain.
We ran two error models.
In the first, Hi-C reads (250\,bp in HG002, 150\,bp in CHM13, as the real reads) carried 0.2\% substitutions, Pore-C reads 0.5\% substitutions, 0.25\% insertions and 0.25\% deletions, and CiFi reads one tenth of these rates.
In the second, each simulated Hi-C read took the length, base qualities and mismatch positions of a random real read of the same mate, with random bases over its 5$'$ soft clip and the sequence of the partner locus over its 3$'$ soft clip, and each Pore-C read took the substitution, insertion and deletion rates of a random real read and each monomer the soft-clipped fractions of a random real monomer.
We sampled the real reads at regular intervals along every input file (every 3,000th Hi-C read in HG002, every 1,500th in CHM13, every 2,000th Pore-C read), because the first reads of a file carry more errors than the rest.
Mapped with the same options, the second model matched the real reads in mismatches per aligned base (HG002 Hi-C read 1: 1.17\% simulated, 1.15\% real) and in the fraction of pairs passing MAPQ $\geq1$ and $\geq2$ (58.5\% and 7.2\% simulated, 54.7\% and 7.48\% real).
We processed the simulated reads with the pipelines above and counted the contacts placed in the PHR of another end.

\subsection*{Depth calibration and subtraction of expected errors}

For each end we computed a scale factor, the ratio of real to simulated raw cis contacts inside its PHR window at 5\,kb.
For a pair of ends, the expected false contacts are the simulated contacts multiplied by a scale factor, and we took the smaller and the larger of the two ends' factors as lower and upper bounds.
We subtracted the upper bound from each raw matrix pixel by pixel at the resolution of each statistic (50 and 20\,kb), set negative pixels to zero, balanced the matrix as for the original one and recomputed the arm-level Mantel test.
As controls, the subtraction applied to a simulated matrix in place of a real one left no contact between ends, and without subtraction the procedure reproduced the original tables exactly.

\subsection*{Hi-C with uniquely mapped reads}

From the HiC-Pro valid pairs we kept the pairs with MAPQ $\geq1$ or $\geq2$ on both reads and built the matrices as above.
In this setup Bowtie2 gives multi-mapping reads MAPQ 0 or 1, so MAPQ $\geq2$ keeps uniquely mapped reads.
Because the default \texttt{cooler balance} filter masks low-coverage bins, which in CHM13 at MAPQ $\geq2$ removed the bins of 19 of 38 PHR windows, we balanced the CHM13 matrices with \texttt{--mad-max 0}, which kept all 38.

\subsection*{Matched segmental duplication control}

We took the SEDEF segmental duplications of CHM13 and HG002 whose two copies lie on different chromosomes, merged the copies on each chromosome into loci, and described each pair of loci by identity (weighted by aligned length), total aligned length and copy number (the mean number of other chromosomes the two loci share duplications with).
A PHR pair has both loci overlapping a PHR, and a control candidate has both loci at least 1\,Mb from either chromosome end.
For each PHR pair we chose up to three candidates nearest in identity, log length and log copy number within 0.5 standard deviations, each used once (HG002: 380 of 734 PHR pairs matched, 1,069 controls).
We built one pggb graph per genome from all PHR and control loci with the parameters used for the PHRs, computed Jaccard similarity with \texttt{odgi similarity}, computed contact per pair as above, and compared Spearman $\rho$ between similarity and contact for PHR and control pairs, with 95\% intervals from 1,000 bootstrap resamples of loci.
```

### Extended Data: una tabella (`tab:ed_simcorr`)

Colonne: saggio, Mantel del paper, Mantel nulla (simulazione realistica), contatti tra PHR spiegati dagli errori (min-max), Mantel dopo la sottrazione a 50 e 20 kb. Righe: Hi-C HG002, Hi-C CHM13, Pore-C, CiFi con read multiple, Hi-C MAPQ >= 2 HG002 e CHM13. Una seconda tabella con il controllo SD (righe della sezione 4 di `similarity_contact_numbers.md`).

### Discussion: una frase

```latex
In short-read Hi-C, placement errors explain up to 55.7\% of the contact between PHRs in HG002 and up to 85.9\% in CHM13, and they follow sequence similarity. The association between similarity and contact holds after we remove the expected errors and in long-read Pore-C and CiFi, where PHRs exceed matched segmental duplications.
```
