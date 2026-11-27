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

- **Righe 344-353 (Results) e 1130-1133 (Methods): togliere l'argomento del "limite inferiore".** Il testo dice che le read multiple, messe a caso, attenuano i contatti in modo simmetrico, spingono rho e B/W verso il nullo, e che un artefatto metterebbe le PHR sotto i fianchi. Le simulazioni nulle dicono il contrario: senza nessun contatto vero tra estremità, i soli errori di mappatura danno un Mantel di 0,76-0,82 in Hi-C e 0,64 in Pore-C, perché una read finisce sulle PHR più simili (sezioni 6-7). Anche l'ordinamento PHR 0,657 sopra fianchi 0,520 non esclude più l'artefatto.
- **Dopo la riga 343 (Results): nuovo paragrafo.** Simulazioni nulle a livello di read con profili di errore presi dalle read vere; calibrazione sulla profondità (gli errori spiegano fino al 56% dei contatti tra PHR in Hi-C HG002, fino all'86% in CHM13, dal 9,0% all'82,0% in Pore-C); Mantel del paper sulle matrici senza i contatti falsi attesi (0,40-0,66 in Hi-C, 0,47-0,48 in Pore-C, p < 0,0001) (sezioni 6-8).
- **Results: una frase sul controllo con duplicazioni segmentali appaiate.** In Pore-C e CiFi le PHR stanno sopra i controlli (0,420 contro 0,202; 0,191 contro 0,071), in Hi-C no (0,488 contro 0,493), per gli errori di mappatura (sezione 4).
- **Results o Extended Data: Hi-C con sole read uniche.** Mantel 0,487 in HG002 (MAPQ >= 2) e 0,503 in CHM13 con il bilanciamento che non scarta i bin delle PHR (`cooler balance --mad-max 0`); con il bilanciamento di default CHM13 dà -0,023 perché scarta 19 delle 38 PHR (sezione 3).
- **Methods: nuove sottosezioni.** Simulazioni nulle (Hi-C con HiC-Pro, Pore-C e CiFi con la catena del paper; profili di errore campionati lungo tutti i file), calibrazione, sottrazione degli artefatti alla risoluzione di ogni statistica; controllo con duplicazioni (SEDEF, appaiamento su identità, lunghezza e copie, nuovo grafo pggb, rho per coppia con bootstrap); bilanciamento delle matrici con sole read uniche.
- **Extended Data: tabelle** delle sezioni 4, 7 e 8 del file di prove.
- **Discussion: una frase.** Il contatto Hi-C tra PHR è gonfiato dagli errori di mappatura; la conclusione poggia sulle statistiche corrette e sui saggi con sole read uniche.
- **Fig. 4A: dati invariati.** Nella legenda o nel testo si può dire quale quota dei contatti Pore-C tra PHR può venire dagli errori (9,0-82,0%).

