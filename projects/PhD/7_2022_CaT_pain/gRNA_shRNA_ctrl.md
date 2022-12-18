sgRNA design & quality control
==============================
*Borys Olifirov, 19.09.2022*

Target - Rattus norvegicus **T-type Cav3.2 subunit alpha1 H**, *Cacna1h* gene (transcript ID NM153814.2, Gene ID: 114862) 

## sgRNA
Delivery of dual sgRNA cassettes and Nme2Cas9 in a single AAV vector, mouse targeting (addgene plasmid **#159537**).

Existing gRNA's spacers were designed by rat *Cacna1h* mRNA sequence.

Instead of a check of single spacer sequences I prefer to design a new set of spacers. I used a highly recommended online tool for sgRNA design CRISPick (https://portals.broadinstitute.org/gppx/crispick/public).

The top 10 results are presented in **sgRNA_top_10.xlsx** file with all supplementary information including potential cut efficacy and off-target cutting sites (columns legend you may find by link https://portals.broadinstitute.org/gppx/crispick/public/how-to-use).

Was selected sequence 1 (target exon **17**), 2 (target exon **6**) and 4 (target exon **8**). To finalize sgRNA design I added BspQI or BsmBI cohesive ends to double-strand sgRNA sequences. 

- **gRNA I** (BspQI)
    5'-AAC / CCA-5' << U6
- **gRNA II** (BsmBI)
    U6 >> 5'-CACC / CAAC-5'

Deletion/inversion region size with different sgRNAs pairs:

sgRNA pair|Flanked region size|Size of CDS with deletion, bp
----------|-------------------|-----------------------------
6+8       |1058 bp            |6023 bp
6+17      |2984 bp            |4096 bp
8+17      |1926 bp            |5155 bp

As a negative control, it is recommended to use the sgRNA for a gene absent in the model organism (GFP, luciferase, etc.). All sgRNA sequences with cohesive ends including negative control with GFP saved in **sgRNA.fasta**.