# aav-gene-therapy
This repo contains the code associated with the paper "Contaminating manufacturing plasmids and disrupted vector genomes present in liver tissue
following adeno-associated virus gene therapy".

## Short read metagenomic sequencing
Data was first analysed with metaMix using the pipeline described at https://github.com/smorfopoulou/clinical_metagenomics.
The raw (pSMN plasmid) or human-filtered (all other reference sequences) reads were aligned to the reference genomes using gt_alignments_illumina.sh and plots were produced using alignment_depth.R.
Taxonomic classification was also run with the nf-core pipeline taxprofiler using Kraken2 and Bracken, with the script taxprofiler_kraken2.sh.

## Long read metagenomic sequencing
Following basecalling with Minknow, data was preprocessed using preprocess_nanopore.sh.
The trimmed reads were aligned to the plasmid sequences using gt_alignments_nanopore.sh and filtered using filter_sam.R.
Alignment dot plots were produced from teh outputs using gt_dot_plots.sh.
Taxonomic classification was also run with the nf-core pipeline taxprofiler using Kraken2 and Bracken, with the script taxprofiler_kraken2.sh.

## Version of main tools used in analysis
bowtie2 v2.4.1

samtools v1.14

minimap2 v2.17

redotable v1.1

porechop 0.2.4

nextflow 23.04.3

taxprofiler v1.1.0 

bracken v2.8

R v4.2.2

tidyverse v2.0.0

GenomicAlignments v1.32.1

metaMix v0.4

## Reference genomes
Human genome GRCh38 was used for all analyses.

Viral reference genome accession numbers can be found in reference_accessions.csv.

The manufacturing plasmid sequences are not available on RefSeq/GenBank to the best of our knowledge.

pSMN:  Kaspar, B. K., Burghes, A. & Porensky, P. Intrathecal delivery of recombinant Adeno-associated virus 9. (2022). AMR gene from patent sequence replaced by KanR (see methods section).

pAAV2/9: Gao, G., Wilson, J. & Alvira, M. Adeno-associated virus (aav) serotype 9 sequences, vectors containing same, and uses therefor. (2005).

pHelper: Gray, J. Molecule Information, pHGTI-Adeno1, Harvard Gene Therapy Initiative. (2004).







