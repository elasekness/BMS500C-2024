# BMS500C-2024
Tutorials for 2024 BMS500C Bioinformatics

The overall objectives of the bioinformatics modules are to:

1. become familiar with the command-line (Linux) environment, including navigating the directory structure, making, viewing and editing files, and learning other basic BASH commands.

2. perform more sophisticated command-line operations such as using regular expressions (regex), piping, and for-loops to efficiently manipulate and reformat data.

3. employ commonly-used software packages to process fastq files (including adapter removal and quality filtering), perform both reference-based and de-novo assemblies, and generate summary statistics for each.

4. understand the basics of BLAST and homology searches by querying your assembled contigs against NCBI's core nucleotide database

5. generate genealogies to assess the taxonomic identity of your contigs

We will accomplish these objectives by working with _Cryptosporidium_ next-generation sequencing data.

_Cryptosporidium_ is a genus of protistan parasites belonging to the phylumm Apicomplexa, which includes _Plasmodium_, the etiological agent of malaria. It is considered a zoonotic disease (capable of transfer
between humans and other animals) with a broad range of hosts.  _Cryptosporidium_ resides in the gastrointestinal tracts of humans where it can cause severe diarrhoea, particularly in children and immunocompromised 
individuals.  Along with _Giardia_, it is a leading cause of diarrhoeal disease worldwide.

There are almost 50 validated species of _Cryptosporidium_ and numerous genotypes.  Some species, such as _Cryptosporidium parvum_, are capable of infecting numerous vertebrates
while others appear more species-specific.

Delimitation of _Cryptosporidium_ species relies on both morphological and molecular evidence.
Genotyping typically occurs by sequencing one or more of the following genes:
1. Actin
2. 18S small subunit ribosomal gene
3. HSP 70 (heat shock protein)

During this module, we will perform _Cryptosporidium_ species identification for various _Cryptosporidium_ positive samples.
Pooled amplicons of actin, hsp70, and ssu from these samples were sequenced on an Illumina NextSeq machine to generate ~1/2 million paired-end reads.
We will assemble these reads for each amplicon and determine the phylogenetic affiliation of the isolate.
