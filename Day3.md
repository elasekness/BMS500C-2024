## Tutorial 3

**Objective:** To determine the _Cryptosporidum_ species present in your sample by performing BLAST homology searches and generating species-level phylogenies of _Cryptosporidium_.

BLAST is a tool to find regions of similarity between sequences in either protein or nucleotide format.  We can infer homology between sequences
based on the E-value, which estimates how many random alignments in the database score as good or better than the match (alignment) in question.
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is available as a webtool but is also installed on our VMs.
We will use both web-based BLAST and command-line BLAST functions to infer the _Cryptosporidium_ species present in your sample and compare the results.

Phylogenies depict relationships among taxa.  They are inferences about the evolutionary history of a group of organisms.
Phylogenies can be generated from morhpological or molecular data or a combination of both.  Methods for reconstructing
relationships include distance-based methods, parsimony, maximum likelihood (ML), and Bayesian inference. 
Tree reconstruction with molecular data requires a multi-sequence alignment of homologous genes (for reconstructing species relationships we want to align orthologs).
With ML and Bayesian methods, a model of protein or nucleotide substitution, which determine how character states change over time, must also be specified.
In ML methods, we determine the probablity of the data (the alignment) given the tree (including the tree-topology and branch lengths) and subsituttion model. 
The best tree is considered the tree that maximizes the likelihood of the data.  We will use [IQ-Tree](http://www.iqtree.org/) to infer phylogenies
with our actin, Hsp70, and ssu datasets.  The results will help us to phylogenetically place our unkown _Cryptosporidium_ species.

<br>

First let's use the command-line BLAST functions to determine the best matches between our assembled genes and those of the reference database.
We will create a blast database using our reference.fasta file, which contains the _Cryptosporidium_ sequences we downloaded from NCBI.

## Make a blast database

	makeblastdb -in reference.fasta -dbtype nucl

## BLAST your contigs against the reference database

	blastn -query final.contigs.fa -db reference.fasta -max_target_seqs 5 -outfmt '6 qaccver saccver pident qcovs length qlen mismatch gapopen qstart qend sstart send' -out contigs.br

> `blastn` searches are for nucleotide queries against a nucleotide database
> `-max_target_seqs` specifies the maximum number of hits to return
> `-outfmt 6` returns the results in tab-delimited format. The information proceeding the 6 specifies the the fields to return.
> We are outputting the query accession, the subject accession, the percent identity between query and subject, the percent coverage of the query by the hit, the alignment legnth, the query length,
the number of mismatches between query and subject, the number of gaps, and the start and end of the alignment in the query and subject, respectively.
* Based on your BLAST results:
*  - are all genes represented in your contigs file?
*  - is the tophit from the same species if you have multiple genes present?
*  - can you confidently assign a species to your sample given the PID, qcov, qlen, etc?


## BLAST your contigs using the BLAST webpage and the nucleotide core (core_nt) database

Navigate to the [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) webpage.
Notice that nucleotide core is the default database.
Copy and paste your contigs into the the Query sequence section and execute the blast search.  
You can print your contigs to STDOUT with the `cat` command.

	cat final.contigs.fa

* Based on your BLAST results:
*  - are all of your contigs from _Cryptosporidium_?
*  - do these results agree with your command-line blast results?
 
 
