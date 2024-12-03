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

