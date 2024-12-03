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
With ML and Bayesian methods, a model of protein or nucleotide substitution, which determines how character states change over time, must also be specified.
In ML methods, we determine the probablity of the data (the alignment) given the tree (including the tree-topology and branch lengths) and subsituttion model. 
The best tree is considered the tree that maximizes the likelihood of the data.  We will use [IQ-Tree](http://www.iqtree.org/) to infer phylogenies
with our actin, Hsp70, and ssu datasets.  The results will help us to phylogenetically place our unkown _Cryptosporidium_ species.

You will also need to install multi-sequence alignment ([SeaView](https://doua.prabi.fr/software/seaview) and tree-viewing [FigTree](https://github.com/rambaut/figtree/releases) software on your personal computer.

<br>

First let's use the command-line BLAST functions to determine the best matches between our assembled genes and those of the reference database.
We will create a blast database using our reference.fasta file, which contains the _Cryptosporidium_ sequences we downloaded from NCBI.

## Make a blast database

	makeblastdb -in reference.fasta -dbtype nucl

## BLAST your contigs against the reference database

	blastn -query final.contigs.fa -db reference.fasta -max_target_seqs 5 -outfmt '6 qaccver saccver pident qcovs length qlen mismatch gapopen qstart qend sstart send' -out contigs.br

> `blastn` searches are for nucleotide queries against a nucleotide database. <br>
> `-max_target_seqs` specifies the maximum number of hits to return. <br>
> `-outfmt 6` returns the results in tab-delimited format. The information proceeding the 6 specifies the the fields to return. <br>
> We are outputting the query accession, the subject accession, the percent identity between query and subject, the percent coverage of the query by the hit, the alignment legnth, the query length,
the number of mismatches between query and subject, the number of gaps, and the start and end of the alignment in the query and subject, respectively.<br>

Based on your BLAST results:
*  are all genes represented in your contigs file?
*  is the tophit from the same species if you have multiple genes present?
*  can you confidently assign a species to your sample given the PID, qcov, qlen, etc?


## BLAST your contigs using the BLAST webpage and the nucleotide core (core_nt) database

Navigate to the [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) webpage.<br>
Notice that nucleotide core is the default database.<br>
Copy and paste your contigs into the the Query sequence section and execute the blast search.<br> 
You can print your contigs to STDOUT with the `cat` command.

	cat final.contigs.fa

Based on your BLAST results:
*  are all of your contigs from _Cryptosporidium_?
*  do these results agree with your command-line blast results?

<br>

Let's see whether the phylogenetic placement of our sequences agrees with the results from BLAST searches.

<br>

## Align your contigs with curated _Cryptosporidum_ actin, Hsp70, and ssu genes.

For alignment purposes, we will copy pre-constructed actin, Hsp70, and ssu datasets from the BMS500-2024 directory to your home directories.
This selection of _Cryptosporidium_ sequences was filtered for length and mis-identified sequences (typed as the incorrect species) were removed.

Make a directory for your alignments (assuming you are in your home directory).

	mkdir alignments

Copy the actin, hsp70, and ssu fasta files from BMS500-2024/alignments to your alignments directory (assuming you're in your home directory).

	cp ../BMS500-2024/alignments/*fasta alignments

 * What is the `*` doing here?

Extract the contigs assigned as _Cryptosporidium_ actin, Hsp70, and ssu (based on your BLAST results) from your `contigs.final.fa` file and append them to their respective multi-fasta file in the `alignments` directory. <br>
We will use a small Python script that I wrote to extract each contig and append it to the proper fasta file.  The script requires a file with the name of the contig to extract and the fasta file, from which to extract it.<br>
An example of the workflow is provided below, assuming that the contig 'k119_0' was identified as a _Cryptosporidium_ actin gene.

	cd assembly
	echo 'k119_0' > seqlist.txt
 	gi_fastasampler.py seqlist.txt final.contigs.fa >> ../alignments/actin.fasta

> `echo` simply returns the thing that we've echoed to STDOUT.  Here we are echoing the name of the 'k119_0' contig and writing it to a file called 'seqlist.txt'<br>
> `gi_fastasampler.py` is the Python script that will extract the fasta entry corresponding to the accession(s) in the seqlist.txt file.<br>
>  Remember that `>>` adds/appends information to a file instead of overwriting it.<br>
* Perform the same set of operations for your assembled Hsp70 and ssu contigs (if you have them).

Align the sequences in each fasta file with [mafft](https://mafft.cbrc.jp/alignment/software/).

	cd ~/alignments
 	mafft actin.fasta > actin.aln.fasta
  	mafft hsp70.fasta > hsp70.aln.fasta
   	mafft ssu.fasta > ssu.aln.fasta

> The tilde (~) is a shorthand way of representing your home directory.  It stands for `/home/your_account`.  Previously, we specified a path relative to where we were.  Here, we are specifying the absolute path to `alignments`
> Instead of issuing three separate `mafft` commands, we could be more efficient by using a BASH for-loop.

## Bash for loops


Bash for loops are basically little shell scripts that can be excecuted from the command line (Bash is the command-line language we are using). Like all loops, they allow you to automate iterative processes. For example, instead of opening 200 hundred fasta files and manually changing the definition lines in each, I can run a for loop that will open each fasta file and make the changes that I specify.

<br>

The basic syntax is:

	for FILE in *common_file_ending; do command $FILE; done

> The interpretation of this code is: <br>
> For every file that ends in some common ending (such as .txt or .gz), perform (do) some command on that file until there are no more files on which to operate, whereby “done” will exit us from the loop. <br>
> The $ in front of FILE indicates that $FILE is a variable, a placeholder which is referring to each file that enters the loop, just as x is a variable that represents 2 in the equation x = 2. <br>
> The `for`, `in`, `do`, and `done` are required parts of the for-loop syntax.

<br>

Use a 'for loop' to count the number of sequences in each alignment file:

	for FILE in *aln.fasta; do grep -c ">" $FILE; done


## Generate ML trees from your actin, Hsp70 and ssu alignments

We are technically generating geneologies and assuming they will be equivalent to a species phylogeny.  However, there are several reasons a gene tree may not recapitulate a species tree, including introgression, gene duplication (paralogy), and incomplete lineage sorting.

The basic command for running iqtree requires an alignment file and the specification of a [substitution model](http://www.iqtree.org/doc/Substitution-Models).<br>
If we do not know the best-fitting substituion model to employ, we can have iqtree determine it with ModelFinder.

	iqtree -s actin.aln.fasta -bb 1000

> Specifying only the alignment file will run ModelFinder first, select the best-fitting model, and then use it for phylogenetic reconstruction. <br>
> Provide a model with the `-m` argument. <br>
> `bb` indicates the number of [ultrafast bootstrap](http://www.iqtree.org/doc/Tutorial#assessing-branch-supports-with-ultrafast-bootstrap-approximation) replicates to perform. Bootstrapping resamples the data with replacement to determine how many times the relationships observed in the ML tree are recovered, which is a measure of support. <br>
> There are several files produced from an iqtree analysis.  The phylogeny is represented in Newick format (parenthetical notation) and has the '.treefile' ending. <br>
> FigTree will convert the newick text file into a graphical representation.
* Generate trees for hsp70 and ssu.
* Import your trees into FigTree (this will require me to transfer your files to the bucket for download to your computers).
* Does the phylogenetic placement and branch lengths of your sequences support a species identification? Do the bootstrap values for these relationships lend confidence to your conclusion?
* Are your results consistent with those from blast?


