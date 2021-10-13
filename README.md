# genome_variation_tracker
Program to track changes in genomes of rapidly evolving species. Specifically made to track changes in Covid 19 genomes from NCBI.

The program currently expects two files:

* sequences.csv
* sequences.fasta

The sequences.csv file is a comma separated file where the first line is a
heading line that that defines the columns for each of the following
information. One column must act as a unique identifier for the data - use
Accession for this column. This file is used to define metadata that may be of
interest in tracking the genomes. For example, if you want to track location
and time, then you need to include the Accession, Geo_Location, and Collection_Date in this line.
The line of sequences.csv that I have used is:

Accession,Length,Sequence_Type,Nuc_Completeness,Geo_Location,Collection_Date

The other file is sequences.fasta. This is a fasta file collected from NCBI.
The genomes from this file are connected to the metadata via the Accession
number.

genome_track.py reads the sequences.csv file and the sequences.fasta file. It
then creates a phylogenetic tree of the sequences. and saves this tree. Genomes
are added to the tree one at a time, so if the tree file exists, rather than
recreating the tree from scratch, the program loads the current tree and then
adds missing genomes to it. This way updating sequences.fasta and sequences.csv
does not make you start the computations over (as long as no changes have been
made to genomes already in the tree.)

Once the tree is produced, you can follow paths within the tree and call out
specific metadata about the genomes in that path. For example, if you set a
root and trace from a leaf to a root, you can see the Geo_Locations that the
genomes passed through to get from the root to that individual genome. Thus
this program could be used to track how a genome spread geographically across
the country.
