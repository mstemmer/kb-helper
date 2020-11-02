Please install kb-python and kallisto before you proceed.
https://www.kallistobus.tools/downloads
https://pachterlab.github.io/kallisto/download

Just run kb_wrapper.py and it will create the file structure for you in your home directory.

Paste (or symlink) fastq files in /kb_data/fastq/<sample_id>/
Make sure that sequencing lanes are merged. So there is only one R1 and one R2 file in the folder.

Create a samplesheet.tsv and place it into kb_data.
sample = provide name of sequenced sample
reference = specify the reference you want to align to
version = specify scSeq experiment version (see kallisto manual for more info)

Example of samplesheet.tsv:
sample  dr_RGC_adult_s17    dr_pineal_s1
reference       D_rerio.GRCz11.101_mt      D_rerio.GRCz11.101
version 10xv3   10xv3

Index generation:
Custom index from BioMart:
--> Place <reference>.fa.gz into kb_data/ref-seqs/ & add <reference> to reference list
--> this tool will then create the index for you
--> only works correctly, if fasta headers follow this BioMart structure:
>transcriptID|geneID|geneName|cDNA
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

Use genomic data:
To create index and tr2g files from genomic data, use 'kb ref' (see kallisto manual for more info)
and place the generated files into kb_data/index/

Make sure that the index files are in this format: <reference>.idx (e.g. D_rerio.GRCz11.101.idx)

