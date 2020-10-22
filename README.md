Please install kb-python and kallisto before you proceed.
https://www.kallistobus.tools/downloads
https://pachterlab.github.io/kallisto/download

Run kb_wrapper.py and it will create the folder structure for you in your home directory.

Paste fastq files in /kb_data/fastq/<sample_name>/. 
Make sure that sequencing lanes are merged. So there is only one R1 and one R2 file in the folder.

Create a samplesheet.tsv and place it into kb_data.  
sample = provide name of sequenced sample <sample_name>.  
reference = specify the reference you want to align to.  
version = specify scSeq experiment version (see kallisto manual for more info).  
--> See samplesheet.tsv in kb_wrapper as an example  
--> Feed the tsv file into kb_wrapper with --samplesheet

Custom index from BioMart:  
--> Place <reference>.fa.gz into kb_data/ref-seqs/ & add <reference> to reference list.  
--> this tool will then create the index for you.  
--> only works correctly, if fasta headers follow BioMart structure:   
Fasta Header:  
\>transcriptID|geneID|geneName


To create index and tr2g files from Genome data, use kb ref (see kallisto manual for more info)
and place the generated files into kb_data/index/




