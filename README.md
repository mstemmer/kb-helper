# kb-helper

This little helper was created to further simplify scSeq data alignment with kallisto and bustools. 
It runs on pre-made index files or creates the necessary index and transcript_to_gene (tr2g) conversion table for you based on a custom fasta input.
In addition, it batch processes several samples for you, based on a samplesheet and outputs gene names instead of gene IDs for more intuitive downstream processing.  

## How to install and configure:  
Please install kb-python and kallisto before you proceed:  
https://www.kallistobus.tools/downloads  
https://pachterlab.github.io/kallisto/download

Now just run kb-helper.py and it will create the file structure for you in your home directory.  
See also `kb-helper.py --help`. 


## How to prepare the samplesheet:  
Simply create a TSV file following this structure:  
```
sample  <sample_id> # provide name of sequenced sample: <sample_id>
reference <reference> # specify the reference you want to align to
version <version> # specify scSeq experiment version (see kallisto manual for more info)
``` 
Specify as many experiments as you want, just by adding new columns.  
An example of a samplesheet.tsv can be found in the repository.

## How to prepare the index:  
**Custom index from BioMart:**  
* Place \<reference>.fa.gz into kb_data/ref-seqs/ & add <reference> to reference list
* this tool will then create the index and tr2g table for you
* only works correctly, if fasta headers follow this BioMart structure:
`>transcriptID|geneID|geneName|cDNA`

**Use genomic data:**  
To create index and tr2g files from genomic data, use 'kb ref' (see kallisto manual for more info)
and place the generated files into kb_data/index/

Make sure that the index files are in this format: \<reference>.idx (e.g. D_rerio.GRCz11.101.idx)  
  and the tr2g table in this: tr2g_\<reference>.tsv (e.g. tr2g_D_rerio.GRCz11.101.tsv) 

## How to prepare the input:  
Paste (or symlink) fastq files in /kb_data/fastq/<sample_id>/
Make sure that sequencing lanes are merged. So there is only one R1 and one R2 file in the folder (e.g. use cat).

## How to run:  

```
python kb-helper.py --samplesheet samplesheet.tsv
```

## Options:
By default, this helper changes gene IDs to gene names at the end of the pipeline. This can be switched off:

```
python kb-helper.py --help
usage: kb-helper [-h] [--samplesheet] [-t] [-m] [--gene_names_off]

optional arguments:
  -h, --help        show this help message and exit
  --samplesheet     Please provide path to samplesheet (tsv format)
  -t , --threads    Number of threads to use. Defaults to 32
  -m , --memory     Memory to use in GB. Defaults to 16
  --gene_names_off  Do not change gene_ID with gene_name.

```

## Importing the output into Seurat

The repository also provides an R notebook to conveniently import the alignment into Seurat.  
In the notebook, just specify there again the sample_ids and references like this:

```
expers <- ("<sample_id>")
index <- ("<reference>")
```

In the second part:  
* the data will be imported (using the same folder structure as above) 
* empty droplet removal will be performed (with dropletutils and defining the inflection point in the knee plot)
* a Seurat object will be created and saved in your kb_data/out/ folder under the respective reference folder

The code on empty droplet removal was imported from the [kallisto vignettes](https://www.kallistobus.tools/tutorials)


## Notes
* kb-helper is in an early state and was primarily written to speed up and to channel multiple scSeq experiments, but I hope that it can be useful to a wider community
* so far it was only tested on Linux (Ubuntu).
