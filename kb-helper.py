from kallisto_bustools import KallistoBustools
import time
import os
import sys
import os.path
import pandas as pd
import argparse

''' Please install kb-python and kallisto before you proceed.
https://www.kallistobus.tools/downloads
https://pachterlab.github.io/kallisto/download

Just run kb-helper.py and it will create the file structure for you in your home directory.

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
'''

class RunAlignment():
    def __init__(self):
        parser = argparse.ArgumentParser(prog='kb-helper')
        parser.add_argument('--samplesheet', dest='samples', metavar='', help='Please provide path to samplesheet (tsv format)')
        parser.add_argument('-t', '--threads', metavar='', dest='t', default=32, type=int, help='Number of threads to use. Defaults to 32')
        parser.add_argument('-m', '--memory', metavar='', dest='m', default=16, type=int, help='Memory to use in GB. Defaults to 16')
        parser.add_argument('--gene_names_off', dest='gn_off', action='store_true', default=False, help='Do not change gene_ID with gene_name.')
        self.args = parser.parse_args()
        # print(self.args.gn_off)

    def setup_dir(self): # Setup working directory
        os.chdir(os.path.expanduser("~")) #change to /home/$USER/
        os.makedirs('kb_data/index', exist_ok=True)
        os.makedirs('kb_data/ref-seqs', exist_ok=True)
        os.makedirs('kb_data/fastq', exist_ok=True)
        os.makedirs('kb_data/out/', exist_ok=True)
        os.chdir('kb_data')
        self.path = os.getcwd()
        print('--> Set up working directory: ' + self.path)

    def read_samplesheet(self):
        print('--> Reading samplesheet')
        try:
            self.sample_sheet = pd.read_table(self.args.samples, header=None, index_col=0)
        except (FileNotFoundError, AttributeError, ValueError):
            print('--> Samplesheet not found. Please provide a valid samplesheet in tsv format with --samplesheet <path>')
            sys.exit("--> Exiting program")
        print(self.sample_sheet)
        print('--> Number of samples: ' + str(len(self.sample_sheet.columns)))

    def align(self):
        for s in range(len(self.sample_sheet.columns)):
            print('--> Working on sample ' + str(s+1))
            sample = self.sample_sheet.iloc[0,s]
            reference = self.sample_sheet.iloc[1,s]
            version = self.sample_sheet.iloc[2,s]

            os.makedirs('out/' + reference, exist_ok=True)

            kb = KallistoBustools(sample, reference, version)
            kb.check_ref()
            kb.mod_tr2g() # checks trg2 file, if there are NA or blank fields in gene_name column of tr2g file
            kb.count(self.args.t, self.args.m) # pseudoalignemnt; default is threads=16, memory=16G
            if self.args.gn_off == False:
                kb.mod_genes() # OPTIONAL, switches gene_name with gene_id in the matrix files. If used, gene_names will be seen in seurat instead of ENSEMBL IDs


if __name__ == '__main__':
    start_time = time.time()

    align = RunAlignment()
    align.setup_dir()
    align.read_samplesheet()
    align.align()

    print("Runtime: {:.2f} minutes".format((time.time()-start_time)/60))
