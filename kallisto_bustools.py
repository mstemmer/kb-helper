import time
import os
import os.path
import sys
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

class KallistoBustools():
    def __init__(self, sample_id, ref_species, version):
        self.ref_name = ref_species
        self.sample_id = sample_id
        self.version = version
        print('--> Initializing kallisto')
        print(f'--> Parameters: Sample= {self.sample_id} Index= {self.ref_name} Version= {self.version}')

    def check_ref_exists(self):
        if os.path.exists(f'index/{self.ref_name}.idx') and os.path.exists(f'index/tr2g_{self.ref_name}.tsv'):
            return True

    def prep_tr2g_fasta(self): # only works on BioMart exports with this header >transcriptID|geneID|geneName
        print('--> Reference fasta file found: adapting headers')
        os.system(f'gunzip -f ref-seqs/{self.ref_name}.fa.gz')
        os.system(f'sed -i \'s/|/ /g\' ref-seqs/{self.ref_name}.fa') #change delimiter from '|' to ''
        print('--> Extracting tr2g file from fasta headers')
        os.system(f'grep -e ">" ref-seqs/{self.ref_name}.fa | cut -c2- > ref-seqs/{self.ref_name}_head.tsv') # extract all fasta headers
        tr2g = pd.read_table(f'ref-seqs/{self.ref_name}_head.tsv', header=None) # Read header table
        tr2g = tr2g[0].str.split(' ', expand=True)
        tr2g = tr2g[[0, 1, 2]]
        print('--> Inspecting tr2g file')
        print(tr2g)
        tr2g.to_csv(f'index/tr2g_{self.ref_name}.tsv', sep='\t', index=False, header=False)
        os.system(f'gzip -f ref-seqs/{self.ref_name}.fa')

    def check_ref(self):
        if self.check_ref_exists() != True:
            if os.path.exists(f'ref-seqs/{self.ref_name}.fa.gz'):
                print('--> No index found. Checking for BioMart compatible reference fasta file...')
                self.prep_tr2g_fasta()
                print('--> Creating index')
                os.system(f'kallisto index \
                -i index/{self.ref_name}.idx \
                ref-seqs/{self.ref_name}.fa.gz --make-unique')
            else:
                print('--> No reference index found')
                print('--> No reference fasta file found')
                print('--> Please provide either <reference>.fa.gz file in ref-seq/')
                print('    OR pre-built <reference>.idx and tr2g_<reference>.tsv in index/')
                print('--> See also readme file for help')
                sys.exit("--> Exiting program")
        else:
            print('--> Reference index and tr2g file found. Continue to alignment.')

    def mod_tr2g(self):
        ''' When gene name is NA or empty, Gene ID is copied to this field.'''
        print('--> Checking tr2g file. If name = NA or blank, then gene_name = gene_id')
        tr2g = pd.read_table(f'index/tr2g_{self.ref_name}.tsv', header=None)
        tr2g_na = tr2g[tr2g[2].str.match('NA', na=True)] # na=True is needed because 'NA'doesnt seem to be recognied...but then also works with blank fields
        tr2g_na[2] = tr2g_na[1] # Set column 2 equal to column 1
        tr2g_rest = tr2g[~tr2g[2].str.match('NA', na=True)] # Extract all non_NA rows
        tr2g_mod = pd.concat([tr2g_rest, tr2g_na]) # Combine the two dfs
        tr2g_mod = tr2g_mod.sort_index() # sort them according to their index
        print(tr2g_mod)
        tr2g_mod.to_csv(f'index/tr2g_{self.ref_name}_mod.tsv', sep='\t', index=False, header=False)

    def run_kb(self, threads, memory):
        os.system(f'kb count -t {threads} -m {memory} \
        -i index/{self.ref_name}.idx \
        -g index/tr2g_{self.ref_name}_mod.tsv \
        -x {str(self.version)} \
        -o out/{self.ref_name}/{self.sample_id} \
        --h5ad \
        fastq/{self.sample_id}/*R1*.fastq.gz fastq/{self.sample_id}/*R2*.fastq.gz')

    def count(self, threads, memory):
        try:
            if os.path.exists(os.path.join('out', self.ref_name, self.sample_id)):
                print('--> Alignment already exists.')
                print('Do you want to run it again? (yes/no)')
                run_input = input()
                if run_input == 'yes':
                    self.run_kb(threads, memory)
                elif run_input == 'no':
                    print('--> Continue')
            else:
                print('--> No alignment found. Creating new alignment')
                self.run_kb(threads, memory)
        except FileNotFoundError:
            print(r'--> Please check your fastq files. Are they named correctly and in the correct folder? "<sample_name>/*R1*.fastq.gz"')
            print(r'                                                                                       "<sample_name>/*R2*.fastq.gz"')
            sys.exit("--> Exiting program")

    def mod_genes(self):
        ''' Changes gene_name and gene_id fields to work only with gene_names in downstream processes
        This is optional and can be switched off.
        '''
        try:
            print('--> Setting gene_id to gene_name in matrix files.')
            tr2g_mod = pd.read_table(f'index/tr2g_{self.ref_name}_mod.tsv', header=None)
            tr2g_mod.drop_duplicates(subset = 1, keep = "first", inplace= True)
            tr2g_mod[2].to_csv(os.path.join('out', self.ref_name, self.sample_id, 'counts_unfiltered', 'cells_x_genes.genes.txt'), index=False, header=False)
        except FileNotFoundError:
            print(r'--> Please check your fastq files. Are they named correctly and in the correct folder? "<sample_name>/*R1*.fastq.gz"')
            print(r'                                                                                       "<sample_name>/*R2*.fastq.gz"')
            sys.exit("--> Exiting program")
