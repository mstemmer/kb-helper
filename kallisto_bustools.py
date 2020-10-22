import time
import os
import os.path
import sys
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

class KallistoBustools():
    ''' This class is meant to run kallisto and bustools in a more interactive manner.
    It will create the working directories and index for you.
    In a second step, specific chromium samples will be pseudoaligned to the specified reference.
    The tr2g gene file will be generated using the fasta headers, if provided as BioMart exports.
    The tool also detects NA genes and uses IDs instead. If wanted, the mod_genes method adjusts the matrix files
    to see gene names instead of IDs in downstream processing tools.
    '''
    def __init__(self, sample_id, ref_species, version):
        self.ref_name = ref_species
        self.sample_id = sample_id
        self.version = version
        print('--> Initializing kallisto')
        print('--> Parameters: Sample= ' + self.sample_id + ' Index= ' + self.ref_name + ' Version= ' + self.version)

    def show_input(self):
        return self.ref_name

    def check_ref_exists(self):
        if os.path.exists('index/' + self.ref_name + '.idx') and os.path.exists('index/tr2g_' + self.ref_name + '.tsv'):
            return True

    def prep_tr2g_fasta(self): # only works on BioMart exports with this header >transcriptID|geneID|geneName
        print('--> Reference fasta file found: adapting headers')
        os.system('gunzip -f ref-seqs/' + self.ref_name + '.fa.gz')
        os.system('sed -i \'s/|/ /g\' ref-seqs/' + self.ref_name + '.fa')
        print('--> Extract tr2g file from fasta headers')
        os.system('grep -e ">" ref-seqs/' + self.ref_name + '.fa | cut -c2- > ref-seqs/' + self.ref_name + '_head.tsv') # extract all fasta headers
        tr2g = pd.read_table('ref-seqs/' + self.ref_name + '_head.tsv', header=None) # Read cdna table and output only transcript ids
        tr2g = tr2g[0].str.split(' ', expand=True)
        tr2g = tr2g[[0, 1, 2]]
        print('--> Inspecting tr2g file')
        print(tr2g)
        tr2g.to_csv('index/tr2g_' + self.ref_name + '.tsv', sep='\t', index=False, header=False)
        os.system('gzip -f ref-seqs/' + self.ref_name + '.fa')

    def check_ref(self):
        if self.check_ref_exists() != True:
            if os.path.exists('ref-seqs/' + self.ref_name + '.fa.gz'):
                print('--> No index found. Checking for BioMart compatible reference fasta file...')
                self.prep_tr2g_fasta()
                print('--> Creating index')
                os.system('kallisto index -i index/' + self.ref_name + '.idx' + ' ref-seqs/' + self.ref_name + '.fa.gz' + ' --make-unique')
            else:
                print('--> No reference index found')
                print('--> No reference fasta file found')
                print('--> Please provide either <reference>.fa.gz file in ref-seq/')
                print('    OR pre-built <reference>.idx and tr2g_<reference>.tsv in index/')
                print('--> See also readme file for help')
                sys.exit("--> Exiting program")
        else:
            print('--> Reference index and tr2g file found. Continue to alignment.')

    def remove_version(self):
        pass

    def mod_tr2g(self):
        print('--> Checking tr2g file. If name = NA or blank, then gene_name = gene_id')
        tr2g = pd.read_table('index/tr2g_' + self.ref_name + '.tsv', header=None)
        # na=True is needed because 'NA'doesnt seem to be recognied...but then also works with blank fields
        tr2g_na = tr2g[tr2g[2].str.match('NA', na=True)]
        tr2g_na[2] = tr2g_na[1] # Set column 2 equal to column 1
        tr2g_rest = tr2g[~tr2g[2].str.match('NA', na=True)] # Extract all non_NA rows
        tr2g_mod = pd.concat([tr2g_rest, tr2g_na]) # Combine the two dfs
        tr2g_mod = tr2g_mod.sort_index() # sort them according to their index
        print(tr2g_mod)
        tr2g_mod.to_csv('index/tr2g_' + self.ref_name + '_mod.tsv', sep='\t', index=False, header=False)

    def count(self, threads, memory):
        try:
            if os.path.exists('out/' + self.ref_name + '/' + self.sample_id):
                print('--> Alignment already exists.')
                print('Do you want to run it again? (yes/no)')
                run_input = input()
                if run_input == 'yes':
                    os.system('kb count -t ' + str(threads) + ' -m ' + str(memory) + ' -i index/' + self.ref_name + '.idx' + ' -g index/tr2g_' + self.ref_name + '_mod.tsv' \
                    + ' -x '+ str(self.version) + ' --h5ad --overwrite' + ' -o out/' + self.ref_name + '/' + self.sample_id \
                    + ' fastq/' + self.sample_id + '/' + '*R1*.fastq.gz' + ' fastq/' + self.sample_id + '/' + '*R2*.fastq.gz')
                elif run_input == 'no':
                    print('--> Continue')
            else:
                print('--> No alignment found. Creating new alignment')
                os.system('kb count -t ' + str(threads) + ' -m ' + str(memory) + ' -i index/' + self.ref_name + '.idx' + ' -g index/tr2g_' + self.ref_name + '_mod.tsv' \
                + ' -x '+ str(self.version) + ' --h5ad --overwrite' + ' -o out/' + self.ref_name + '/' + self.sample_id \
                + ' fastq/' + self.sample_id + '/' + '*R1*.fastq.gz' + ' fastq/' + self.sample_id + '/' + '*R2*.fastq.gz')
        except FileNotFoundError:
            print('Please check your fastq files. Are they named correctly and in the correct folder? <sample_name>')
            sys.exit("--> Exiting program")

    def mod_genes(self):
        try:
            print('--> Setting gene_id to gene_name in matrix files.')
            tr2g_mod = pd.read_table('index/tr2g_' + self.ref_name + '_mod.tsv', header=None)
            tr2g_mod.drop_duplicates(subset = 1, keep = "first", inplace= True)
            tr2g_mod[2].to_csv('out/' + self.ref_name + '/' + self.sample_id + '/' + 'counts_unfiltered/cells_x_genes.genes.txt', index=False, header=False)
        except FileNotFoundError:
            print('--> Please check your fastq files. Are they named correctly and in the correct folder? <sample_name>')
            sys.exit("--> Exiting program")
