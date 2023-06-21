import time
import logging
import os
from Bio import SeqIO

from installed_clients.AssemblyUtilClient import AssemblyUtil

class CreateFasta:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']

    # -----------------------------------------------------------------
    #   Split a Sequence into 50 column chunks
    #   Return a list of lists called fasta (not the string used in the past)
    #
    def splitSequence(self, seq):
        colsz = 50
        start = 0
        lenseq = len(seq)
        fasta = []

        while True:
            end = start + colsz
            if end > lenseq:
                end = lenseq
            fasta.append([seq[start:end]])
            start += colsz
            if start > lenseq:
#                False
                break
        return fasta

    # -----------------------------------------------------------------
    #    Get the assembly sequence associated with the input genome
    #    Return a list called rpt_list (not the string used in the past)
    #
    def get_assembly_sequence(self,input_ref):
        # Download the input data as a Fasta
        # We can use the AssemblyUtils module to download a FASTA file from our Assembly data object.
        # The return object gives us the path to the file that was created.
        logging.info('Downloading Assembly data as a Fasta file.')
        assemblyUtil = AssemblyUtil(self.callback_url)
        fasta_file = assemblyUtil.get_assembly_as_fasta({'ref': input_ref})

 #       rpt_string = ''
        rpt_list = []
        for seq_record in SeqIO.parse(fasta_file['path'], 'fasta'):
            rpt_list.extend([[">" + seq_record.id ]])
            rpt_list.extend(self.splitSequence(str(seq_record.seq)))
            
#       Return a list of lists and not a string
        return rpt_list
        
    # -----------------------------------------------------------------
    #    Create a protein Fasta file for a genome
    #    Return a list called rpt_list (not the string used in the past)
    #
    def create_fasta_from_features(self, pyStr):
        myFeat = pyStr
        rpt_list = []
        for feat in myFeat:
            if 'function' not in feat:
                feat['function'] = 'unknown'

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            if 'type' in feat and feat['type'] not in ['CDS', 'gene']:
                continue

            if ('protein_translation' in feat):
                rpt_list.extend([[">" + feat['id'] + " " + feat['function'] +" (len=" + str(feat['protein_translation_length']) + ")" ]])
                rpt_list.extend(self.splitSequence(feat['protein_translation']))
                
 #       Return a list of lists and not a string
        return rpt_list

    # -----------------------------------------------------------------
    #    Create a Fasta file of the features/cdss dna_sequences in the genome
    #    Before version 9 genomes there was a single function. V9 and later it was a list called functions
    #    Before version 9 there was no cdss and the features needed a type. Values of interest are 'CDS' or 'gene'
    #    Return a list called rpt_list (not the string used in the past)
    #
    def create_fasta_from_mRNA(self, myFeat):
        rpt_list = []
        for feat in myFeat:
            if 'function' not in feat:
                feat['function'] = 'unknown'
                
            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            if 'type' in feat and feat['type'] not in ['CDS', 'gene']:
                continue

            if ('dna_sequence' in feat):
                rpt_list.extend([[">" + feat['id'] + " " + feat['function'] + " (len=" + str(feat['dna_sequence_length']) + ")" ]])
                rpt_list.extend(self.splitSequence(feat['dna_sequence']))

 #       Return a list of lists and not a string
        return rpt_list


    # -----------------------------------------------------------------
    #    Create a Fasta file for a genome
    # ######## NOT WRITTEN YET ######################
    def create_fasta_from_assembly(self, pyStr):
        rpt_list = []
        for feat in myFeat:
            if 'function' not in feat:
                feat['function'] = 'unknown'

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            if 'type' in feat and feat['type'] not in ['CDS', 'gene']:
                continue

            if ('protein_translation' in feat):
                rpt_list.extend([">" + feat['id'] + " " + feat['function'] + " (len=" + str(feat['protein_translation_length']) + ")" ])
                rpt_list.extend(self.splitSequence(feat['protein_translation']) )
                
        return rpt_list
