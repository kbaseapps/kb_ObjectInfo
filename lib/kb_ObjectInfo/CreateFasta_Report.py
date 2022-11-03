import time

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class CreateFasta:
    def __init__(self, config):
        self.config = config

    # -----------------------------------------------------------------
    #   Split a Sequence into 50 column chunks
    #
    def splitSequence(self, seq):
        colsz = 50
        start = 0
        lenseq = len(seq)
        line = ""

        while True:
            end = start + colsz
            if end > lenseq:
                end = lenseq
            # print seq[start:end]
            line += seq[start:end] + "\n"
            start += colsz
            if start > lenseq:
#                False
                break
        return line


    # -----------------------------------------------------------------
    #    Create a protein Fasta file for a genome
    #
    def create_Fasta_from_features(self, pyStr):
        myFeat = pyStr
        line = ""
        for feat in myFeat:
            if 'function' not in feat:
                feat['function'] = 'unknown'

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            if 'type' in feat and feat['type'] not in ['CDS', 'gene']:
                continue

            if ('protein_translation' in feat):
                line += ">" + feat['id'] + " " + feat['function']
                line += " (len=" + str(feat['protein_translation_length']) + ")" + "\n"

                # print line
                line += self.splitSequence(feat['protein_translation']) + "\n"
        return line

    # -----------------------------------------------------------------
    #    Create a protein Fasta file for a genome
    #
    def create_Fasta_from_mRNA(self, myFeat):
        line = ""
        for feat in myFeat:
            if 'function' not in feat:
                feat['function'] = 'unknown'

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            if 'type' in feat and feat['type'] not in ['CDS', 'gene']:
                continue

            if ('dna_sequence' in feat):
                line += ">" + feat['id'] + " " + feat['function']
                line += " (len=" + str(feat['dna_sequence_length']) + ")" + "\n"

                # print line
                line += self.splitSequence(feat['dna_sequence']) + "\n"
        return line


    # -----------------------------------------------------------------
    #    Create a Fasta file for a genome
    # ######## NOT WRITTEN YET ######################
    def create_Fasta_from_assembly(self, pyStr):
        line = ""
        for feat in myFeat:
            if 'function' not in feat:
                feat['function'] = 'unknown'

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            if 'type' in feat and feat['type'] not in ['CDS', 'gene']:
                continue

            if ('protein_translation' in feat):
                line += ">" + feat['id'] + " " + feat['function']
                line += " (len=" + str(feat['protein_translation_length']) + ")" + "\n"

                # print line
                line += self.splitSequence(feat['protein_translation']) + "\n"
        return line
