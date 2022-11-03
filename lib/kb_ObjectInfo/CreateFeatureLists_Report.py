import time
import os
import re
from .CreateFasta_Report import CreateFasta
from pprint import pprint, pformat

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class CreateFeatureLists:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        (self.cats, self.cat2name, self.cat2group, self.domfam2cat, self.domfam2name, self.domfam2ns) = self._configure_categories()

    def _configure_categories(self):

        domain_desc_basepath = os.path.abspath('/kb/module/data')
        domain_to_cat_map_path = dict()
        domain_cat_names_path = dict()
        domain_fam_names_path = dict()
        domain_to_cat_map_path['COG'] = os.path.join(domain_desc_basepath, 'COG_2014.tsv')
        domain_cat_names_path['COG'] = os.path.join(domain_desc_basepath, 'COG_2014_funcat.tsv')
        domain_fam_names_path['COG'] = os.path.join(domain_desc_basepath, 'COG_2014.tsv')
        domain_to_cat_map_path['PF'] = os.path.join(domain_desc_basepath, 'Pfam-A.clans.tsv')
        domain_cat_names_path['PF'] = os.path.join(domain_desc_basepath, 'Pfam-A.clans_names.tsv')
        domain_fam_names_path['PF'] = os.path.join(domain_desc_basepath, 'Pfam-A.clans.tsv')
        domain_to_cat_map_path['TIGR'] = os.path.join(domain_desc_basepath, 'TIGRInfo.tsv')
        domain_cat_names_path['TIGR'] = os.path.join(domain_desc_basepath, 'tigrrole2go.txt')
        #domain_fam_names_path['TIGR']  = os.path.join(domain_desc_basepath, 'tigrfams2go.txt')
        domain_fam_names_path['TIGR'] = os.path.join(domain_desc_basepath, 'TIGRInfo.tsv')
        domain_to_cat_map_path['SEED'] = os.path.join(domain_desc_basepath, 'SEED_subsys.txt')
        domain_cat_names_path['SEED'] = os.path.join(domain_desc_basepath, 'SEED_funcat.txt')
        #domain_cat_names_path['SEED']  = os.path.join(domain_desc_basepath, 'SEED_subsys.txt')
        domain_fam_names_path['SEED'] = os.path.join(domain_desc_basepath, 'SEED_subsys.txt')
        
        cats = []
        cat2name     = dict()
        cat2group    = dict()
        domfam2cat   = dict()
#        cat2domfams  = dict()
#        domfam2group = dict()
        domfam2name  = dict()
        domfam2ns = dict()

        # read all mappings between groups and domfams
        for namespace in ['COG', 'PF', 'TIGR', 'SEED']:

            cat2name[namespace] = dict()
            cat2group[namespace] = dict()
#            domfam2cat[namespace] = dict()
#            cat2domfams[namespace] = dict()
#            domfam2group[namespace]  = dict()
#            domfam2name[namespace]  = dict()
            
            # get high-level cats
            tigrrole_id2cat = dict()
            with open(domain_cat_names_path[namespace], 'r') as dom_cat_handle:
                for line in dom_cat_handle.readlines():
                    line = line.strip()

                    if namespace == 'COG':
                        [cat, cat_group, cat_name] = line.split("\t")[0:3]
                        if  cat not in cats:
                            cats.append(cat)
                        cat2name[namespace][cat] = cat_name
                        cat2group[namespace][cat] = cat_group

                    elif namespace == 'PF':
                        [cat, cat_name] = line.split("\t")[0:2]
                        if cat not in cats:
                            cats.append(cat)
                        cat2name[namespace][cat] = cat_name
                        cat2group[namespace][cat] = None

                    elif namespace == 'TIGR':
                        if line.startswith('!'):
                            continue
                        [cat, cat_id, cat_group, cat_name_plus_go_terms] = line.split("\t")[0:4]
                        tigrrole_id2cat[cat_id] = cat
                        if cat not in cats:
                            cats.append(cat)
                        cat_name = re.sub(' *\> GO:.*$', '', cat_name_plus_go_terms)
                        cat2name[namespace][cat] = cat_name
                        cat2group[namespace][cat] = cat_group

                    elif namespace == 'SEED':
                        #[cat_group, cat_subgroup, cat, domfam] = line.split("\t")[0:4]
                        [cat_group, cat] = line.split("\t")[0:2]
                        if cat not in cats:
                            cats.append(cat)
                        #cat_disp = re.sub('_', ' ', cat)
                        #cat2name[namespace][cat] = cat_disp

            # get domfam to cat map, and vice versa
            with open(domain_to_cat_map_path[namespace], 'r') as dom2cat_map_handle:
                for line in dom2cat_map_handle.readlines():
                    line = line.strip()

                    if namespace == 'COG':
                        [domfam, cat_str, dom_name] = line.split("\t")[0:3]
                        cat = cat_str[0]  # only use first cat
                        
                    elif namespace == 'PF':
                        [domfam, cat, cat_name, dom_id, dom_name] = line.split("\t")[0:5]

                    elif namespace == 'TIGR':
                        if line.startswith('!'):
                            continue
                        [domfam_id, domfam, cat_group, cat_id, domfam_name, ec_id, dom_name] = line.split("\t")[0:7]
                        if cat_id != '' and int(cat_id) != 0 and cat_id in tigrrole_id2cat:
                            cat = tigrrole_id2cat[cat_id]
                        else:
                            continue

                    elif namespace == 'SEED':
                        [cat_group, cat_subgroup, cat, domfam] = line.split("\t")[0:4]
                        domfam = domfam.strip()
                        domfam = re.sub(' *\#.*$', '', domfam)
                        domfam = re.sub(' *\(EC [\d\.\-\w]*\) *$', '', domfam)
                        domfam = re.sub(' *\(TC [\d\.\-\w]*\) *$', '', domfam)
                        domfam = re.sub('_', ' ', domfam)
                        #domfam = 'SEED ' + domfam
                        dom_name = domfam
                        dom_name = cat_subgroup
                        cat2group[namespace][cat] = cat_group
                        
                        #if 'DNA-directed RNA polymerase' in domfam:
                        #    print ("DOMFAM ", domfam + " SUBGROUP ", cat_subgroup, " CAT ", cat)

                    domfam2ns[domfam] = namespace
                    domfam2cat[domfam] = cat
                    domfam2name[domfam]  = dom_name
                    
 #                   if cat not in cat2domfams[namespace]:
 #                       cat2domfams[namespace][cat] = []
 #                   cat2domfams[namespace][cat].append(domfam)
                    
 #                   if cat in cat2group[namespace]:
 #                       domfam2group[namespace][domfam]  = cat2group[namespace][cat]
 #                   else:
 #                       domfam2group[namespace][domfam]  = None
 
        return(cats, cat2name, cat2group, domfam2cat, domfam2name, domfam2ns)
        
    # -----------------------------------------------------------------
    #    Create a Delimited Table version of the genes in a genome
    #

    def delimitedTable(self, genome, format, features):

#        seed_basepath = os.path.abspath('/kb/module/data')
#        seed_subsys   = os.path.join(seed_basepath, 'subsys.txt')

#        seed_cat = dict()

#        with open(seed_subsys, 'r') as seed_handle:
#            for line in seed_handle.readlines():
#
#                line = line.strip()
#                [cat_group, cat_subgroup, cat, seedfam] = line.split("\t")[0:4]
#                if seedfam in seed_cat:
#                    seed_cat[seedfam] += "; " + cat
#                else:
#                    seed_cat[seedfam] = cat
        seed_cat = self.domfam2cat
        
        line = ""
        lineList = ["Feature ID", "Feature type", "Contig", "Location", "Strand", "Feature function", "Aliases",
                    "RAST Functional Assignment", "RAST Functional Group 1", "RAST Functional Group 2"]
        if format == 'tab':
            line += "\t".join(lineList) + "\n"
        else:
            line += ",".join(lineList) + "\n"
            
        cat = ''
        domfam = ''
        for feat in genome[features]:
            if 'function' not in feat:
                feat['function'] = 'unknown'
            else:
                if feat['function'] in seed_cat:
                    domfam = feat['function']
                    cat = seed_cat[domfam]
                #if 'DNA-directed RNA polymerase' in feat['function']:
                #                    print ("DOMFAM ", domfam, " CAT", cat)                    

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])
                for func in feat['functions']:
                    if func in seed_cat:
                        domfam = func
                        cat = seed_cat[domfam]
                    #if 'DNA-directed RNA polymerase' in feat['function']:
                    #    print ("DOMFAM ", domfam, " CAT", cat)
                    break # Taking just the first


            aliases = ''
            if 'aliases' in feat:
                for al in feat['aliases']:
                    if isinstance(al, (str)):
                        aliases = ', '.join(feat['aliases'])
                        break
                    elif isinstance(al, (list)):
                        if ('synonym' in al[0]  or 'gene' in al[0] or 'protein_id' in al[0] ) and aliases > '     ':
                            aliases += ", " + al[1]
                        elif ('synonym' in al[0]  or 'gene' in al[0] or 'protein_id' in al[0] )  :
                            aliases += al[1]

                            
            if 'type' not in feat:
                feat['type'] = features

            location = ''
            contig = ''
            strand = ''
            if len(feat['location']) > 0:
                locList = []
                # For those REALLY rare occassions when there is more than one location in Prokaryotes
                for loc in feat['location']:
                    contig = loc[0]
                    strand = loc[2]
                    if strand == '+':
                        start = loc[1]
                        stop = loc[1] + loc[3] - 1
                    else:
                        start = loc[1]
                        stop = loc[1] - loc[3] + 1
                    locList.append(str(start) + '..' + str(stop))

                location = ", ".join(locList)
                
            catgroup = ''
            subgroup = ''
            if cat > ' ':
                if cat in self.cat2group['SEED']:
                    catgroup = self.cat2group['SEED'][cat] 
                if domfam in self.domfam2name:
                    subgroup = self.domfam2name[domfam]
 
            if format == 'tab':
                lineList = [feat['id'], feat['type'], contig, location, strand, feat['function'], aliases, cat, subgroup, catgroup]
                line += "\t".join(lineList) + "\n"
            else:
                feat['function'] = '"' + feat['function'] + '"'
                aliases = '"' + aliases + '"'
                location = '"' + location + '"'
                lineList = [feat['id'], feat['type'], contig, location, strand, feat['function'], aliases, cat, subgroup, catgroup]
                line += ",".join(lineList) + "\n"

        return line


        # -----------------------------------------------------------------
        #    Create a GFF3 version of the features in a genome
        #

    def gff3(self, genome, features):
        line = ""
        if features not in genome:
            return line
        
        for feat in genome[features]:
            if 'function' not in feat:
                feat['function'] = 'unknown'

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])

            aliases = ''
            if 'aliases' in feat:
                for al in feat['aliases']:
                    if isinstance(al,list) and len(al) == 2:
                        aliases += ":" + al[0] + "-" + al[1]
                    elif isinstance(al,list):
                        for i in al:
                            aliases += ":" + i
                    else:
                        aliases = ':'.join(al)
            if 'type' not in feat:
                feat['type'] = features

            location = ''
            contig = ''
            strand = ''
            start = 0
            stop = 0
            if len(feat['location']) > 0:
                locList = []
                # For those REALLY rare occassions when there is more than one location in Prokaryotes
                for loc in feat['location']:
                    contig = loc[0]
                    strand = loc[2]
                    if strand == '+':
                        start = loc[1]
                        stop = loc[1] + loc[3] - 1
                        locList.append(str(start) + '..' + str(stop))
                    else:
                        start = loc[1]
                        stop = loc[1] - loc[3] + 1
                        locList.append(str(start) + '..' + str(stop))

                location = ", ".join(locList)

            ph = "."  # Placeholder for missing data
            attrib = "ID=" + feat['id']
            if feat['function'] != 'unknown':
                attrib += ";FUNCTION=" + feat['function']
            if aliases > '     ':
                attrib += ";ALIASES=" + aliases

            lineList = [contig, ph, feat['type'], str(start), str(stop), ph, strand, ph, attrib]
            line += "\t".join(lineList) + "\n"

        return line



    # -----------------------------------------------------------------
    #   Domain Annotation Reports
    #

    #   OBJECT: DomainAnnotation
    #   FUNCTION: User-defined function to format all the domains for a gene
    #
    def printGeneDomain(self, contig, geneName, geneDomain, format, cutoff):
        line = ""
        lineList = ""
        for domain in geneDomain:
            list = geneDomain[domain]
            if list[0][2] < cutoff:
                if '.' in domain and len(domain) < 15:
                    domain = domain.split('.')[0]
                #print ("DOMAIN", domain)
                if domain in self.domfam2ns:
                    namespace = self.domfam2ns[domain]
                else:
                    namespace = ''
                #print ("NAMESpaCE", namespace)
                if domain in self.domfam2name:
                    dom_name = self.domfam2name[domain]
                else:
                    dom_name = ''
                #print ("DOMNAME", dom_name)
                if domain in self.domfam2cat:
                    cat = self.domfam2cat[domain]
                else:
                    cat = ''
                #print ("CAT", cat)
                
                if cat > ' ' and cat in self.cat2name[namespace]:
                    cat_name = self.cat2name[namespace][cat]
                else:
                    cat_name = ''
                #print ("CATNAME", cat_name)
                
                if cat > ' ' and cat in self.cat2group[namespace]:
                    cat_group = self.cat2group[namespace][cat]
                else:
                    cat_group = ''
                #print ("CATGROUP", cat_group)
                
                lineList = [contig, geneName, domain, str(list[0][2]), str(list[0][0]), str(list[0][1]), dom_name, cat, cat_name, str(cat_group)]
                if format == 'tab':
                    line += "\t".join(lineList)
                elif format == 'csv':
                    line += ",".join(lineList)
                #            print line
                line += "\n"
        return line


    #
    #   OBJECT: DomainAnnotation
    #   FORMAT: tab or comma delimited list of the genes, domains, e-values, and start/stop of domain hit
    #   Loop through all of the contigs and get all of the genes
    #   Uses printGeneDomain to print out individual lines
    #
    def readDomainAnnList(self, pyStr, format, cutoff):
        #   Make sure the cutoff is a number
        if not isinstance(cutoff, (int, float, complex)):
            print ("Cutoff Value must be numeric.")
            return

        # Header
        line = ""
        lineList = ["Contig", "Gene ID", "Domain", "Evalue", "Start", "Stop","Domain Name", "Category", "Category Name", "Category Group"]

        #   Check for valid formats
        if format not in ['tab', 'csv']:
            print ("Invalid format. Valid formats are tab and csv")
            return
        elif format == 'tab':
            line = "\t".join(lineList)
        elif format == 'csv':
            line = "'" + ",".join(lineList)

        # Add line-end to the header
        line += "\n"

        myData = pyStr['data']

        for contig in myData:
            contigData = myData[contig]
            for gene in contigData:
                if (gene[4]):
                    domain = gene[4]
                    line += self.printGeneDomain(contig, gene[0], domain, format, cutoff)

        return line


    #
    #   OBJECT: DomainAnnotation
    #   FUNCTION: User-defined function to count the domains for a gene
    #
    def countGeneDomain(self, contig, geneName, geneDomain, format, cutoff, myDict):
        for domain in geneDomain:
            list = geneDomain[domain]
            if list[0][2] < cutoff:
                if domain in myDict:
                    myDict[domain] += 1
                else:
                    myDict[domain] = 1

        return myDict


    #
    #   OBJECT: DomainAnnotation
    #   FORMAT: List of the domains and number of occurrences in the genome
    #   Uses countGeneDomain to get the statistics for an individual gene
    #
    def readDomainAnnCount(self, pyStr, format, cutoff):
        #   Make sure the cutoff is a number
        if not isinstance(cutoff, (int, float, complex)):
            print ("Cutoff Value must be numeric.")
            return

        # Header
        line = ""
        lineList = ["Domain", "Count"]

        #   Check for valid formats
        if format not in ['tab', 'csv']:
            print ("Invalid format. Valid formats are tab and csv")
            return
        elif format == 'tab':
            line = "\t".join(lineList)
        elif format == 'csv':
            line = "'" + ",".join(lineList)

        # Add line-end to the header
        line += "\n"

        myData = pyStr['data']
        count = 0
        myDict = {}
        for contig in myData:
            contigData = myData[contig]

            for gene in contigData:
                if (gene[4]):
                    myDict = self.countGeneDomain(contig, gene[0], gene[4], format, cutoff, myDict)

        domainList = list(myDict.keys())
        domainList.sort()
        for domain in domainList:
            lineList = [domain, str(myDict[domain])]
            if format == 'tab':
                line += "\t".join(lineList)
            elif format == 'csv':
                line += ",".join(lineList)
            line += "\n"

        return line


    #
    #   OBJECT: FeatureSet or SequenceSet
    #   FORMAT: List of the contents of the object
    #
    def readFeatSeq(self, pyStr, format):
        
        # Header
        line = ""
        lineList = list()

        cf = CreateFasta(self.config)
#
#   Type 1 - Order matters
#
        if 'description' in pyStr and 'elements' in pyStr and 'element_ordering' in pyStr:
            lineList.append(['Description', str(pyStr['description']), "\nOrdered Elements:"])
            eleOrder = pyStr['element_ordering']
            count = 1
            for index in eleOrder:
                lineList.append([str(count), index, ",".join(pyStr['elements'][index]) ])
                count += 1

#
#   Type 2 - Unordered
#
        elif 'description' in pyStr and 'elements' in pyStr:
            lineList.append(['Description', pyStr['description']])
            lineList.append(["\nUnordered Elements:"])
            myElements = pyStr['elements']
            count = 0
            for element in myElements:
                if isinstance(myElements[element],list):
                    for i in myElements[element]:
                        lineList.append([element, i])
                else:
                    lineList.append([element, myElements[element]])
                count += 1

#
#   Type 3 - With Sequences
#
        elif 'description' in pyStr and 'sequences' in pyStr and 'sequence_set_id' in pyStr:
            lineList.append(['Set Description', pyStr['description']])
            lineList.append(["Sequence Set ID", pyStr['sequence_set_id']])
            lineList.append(["Sequences:"])

            mySequences = pyStr['sequences']
            count = 0
            for seq in mySequences:
                seqline = cf.splitSequence(seq['sequence'])
                lineList.append([">" + seq['sequence_id'], seq['description']])
                lineList.append([seqline])
                count += 1

#
#   Type Unknown
#
        else:
            print ("This type of FeatureSet has not been described yet")

        #   Check for valid formats
        if format not in ['tab', 'csv']:
            print ("Invalid format. Valid formats are tab and csv")
            return
        elif format == 'tab':
            for row in lineList:
                line += "\t".join(row) + "\n"
        elif format == 'csv':
            line = "'"   ## Needed for Excel to recognize comma-delimited
            for row in lineList:
                line += ",".join(row) + "\n"

        # Add line-end to the header
        line += "\n"

        #print ("LINE: ", line)
            
        return line

    def readProtComp(self, pyStr,format):
# Header
        line = ""
        lineList = ["Genome-name1", "Genome-name2", "bit-score", "bbh-percent"]
    
    #   Check for valid formats
        if format not in ['tab','csv']:
            print ("Invalid format. Valid formats are tab and csv")
            return
        elif format == 'tab':
            line = "\t".join(lineList)
        elif format == 'csv':
            line = "'" + ",".join(lineList)
        # Add line-end to the header
        line += "\n"
        #print (pyStr)
        
        names1 = pyStr["proteome1names"]
        names2 = pyStr["proteome2names"]
        pairs1 = pyStr["data1"]
        foundPairs = {}
    
        count = 0
        for pos1, name1 in enumerate(names1):
            if not pairs1[pos1]:
                # The list is empty, no genome2 gene
                lineList = [name1, '--']
                if format == 'tab':
                    line += "\t".join(lineList) + "\n"
                elif format == 'csv':
                    line += ",".join(lineList) + "\n"
                continue
            for pair in pairs1[pos1]:
                pos2 = pair[0]
                loc = str(pos1) + ".." + str(pos2)
                foundPairs[loc] = 'Y'
                bit_score = pair[1]
                bbh_percent = pair[2]
                name2 = names2[pos2]
                lineList = [name1, name2, str(bit_score), str(bbh_percent)]
                if format == 'tab':
                    line += "\t".join(lineList) + "\n"
                elif format == 'csv':
                    line += ",".join(lineList) + "\n"
                count += 1
                
        pairs2 = pyStr["data2"]
        for pos2, name2 in enumerate(names2):
            if not pairs2[pos2]:
                # The list is empty, no genome1 gene
                lineList = ['--', name2]
                if format == 'tab':
                    line += "\t".join(lineList) + "\n"
                elif format == 'csv':
                    line += ",".join(lineList) + "\n"
            for pair in pairs2[pos2]:
                pos1 = pair[0]
                loc = str(pos1) + ".." + str(pos2)
                if loc in foundPairs:
                    continue
                bit_score = pair[1]
                bbh_percent = pair[2]
                name1 = names1[pos1]
                lineList = [name1, name2, str(bit_score), str(bbh_percent)]
                if format == 'tab':
                    line += "\t".join(lineList) + "\n"
                elif format == 'csv':
                    line += ",".join(lineList) + "\n"
                count += 1             
        #print ("DEBUG LINE: ", line)
        return line
