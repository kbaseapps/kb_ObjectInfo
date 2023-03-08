import time
import os
import re
import logging

from .CreateFasta_Report import CreateFasta
from installed_clients.DataFileUtilClient import DataFileUtil
from pprint import pprint, pformat

class CreateFeatureLists:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        (self.cats, self.cat2name, self.cat2group, self.domfam2cat, self.domfam2name, self.domfam2ns) = self._configure_categories()
        logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')

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

    def delimitedTable(self, genome, features):

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
        
        rpt_list = [["Feature ID", "Feature type", "Contig", "Start", "Stop", "Strand", "Feature Function", "Aliases",
                    "RAST Functional Assignment", "RAST Functional Group 1", "RAST Functional Group 2"]]
                    
        cat = ''
        domfam = ''
        for feat in genome[features]:
            if 'function' not in feat:
                feat['function'] = 'unknown'
            else:
                if feat['function'] in seed_cat:
                    domfam = feat['function']
                    cat = seed_cat[domfam]

            if 'functions' in feat:
                feat['function'] = ', '.join(feat['functions'])
                for func in feat['functions']:
                    if func in seed_cat:
                        domfam = func
                        cat = seed_cat[domfam]
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
 
            rpt_list.append([feat['id'], feat['type'], contig, str(start), str(stop), strand, feat['function'], aliases, cat, subgroup, catgroup])

        return rpt_list


        # -----------------------------------------------------------------
        #    Create a GFF3 version of the features in a genome
        #

    def gff3(self, genome, features):
        rpt_list = []
        if features not in genome:
            return rpt_list

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
            rpt_list.append([contig, ph, feat['type'], str(start), str(stop), ph, strand, ph, attrib])

        return rpt_list



    # -----------------------------------------------------------------
    #   Domain Annotation Reports
    #

    #   OBJECT: DomainAnnotation
    #   FUNCTION: User-defined function to format all the domains for a gene
    #
    def printGeneDomain(self, contig, geneName, geneDomain, cutoff):
        rpt_list = []
        
        for domain in geneDomain:
            list = geneDomain[domain]
            if list[0][2] < cutoff:
                if '.' in domain and len(domain) < 15:
                    domain = domain.split('.')[0]

                if domain in self.domfam2ns:
                    namespace = self.domfam2ns[domain]
                else:
                    namespace = ' '

                if domain in self.domfam2name:
                    dom_name = self.domfam2name[domain]
                else:
                    dom_name = ' '

                if domain in self.domfam2cat:
                    cat = self.domfam2cat[domain]
                else:
                    cat = ' '

                if cat > ' ' and cat in self.cat2name[namespace]:
                    cat_name = self.cat2name[namespace][cat]
                else:
                    cat_name = ' '

                if cat > ' ' and cat in self.cat2group[namespace]:
                    cat_group = self.cat2group[namespace][cat]
                else:
                    cat_group = ' '

                rpt_list.append([contig, geneName, domain, str(list[0][2]), str(list[0][0]), str(list[0][1]), dom_name, cat, cat_name, str(cat_group)])
                
        # Returning a list of lists
        return rpt_list

    #
    #   OBJECT: DomainAnnotation
    #   FORMAT: tab or comma delimited list of the genes, domains, e-values, and start/stop of domain hit
    #   Loop through all of the contigs and get all of the genes
    #   Uses printGeneDomain to print out individual lines
    #
    def readDomainAnnList(self, pyStr, cutoff):
        # Header
        rpt_list = [["Contig", "Gene ID", "Domain", "Evalue", "Start", "Stop","Domain Name", "Category", "Category Name", "Category Group"]]

        myData = pyStr['data']

        for contig in myData:
            contigData = myData[contig]
            for gene in contigData:
                if (gene[4]):
                    domain = gene[4]

                    # The return list is a list of lists and this is done many times
                    # Need to append them to our list one at a time, otherwise list of lists of lists
                    rtn_list = (self.printGeneDomain(contig, gene[0], domain, cutoff))
                    for rtn in rtn_list:
                        rpt_list.append(rtn)
  
        return rpt_list

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
    def readDomainAnnCount(self, pyStr, cutoff):

        # Header
        rpt_list = [["Domain", "Count"]]

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
            rpt_list.append([domain, str(myDict[domain])])
    
        return rpt_list

    #
    #   OBJECT: FeatureSet or SequenceSet
    #   FORMAT: List of the contents of the object
    #
    def readFeatSeq(self, pyStr):
        
        # Header
        rpt_list = []
        cf = CreateFasta(self.config)
#
#   Type 1 - Order matters
#
        if 'description' in pyStr and 'elements' in pyStr and 'element_ordering' in pyStr:
            rpt_list = [['Description', str(pyStr['description'])],['Genomes']]
            
            eleOrder = pyStr['element_ordering']
            for index in eleOrder:
                genome_name = self.dfu.get_objects({'object_refs': [pyStr['elements'][index][0]]})['data'][0]['info'][1]
                rpt_list.append([pyStr['elements'][index][0], genome_name ])

            rpt_list += ([" "],["Ordered Elements:"],["Index","Feature ID","Source Genome Object ID"])
            
            count = 1
            for index in eleOrder:
                rpt_list += ([[str(count)] + [index] + pyStr['elements'][index] ])
                count += 1

#
#   Type 2 - Unordered
#
        elif 'description' in pyStr and 'elements' in pyStr:
            rpt_list = [['Description:', pyStr['description']]]
            myElements = pyStr['elements']
            genome_names = {}
            
            for element in myElements:
                if isinstance(myElements[element],list):
                    for gid in myElements[element]:
                        if gid not in genome_names.keys():
                            genome_name = self.dfu.get_objects({'object_refs': [gid]})['data'][0]['info'][1]
                            genome_names[gid] = genome_name
                else:
                    gid = myElements[element]
                    if gid not in genome_names.keys():
                        genome_name = self.dfu.get_objects({'object_refs': [gid]})['data'][0]['info'][1]
                        genome_names[gid] = genome_name

            rpt_list += [[" "],["Unordered Elements:"],["Feature ID","Source Genome Object ID","Genome"]]
            count = 0
            for element in myElements:
                if isinstance(myElements[element],list):
                    for i in myElements[element]:
                        rpt_list.append([element, i, genome_names[i]])
                else:
                    rpt_list.append([element, myElements[element], genome_names[myElements[element]]])
                count += 1

#
#   Type 3 - With Sequences
#
        elif 'description' in pyStr and 'sequences' in pyStr and 'sequence_set_id' in pyStr:

            rpt_list = [['Set Description', pyStr['description']],["Sequence Set ID", pyStr['sequence_set_id']],["Sequences:"]]
            mySequences = pyStr['sequences']
            count = 0
            for seq in mySequences:
                seqline = cf.splitSequence(seq['sequence'])
                rpt_list.append([">" + seq['sequence_id']+'   '+seq['description']])
                rpt_list.extend(seqline)
                count += 1

#
#   Type Unknown
#
        else:
            logging.error("This type of FeatureSet has not been described yet")

        return rpt_list

    def readProtComp(self, pyStr):
# Header
 
        id1 = pyStr["genome1ref"]
        genome1 = self.dfu.get_objects({'object_refs': [id1]})['data'][0]['info'][1]
        id2 = pyStr["genome2ref"]
        genome2 = self.dfu.get_objects({'object_refs': [id2]})['data'][0]['info'][1]
        
        rpt_list = [["Genome1 = "+genome1],["Genome2 = "+genome2],[" "],
                    ["Genome1", "Genome2", "bit-score", "bbh-percent"]]

        names1 = pyStr["proteome1names"]
        names2 = pyStr["proteome2names"]
        pairs1 = pyStr["data1"]
        foundPairs = {}
    
        count = 0
        for pos1, name1 in enumerate(names1):
            if not pairs1[pos1]:
                # The list is empty, no genome2 gene
                rpt_list.append([name1, 'NA'])
                continue
            for pair in pairs1[pos1]:
                pos2 = pair[0]
                loc = str(pos1) + ".." + str(pos2)
                foundPairs[loc] = 'Y'
                bit_score = pair[1]
                bbh_percent = pair[2]
                name2 = names2[pos2]
                rpt_list.append([name1, name2, str(bit_score), str(bbh_percent)])
                
                count += 1
                
        pairs2 = pyStr["data2"]
        for pos2, name2 in enumerate(names2):
            if not pairs2[pos2]:
                # The list is empty, no genome1 gene
                rpt_list.append(['NA', name2])
                
            for pair in pairs2[pos2]:
                pos1 = pair[0]
                loc = str(pos1) + ".." + str(pos2)
                if loc in foundPairs:
                    continue
                bit_score = pair[1]
                bbh_percent = pair[2]
                name1 = names1[pos1]
                rpt_list.append([name1, name2, str(bit_score), str(bbh_percent)])

                count += 1             

        return rpt_list
