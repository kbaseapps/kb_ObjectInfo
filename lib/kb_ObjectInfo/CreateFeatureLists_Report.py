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
        #(self.cats, self.cat2name, self.cat2group, self.domfam2cat, self.domfam2name, self.domfam2ns) = self._configure_categories()
        logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')

    def configure_seed_categories(self):
        domain_desc_basepath = os.path.abspath('/kb/module/data')
        domfam2name  = dict()
        domfam2ns = dict()
        domfam2cat   = dict()
        cats = []
        cat2name     = dict()
        cat2group    = dict()
            
        # Initialize
        for namespace in ['SEED']:
            cat2name[namespace] = dict()
            cat2group[namespace] = dict()
            
        # get high-level cats
        with open(os.path.join(domain_desc_basepath, 'SEED_funcat.txt'), 'r') as dom_cat_handle:
            for line in dom_cat_handle.readlines():
                line = line.strip()
                [cat_group, cat] = line.split("\t")[0:2]
                if cat not in cats:
                    cats.append(cat)

        # get domfam to cat map, and vice versa
        with open(os.path.join(domain_desc_basepath, 'SEED_subsys.txt'), 'r') as dom2cat_map_handle:
            for line in dom2cat_map_handle.readlines():
                line = line.strip()
                [cat_group, cat_subgroup, cat, domfam] = line.split("\t")[0:4]
                domfam = domfam.strip()
                domfam = re.sub(' *\#.*$', '', domfam)
                domfam = re.sub(' *\(EC [\d\.\-\w]*\) *$', '', domfam)
                domfam = re.sub(' *\(TC [\d\.\-\w]*\) *$', '', domfam)
                domfam = re.sub('_', ' ', domfam)
                dom_name = domfam
                dom_name = cat_subgroup
                cat2group['SEED'][cat] = cat_group
                        
                domfam2ns[domfam] = 'SEED'
                domfam2cat[domfam] = cat
                domfam2name[domfam]  = dom_name
        
        
        return(cats, cat2name, cat2group, domfam2cat, domfam2name, domfam2ns)
        
    # -----------------------------------------------------------------
    #    Create a Delimited Table version of the genes in a genome
    #

    def delimitedTable(self, genome, features):

        (cats, cat2name, cat2group, domfam2cat, domfam2name, domfam2ns) = self.configure_seed_categories()
        seed_cat = domfam2cat
        
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
                feat['type'] = 'feature'

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
                if cat in cat2group['SEED']:
                    catgroup = cat2group['SEED'][cat]
                if domfam in domfam2name:
                    subgroup = domfam2name[domfam]
 
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
    def configure_domains(self):

        # https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CDSource_external
        # cd00009   - Conserved Protein Domain Family - curated at NCBI
        # cl21655   - Part of CDD at NCBI but not predicted
        # COG5463   - Clusters of Orthologous Genes (COGs)
        # KOG0989   - Eukaryotic COGs
        # NF011716  - NCBIfams
        # PF00015   - Protein Families - Pfam
        # PRK000000 - Protein Clusters - NCBI - okay to exclue CHL, MTH, PHA, PLN, PTZ, convert pfam to PF
        # sd000043   - Conserved Protein Domain Family - Domain models specifically built to annotate structural motifs;
        # smart00156 - Simple Modular Architecture Research Tool
        # TIGR03904  - The Institute for Genomic Research's database of protein families - TIGRfam
        
        domain_desc_basepath = os.path.abspath('/kb/module/data')
        domfam2name  = dict()
        domfam2ns = dict()
        domfam2cat   = dict()
        domfam2short = dict()
        cats = []
        cat2name     = dict()
        cat2group    = dict()
            
        # Initialize
        for namespace in ['COG','PF','TIGR','PRK','NCBIfams','cdd','Other','smart']:
            cat2name[namespace] = dict()
            cat2group[namespace] = dict()
#
#       Define all the categories and all of the domains
#       All the old and new data sources now combined in combine_data_sources.py
#
        with open(os.path.join(domain_desc_basepath,'all_domains.tsv'), 'r') as dom2cat_map_handle:
            for line in dom2cat_map_handle.readlines():
                [namespace,domfam,short_name,dom_name,cat] = line.split("\t")[0:5]
                domfam2ns[domfam] = namespace
                domfam2name[domfam]  = dom_name
                domfam2cat[domfam] = cat
                domfam2short[domfam] = short_name

        with open(os.path.join(domain_desc_basepath,'all_categories.tsv'), 'r') as dom2cat_map_handle:
            for line in dom2cat_map_handle.readlines():
                [namespace,cat,cat_name,cat_group] = line.split("\t")[0:4]
                if cat not in cats:
                    cats.append(cat)
                cat2name[namespace][cat] = cat_name
                cat2group[namespace][cat] = cat_group
        
        return(cats, cat2name, cat2group, domfam2cat, domfam2name, domfam2ns, domfam2short)
        
    #
    #   OBJECT: DomainAnnotation
    #   FORMAT: tab or comma delimited list of the genes, domains, e-values, and start/stop of domain hit
    #   Loop through all of the contigs and get all of the genes
    #   Uses printGeneDomain to print out individual lines
    #
    def readDomainAnnList(self, pyStr, cutoff):
    
        (cats, cat2name, cat2group, domfam2cat, domfam2name, domfam2ns, domfam2short) = self.configure_domains()
        # Header
        rpt_list1 = [["Contig", "Gene ID", "Domain", "Short Name", "Evalue", "Start", "Stop","Domain Name","Namespace", "Category", "Category Name", "Category Group"]]

        myData = pyStr['data']
        
        for contig in myData:
            contigData = myData[contig]
            for gene in contigData:
                if (gene[4]):
                    geneDomain = gene[4]
                    geneName   = gene[0]
        
                    for domain in geneDomain:
                        list = geneDomain[domain]
                        if list[0][2] < cutoff:
                            if '.' in domain and len(domain) < 15:
                                domain = domain.split('.')[0]

                            if domain in domfam2ns:
                                namespace = domfam2ns[domain]
                            else:
                                namespace = 'Other'
                                if domain.startswith('COG'):
                                    namespace = 'COG'
                                elif domain.startswith('PF') or domain.startswith('pfam'):
                                    namespace = 'PF'
                                elif domain.startswith('TIGR'):
                                    namespace = 'TIGR'
                                elif domain.startswith('PRK'):
                                    namespace = 'PRK'
                                elif domain.startswith('smart'):
                                    namespace = 'smart'
                                elif domain.startswith('NF'):
                                    namespace = 'NCBIfams'
                                elif domain.startswith('cd') or domain.startswith('cl') or domain.startswith('sd'):
                                    namespace = 'cdd'
                                domfam2ns[domain] = namespace

                            if domain in domfam2name:
                                dom_name = domfam2name[domain]
                            else:
                                dom_name = 'Other'
                                domfam2name[domain] = dom_name

                            if domain in domfam2cat:
                                cat = domfam2cat[domain]
                            else:
                                cat = 'Other'
                                domfam2cat[domain] = cat

                            if domain in domfam2short:
                                short_name = domfam2short[domain]
                            else:
                                short_name = 'Other'
                                domfam2short[domain] = cat

                            if namespace not in cat2name:
                                cat2name[namespace]  = dict()
                            if namespace not in cat2group:
                                cat2group[namespace] = dict()
                    
                            if cat > ' ' and namespace in cat2name and cat in cat2name[namespace]:
                                cat_name = cat2name[namespace][cat]
                            else:
                                cat = 'Other'
                                cat_name = 'Other'
                                cat2name[namespace][cat] = cat_name

                            if cat > ' '  and namespace in cat2group and cat in cat2group[namespace] :
                                cat_group = cat2group[namespace][cat]
                            else:
                                cat_group = 'Other'
                                cat2group[namespace][cat] = cat_group
                                        
                            rpt_list1.append([contig, geneName, domain, short_name, str(list[0][2]), str(list[0][0]), str(list[0][1]), dom_name, namespace, cat, cat_name, str(cat_group)])
                
        # Header
        rpt_list2 = [["Count","Domain","Short Name","Namespace","Category","Category Name","Category Group","Domain Name"]]

        myDict = {}
        
        for (contig, geneName, domain, short_name, evalue, start, stop, dom_name, namespace, cat, cat_name, cat_group) in rpt_list1:
            # First line from previous list
            if contig == 'Contig':
                continue
                
            if domain in myDict:
                myDict[domain] += 1
            else:
                myDict[domain] = 1
                
        for domain in sorted(myDict.keys()):
            namespace  = domfam2ns[domain]
            dom_name   = domfam2name[domain]
            short_name = domfam2short[domain]
            cat        = domfam2cat[domain]
            cat_name   = cat2name[namespace][cat]
            cat_group  = cat2group[namespace][cat]
            
            rpt_list2.append([str(myDict[domain]),domain,short_name,namespace,cat,cat_name,cat_group,dom_name])
                                        
        # Header
        rpt_list3 = [["Namespace","Category Group","Category","Category Name","Count"]]
        
        ns2cat         = {}
        ns2cat['COG']  = []
        ns2cat['PF']   = []
        ns2cat['TIGR'] = []
        myDict         = {}
        myDict['COG']  = {}
        myDict['PF']   = {}
        myDict['TIGR'] = {}
        
        for (geneCount,domain,short_name, namespace,cat,cat_name,cat_group,dom_name) in rpt_list2:
                
            # No category. No point in summarizing
            if namespace != 'COG' and namespace != 'PF' and namespace != 'TIGR':
                continue
        
            # First line
            if domain == 'Domain':
                continue
                
            if cat not in cat2group[namespace]:
                cat2group[namespace][cat] = 'Other'
            if cat not in cat2name[namespace]:
                cat2name[namespace][cat] = 'Other'
            
            if cat in myDict[namespace]:
                myDict[namespace][cat] += int(geneCount)
            else:
                myDict[namespace][cat] = int(geneCount)
                    
                if namespace in ns2cat:
                    ns2cat[namespace].append(cat)
                else:
                    ns2cat[namespace] = [cat]
        
        for namespace in sorted(ns2cat.keys()):
            if namespace == 'COG':
                cat_list = ['D','M','N','O','T','U','V','W','Y','Z','A', 'B','J','K','L', 'C', 'E', 'F', 'G', 'H','I','P', 'Q', 'X', 'R', 'S']
            else:
                cat_list = sorted(ns2cat[namespace])
                
            for cat in cat_list:
                if cat in myDict[namespace]:
                    rpt_list3.append([namespace,cat2group[namespace][cat],cat,cat2name[namespace][cat],str(myDict[namespace][cat])])
      
        return (rpt_list1, rpt_list2, rpt_list3)

    #
    #   OBJECT: FeatureSet or SequenceSet
    #   FORMAT: List of the contents of the object
    #
    def readFeatSeq(self, pyStr):
        
        # Header
        desc_list = []
        rpt_list = []
        cf = CreateFasta(self.config)
        seq_list = []
        header = ""
#
#   Type 1 - Order matters
#
        if 'description' in pyStr and 'elements' in pyStr and 'element_ordering' in pyStr:
            header = "Ordered Elements for "+str(pyStr['description'])
            desc_list = [["Genome ID","Genome Name"]]
            genome_names = {}
                        
            eleOrder = pyStr['element_ordering']
            for index in eleOrder:
                gid = pyStr['elements'][index][0]
                genome_name = self.dfu.get_objects({'object_refs': [gid]})['data'][0]['info'][1]
                genome_names[gid] = genome_name
                desc_list.append([pyStr['elements'][index][0], genome_name ])

            rpt_list = [["Index","Feature ID","Source Genome Object ID","Source Genome Name"]]
            
            count = 1
            for index in eleOrder:
                gid = pyStr['elements'][index][0]
                rpt_list.append([str(count), index, gid, genome_names[gid]])
                count += 1

#
#   Type 2 - Unordered
#
        elif 'description' in pyStr and 'elements' in pyStr:
            header = "Unordered Elements for "+str(pyStr['description'])
            desc_list = [['Description:', pyStr['description']]]
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

            rpt_list += [["Feature ID","Source Genome Object ID","Genome"]]
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
            header = "Sequences"
            desc_list = [['Set Description', pyStr['description']],["Sequence Set ID", pyStr['sequence_set_id']]]
            mySequences = pyStr['sequences']
            count = 0
            for seq in mySequences:
                seqline = cf.splitSequence(seq['sequence'])
                seq_list.append([">" + seq['sequence_id']+"  "+seq['description']])
                seq_list.extend(seqline)
                count += 1

#
#   Type Unknown
#
        else:
            logging.error("This type of FeatureSet has not been described yet")

        return (header, desc_list, rpt_list, seq_list)

    def readProtComp(self, pyStr):
        id1 = pyStr["genome1ref"]
        genome1 = self.dfu.get_objects({'object_refs': [id1]})['data'][0]['info'][1]
        id2 = pyStr["genome2ref"]
        genome2 = self.dfu.get_objects({'object_refs': [id2]})['data'][0]['info'][1]
        
        rpt_list1 = [["Genome1", genome1],["Genome2",genome2]]
        rpt_list2 = [["Genome1", "Genome2", "bit-score", "bbh-percent"]]

        names1 = pyStr["proteome1names"]
        names2 = pyStr["proteome2names"]
        pairs1 = pyStr["data1"]
        foundPairs = {}
    
        count = 0
        for pos1, name1 in enumerate(names1):
            if not pairs1[pos1]:
                # The list is empty, no genome2 gene
                rpt_list2.append([name1, 'NA', '0', '0'])
                continue
            for pair in pairs1[pos1]:
                pos2 = pair[0]
                loc = str(pos1) + ".." + str(pos2)
                foundPairs[loc] = 'Y'
                bit_score = pair[1]
                bbh_percent = pair[2]
                name2 = names2[pos2]
                rpt_list2.append([name1, name2, str(bit_score), str(bbh_percent)])
                
                count += 1
                
        pairs2 = pyStr["data2"]
        for pos2, name2 in enumerate(names2):
            if not pairs2[pos2]:
                # The list is empty, no genome1 gene
                rpt_list2.append(['NA', name2, '0', '0'])
                
            for pair in pairs2[pos2]:
                pos1 = pair[0]
                loc = str(pos1) + ".." + str(pos2)
                if loc in foundPairs:
                    continue
                bit_score = pair[1]
                bbh_percent = pair[2]
                name1 = names1[pos1]
                rpt_list2.append([name1, name2, str(bit_score), str(bbh_percent)])

                count += 1             

        return rpt_list1, rpt_list2
