#####################
#
#   Combine all of the data sources.
#   While there is much overlap,each source has one or more unique records.
#   This process makes the import of the definitions cleaner in scripts.
#   Unnecessary data has been eliminated
#
#   Define the domains (namespace, ID, short name, and description/name)
#   Where possible, also define a category (category and category description)
#   In some cases, the categories also have a higher level grouping
#
# NUMBERING SYSTEMS USED IN THE DOMAIN ANNOTATION
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
#####################

import re

domfam2ns = dict()
domfam2short = dict()
domfam2name  = dict()
domfam2cat   = dict()

cats = []
cat2name     = dict()
cat2group    = dict()


path = '/Users/mland/my_KBASE/kb_ObjectInfo/data/'
###############
#  Get the COG functional categories
###############
namespace = 'COG'
cat_names_path      = path + 'COG_2014_funcat.tsv'
if namespace not in cat2name:
    cat2name[namespace]  = dict()
if namespace not in cat2group:
    cat2group[namespace] = dict()
    
with open(cat_names_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        [cat, cat_group, cat_name] = newline.split("\t")[0:3]
        if  cat not in cats:
            cats.append(cat)
        
        #print(cat,cat_name,cat_group)
        cat2name[namespace][cat] = cat_name
        cat2group[namespace][cat] = cat_group
file1.close()

###############
# Get the COG domain names 2014
###############
namespace = 'COG'
domain_to_cat_path  = path + 'COG_2014.tsv'
count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        (domain,cat_str,dom_name) = newline.split("\t")
        cat = cat_str[0]  # only use first cat
    
        #print("A",namespace,cat,domain,dom_name)
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2cat[domain]  = cat
        domfam2short[domain] = 'Other'
        if  cat not in cats:
            cats.append(cat)
        count += 1
file1.close()
print ("FINAL COG 2014 NUMBER",count)

###############
# Get the COG domain updates
###############
domain_to_cat_path  = path + 'COG-20.def.tsv'
count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        # First line header
        if "COG ID" in line:
            continue
        (domain,cat,dom_name,short_name) = line.split("\t")[0:4]
        dom_name = dom_name.replace('"','')
        
        if domain in domfam2ns:
            if short_name > '     ':
                domfam2short[domain] = short_name
            continue
                    
        if short_name < '     ':
            short_name = 'Other'
        if dom_name < '     ':
            dom_name = 'Other'
            
        #print("B",namespace,cat,domain,short_name,dom_name)
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2cat[domain]  = cat
        domfam2short[domain] = short_name
        count += 1
file1.close()
print ("FINAL COG v20 NUMBER",count)

###############
#  Get the Pfam functional categories
###############
namespace = 'PF'
cat_group = 'NA'
cat_names_path     = path + 'Pfam-C.clans_names.v35.tsv'
if namespace not in cat2name:
    cat2name[namespace]  = dict()
if namespace not in cat2group:
    cat2group[namespace] = dict()
    
with open(cat_names_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        [cat, cat_id, cat_name] = newline.split("\t")[0:3]
        if cat < '      ':
            cat = 'Other'
        if cat not in cats:
            cats.append(cat)
        cat2name[namespace][cat] = cat_name
        cat2group[namespace][cat] = 'NA'
            
        #print(cat,cat_name,cat_group)
file1.close()

###############
# Get the PF domain names
###############
namespace = 'PF'
cat2name[namespace]['Other'] = 'Other'
cat2group[namespace]['Other'] = 'NA'
domain_to_cat_path = path + 'Pfam-A.clans.tsv'
count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        (domain, cat, cat_name, short_name, dom_name) = newline.split("\t")[0:5]
        
        if cat < '      ':
            cat = 'Other'
            cat_name = 'Other'
        if  cat not in cats:
            cats.append(cat)
            cat2name[namespace][cat] = cat_name
            cat2group[namespace][cat] = 'NA'
            
        #print("A",namespace,cat,domain,dom_name)
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2short[domain] = short_name
        domfam2cat[domain]  = cat
        count += 1
file1.close()
print ("FINAL PfamA.orig NUMBER",count)

###############
# Get the PF domain name updates
###############
namespace = 'PF'
domain_to_cat_path = path + 'Pfam-A.clans.v35.tsv'
cat = 'Other'
count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        (domain, cat, cat_name, short_name, dom_name) = newline.split("\t")[0:5]

        if domain in domfam2name:
            if domfam2cat[domain] == 'Other' and cat > '     ':
                domfam2cat[domain] = cat
            if domfam2name[domain] < '    ' and dom_name > '     ':
                domfam2name[domain] = dom_name
            if domfam2short[domain] < '    ' and short_name > '     ':
                domfam2short[domain] = short_name

        if cat < '      ':
            cat = 'Other'
            cat_name = 'Other'
        if  cat not in cats:
            cats.append(cat)
            cat2name[namespace][cat] = cat_name
            cat2group[namespace][cat] = 'NA'
            
        #print("A",namespace,cat,domain,dom_name)
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2short[domain] = short_name
        domfam2cat[domain]  = cat
        count += 1
file1.close()
print ("FINAL PfamA.v35 NUMBER",count)

###############
#  Get the TIGRfam functional categories
###############
namespace = 'TIGR'
cat_group = 'NA'
cat_names_path     = path + 'TIGRrole2go.txt'
tigrrole_id2cat = dict()
if namespace not in cat2name:
    cat2name[namespace]  = dict()
if namespace not in cat2group:
    cat2group[namespace] = dict()
    
with open(cat_names_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
                            
        if line.startswith('!'):
            continue
        cat_id = '0'
        [cat, cat_id, cat_group, cat_name_plus_go_terms] = newline.split("\t")[0:4]
        if cat_id != '' and int(cat_id) != 0:
            tigrrole_id2cat[cat_id] = cat
        if cat not in cats:
            cats.append(cat)
        cat_name = re.sub(' *\> GO:.*$', '', cat_name_plus_go_terms)
        cat2name[namespace][cat] = cat_name
        cat2group[namespace][cat] = cat_group
            
        #print(cat,cat_name,cat_group)
file1.close()

###############
# Get the TIGRfam domain names
###############
namespace = 'TIGR'
domain_to_cat_path = path + 'TIGRInfo.tsv'
if namespace not in cat2name:
    cat2name[namespace]  = dict()
if namespace not in cat2group:
    cat2group[namespace] = dict()
    count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        if line.startswith('!'):
            continue
        cat_id = '0'
        [domfam_id, domain, cat_group, cat_id, short_name, ec_id, dom_name] = newline.split("\t")[0:7]
        if cat_id != '' and int(cat_id) != 0 and cat_id in tigrrole_id2cat:
            cat = tigrrole_id2cat[cat_id]
        else:
            cat = 'Other'
            cat_name = 'Other'
            
        #print("A",namespace,cat,domain,dom_name)
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2short[domain] = short_name
        domfam2cat[domain]  = cat
        count += 1
file1.close()
print ("FINAL TIGR NUMBER",count)

###############
# Get the NCBIfams domain names
###############

namespace = 'NCBIfams'
cat = 'Other'
line_split = []
domain_to_cat_path = path + "hmm_PGAP.tsv"
count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        line = line.strip()
        if line.startswith('#'):
            continue
        line_split = line.split("\t")
        domain = line_split[1]
        short_name = line_split[2]
        dom_name  = line_split[10]
        
        if '.' in domain and len(domain) < 15:
            domain = domain.split('.')[0]

        if domain in domfam2name:
            continue

        if domain.startswith('NF'):
            namespace = 'NCBIfams'
        elif domain.startswith('PF'):
            namespace = 'PF'
        elif domain.startswith('TIGR'):
            namespace = 'TIGR'
        elif domain.startswith('PRK'):
            namespace = 'PRK'
        else:
            print(namespace,domain,short_name,dom_name)
        count += 1
                    
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2short[domain] = short_name
        domfam2cat[domain]  = 'Other'
file1.close()

print ("FINAL hmm NUMBER",count)

###############
# Get the cdd domain names
###############

domain_to_cat_path = path + 'cddid_all.tsv'
count = 0
with open(domain_to_cat_path, "r") as file1:
    for line in file1:
        line = line.strip()
        if line.startswith('!'):
            continue
        (pssmID,domain,short_name,dom_name) = line.split("\t")[0:4]
        domain.replace("pfam","PF")
        if domain in domfam2name:
            continue
        count += 1
        
        namespace = 'Other'
        cat = 'Other'
        if short_name < '     ':
            short_name = 'Other'
        if dom_name < '     ':
            dom_name = 'Other'
            
        if domain.startswith('CHL') or domain.startswith('MTH') or domain.startswith('PHA') or domain.startswith('PLN') or domain.startswith('PTZ'):
            # NCBI protein clusters
            continue
        elif domain.startswith('KOG'):
            # Eukaryotic COG
            continue
        elif domain.startswith('COG'):
            namespace = 'COG'
            x = dom_name.find('[')
            y = dom_name.find(']')
            if x > 0 and y > 0:
                cat = dom_name[x+1:y]
                dom_name = dom_name[0:x].replace('"','')
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
        else:
            print (line)
            
        if len(dom_name) > 100:
            dom_name = dom_name[0:100] + '...'
                        
        dom_name = dom_name.replace('"','')
        #print(namespace,cat,domain,short_name,dom_name)
            
        domfam2ns[domain]   = namespace
        domfam2name[domain] = dom_name
        domfam2short[domain] = short_name
        domfam2cat[domain]  = cat
file1.close()
print ("FINAL cdd NUMBER",count)

##############
# Save the Categories to a file
##############

domain_to_cat_path  = path + 'all_categories.tsv'
with open(domain_to_cat_path, "w") as file1:
    for namespace in ['COG','PF','TIGR']:
        for cat in cat2name[namespace]:
            outline = "\t".join([namespace,cat,cat2name[namespace][cat],cat2group[namespace][cat],"\n"])
            file1.write(outline)
file1.close()


##############
# Save the Domains to a file
##############

domain_to_cat_path  = path + 'all_domains.tsv'
with open(domain_to_cat_path, "w") as file1:
    for domain in domfam2ns:
        outline = "\t".join([domfam2ns[domain],domain,domfam2short[domain],domfam2name[domain],domfam2cat[domain],"\n"])
        file1.write(outline)
file1.close()



