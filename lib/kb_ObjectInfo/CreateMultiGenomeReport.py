import time
import os

from installed_clients.DataFileUtilClient import DataFileUtil

class CreateMultiGenomeReport:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)

    # Listing of the Elements of a GenomeSet
    #
    def getGenomeSet(self, obj_name, obj_id, obj_data, format):
        lst = ['unk', 'unk', 'unk', 'unk', 'unk', 'unk', 'unk', 'unk', 'unk']
        domain, size, num_feat, gc_cont, num_ctg, source, gen_code, assembly, sci_name = lst
        features = { 'CDS' : 0, 'gene' : 0, 'non_coding_features': 0, 'other' : 0, 'rRNA' : 0, 'tRNA' : 0 }

#       Metadata available in obj_data['info']
        name = obj_data['info'][1]
        if 'Domain' in obj_data['info'][10]:
            domain = obj_data['info'][10]['Domain']
        if 'Size' in obj_data['info'][10]:
            size = obj_data['info'][10]['Size']
        if 'Number features' in obj_data['info'][10]:
            num_feat = obj_data['info'][10]['Number features']
        if 'GC content' in obj_data['info'][10]:
            gc_cont = str(float(obj_data['info'][10]['GC content']) * 100)
        if 'Number contigs' in obj_data['info'][10]:
            num_ctg = obj_data['info'][10]['Number contigs']
        if 'Source' in obj_data['info'][10]:
            source = obj_data['info'][10]['Source']
        if 'Genetic code' in obj_data['info'][10]:
            gen_code = obj_data['info'][10]['Genetic code']
            
#       Metadata available in obj_data['data']
        if 'scientific_name' in obj_data['data']:
            sci_name = obj_data['data']['scientific_name']
        if 'assembly_ref' in obj_data['info'][10]:
            assembly = obj_data['info'][10]
        elif 'assembly_ref' in obj_data['data']:
            assembly = obj_data['data']['assembly_ref']
        elif 'contigset_ref' in obj_data['data']:
            assembly = obj_data['data']['contigset_ref']
            
#       Genomes v.9 and newer have feature_counts added up already
        if 'feature_counts' in obj_data['data']:
            for feat in obj_data['data']['feature_counts']:
                if feat in features:
                    features[feat] = obj_data['data']['feature_counts'][feat]
                else:
                    features['other'] += 1
                    
#       Older genomes don't have the feature_counts so you have to add them up by type
        elif 'features' in obj_data['data']:
            for feat in obj_data['data']['features']:
                if 'type' in feat:
                    feat_type = feat['type']
                    if feat_type in features:
                        features[feat_type] += 1
                    else:
                        features['other'] += 1

        rpt_list = []

#       For the list format, give the headers and values in two columns
        if format == 'list':
            rpt_list = [["Description for:",obj_name],["Object ID:",obj_id],["Scientific Name:",sci_name], ["Size:",size],["Source:",source],["Domain:",domain],["Assembly Reference:",assembly], ["Features:",num_feat],["Contigs:",num_ctg],["Percent GC:",gc_cont],["Genetic Code:",gen_code]]
            for feat in sorted(features):
                rpt_list.append([feat, str(features[feat])])
                
#       For tab and csv, give all the values as columns
        elif format == 'tab' or format == 'csv':
            rpt_list = [name, obj_id, sci_name, size, source, domain, assembly, num_feat, num_ctg, gc_cont, gen_code]
#           There are extra features that need to be displayed
            for feat in sorted(features):
                rpt_list.append(str(features[feat]))

#       Return a list
        return rpt_list

    # Describe a GenomeSet
    #
    def readGenomeSet(self, obj_name, pyStr, format):
        myGS = pyStr['elements']
        rpt_list = []

        if format == 'tab' or format == 'csv':
            rpt_list = [["Name", "Genome Object ID", "Scientific Name", "Size", "Source", "Domain", "Assembly Object ID", "Features",
                        "Contigs", "Percent GC", "Genetic Code", "Number of CDS", "Number of gene", "Number of Non-coding",
                        "Number of other", "Number of rRNA", "Number of tRNA"]]
 
        for ele in myGS:
            genome = self.dfu.get_objects({'object_refs': [myGS[ele]['ref']]})
            if format == 'tab' or format == 'csv':
                rpt_list.extend([self.getGenomeSet(obj_name,myGS[ele]['ref'], genome['data'][0], format)])
            else:
                rpt_list.extend(self.getGenomeSet(obj_name,myGS[ele]['ref'], genome['data'][0], format))
            
        return rpt_list

    # Metadata for a GenomeSet
    #
    def getGenomeSetMeta(self, obj_data):
        desc = 'None'
        if obj_data['data']['description']:
            desc = obj_data['data']['description']
            
        rpt_list = []
        rpt_list = [["Name",obj_data['info'][1]],["Type",obj_data['info'][2]],["Created By",obj_data['info'][5]],
                    ["Narrative",obj_data['info'][7]],["Description",desc],
                    ["Number of Elements",str(len(obj_data['data']['elements']))]]
        
        for ele in obj_data['data']['elements']:
            gref = obj_data['data']['elements'][ele]['ref']
            genome = self.dfu.get_objects({'object_refs': [gref]})
            name = genome['data'][0]['info'][1]
            sci_name = genome['data'][0]['data']['scientific_name']
            rpt_list.append(["Genome:",ele, name, sci_name])
            
        return rpt_list

    # Return the assembly_refs
    #
    def getAssemblyRef(self, assem):
        myGS = assem['data']['elements']
        assembly_list = []
        for ele in myGS:
            genome = self.dfu.get_objects({'object_refs': [myGS[ele]['ref']]})
            obj_data = genome['data'][0]
            sci_name = obj_data['data']['scientific_name']

            assembly = ""
            if 'assembly_ref' in obj_data['info'][10]:
                assembly = obj_data['info'][10]
            elif 'assembly_ref' in obj_data['data']:
                assembly = obj_data['data']['assembly_ref']
            elif 'contigset_ref' in obj_data['data']:
                assembly = obj_data['data']['contigset_ref']
            assembly_list.append(assembly+':'+sci_name)

        return assembly_list
