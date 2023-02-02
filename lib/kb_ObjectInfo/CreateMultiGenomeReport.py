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
    def getGenomeSet(self, obj_name, obj_id, obj_data, format, rttype):

        lst = ['unk', 'unk', 'unk', 'unk', 'unk', 'unk', 'unk', 'unk', 'unk']
        domain, size, num_feat, gc_cont, num_ctg, source, gen_code, assembly, sci_name = lst
        features = { 'gene' : 0, 'CDS' : 0, 'rRNA' : 0, 'tRNA' : 0, 'other' : 0}
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
        if 'scientific_name' in obj_data['data']:
            sci_name = obj_data['data']['scientific_name']
        if 'assembly_ref' in obj_data['info'][10]:
            assembly = obj_data['info'][10]
        elif 'assembly_ref' in obj_data['data']:
            assembly = obj_data['data']['assembly_ref']
        elif 'contigset_ref' in obj_data['data']:
            assembly = obj_data['data']['contigset_ref']
        if 'features' in obj_data['data']:
            for feat in obj_data['data']['features']:
                if 'type' in feat:
                    type = feat['type']
                    if type in features:
                        features[type] += 1
                    else:
                        features['other'] += 1
        line = ''
        rpt_list = []

        if format == 'list':
            line = "Description for: " + obj_name + "\n"

            line += name + "\n"
            line += "\tObject ID:     {0:s}\n\tScientific Name:   {1:s}\n\tSize:         {2}\n\tSource:       {3:s}\n\tDomain:       {4:s}\n\tAssembly Reference: {5:s}\n".format(
                obj_id, sci_name, size, source, domain, assembly)
            line += "\tFeatures:     {0}\n\tContigs:      {1}\n\tPercent GC:      {2}\n\tGenetic Code: {3}\n".format(num_feat, num_ctg, gc_cont, gen_code)

            rpt_list = [["Description for:",obj_name],["Object ID:",obj_id],["Scientific Name:",sci_name],
                       ["Size:",size],["Source:",source],["Domain:",domain],["Assembly Reference:",assembly],
                       ["Features:",num_feat],["Contigs:",num_ctg],["Percent GC:",gc_cont],["Genetic Code:",gen_code]]

            for feat in sorted(features):
                line += "\t{:8s}      {}\n".format(feat, features[feat])
                rpt_list.append([feat,str(features[feat])])

        elif format == 'tab':
            lst = ["Name", "Genome Object ID", "Scientific Name", "Size", "Source", "Domain", "Assembly Object ID", "Features", "Contigs", "Percent GC",
                   "Genetic Code", "Number of CDS", "Number of gene", "Number of other", "Number of rRNA", "Number of tRNA"]
            line = "\t".join(lst) + "\n"
            rpt_list = [["Name", "Genome Object ID", "Scientific Name", "Size", "Source", "Domain", "Assembly Object ID", "Features",
                        "Contigs", "Percent GC", "Genetic Code", "Number of CDS", "Number of gene", "Number of other", "Number of rRNA",
                        "Number of tRNA"]]
                   
            lst = [name, obj_id, sci_name, size, source, domain, assembly, num_feat, num_ctg, gc_cont, gen_code]
            line = "\t".join(lst)
            
            for feat in sorted(features):
                line += "\t" + str(features[feat])
                lst.append(str(features[feat]))

            line += "\n"
            rpt_list.append(lst)

        elif format == 'csv':
            lst = ["Name", "Genome Object ID", "Scientific Name", "Size", "Source", "Domain", "Assembly Object ID", "Features", "Contigs", "Percent GC",
                   "Genetic Code", "Number of CDS", "Number of gene", "Number of other", "Number of rRNA", "Number of tRNA"]
            line = ",".join(lst) + "\n"
            rpt_list = [["Name", "Genome Object ID", "Scientific Name", "Size", "Source", "Domain", "Assembly Object ID", "Features",
                        "Contigs", "Percent GC", "Genetic Code", "Number of CDS", "Number of gene", "Number of other", "Number of rRNA",
                        "Number of tRNA"]]
                        
            lst = [name, obj_id, sci_name, size, source, domain, assembly, num_feat, num_ctg, gc_cont, gen_code]
            line = ",".join(lst)
            for feat in sorted(features):
                line += "," + str(features[feat])
                lst.append(str(features[feat]))
                
            line += "\n"
            rpt_list.append(lst)
            
        if rttype == 'line':
            return line
        else:
            return rpt_list

    # Metadata for a GenomeSet
    #
    def getGenomeSetMeta(self, obj_data):
        line = ''
        line += "Name         {}\n".format(obj_data['info'][1])
        line += "Type         {}\n".format(obj_data['info'][2])
        line += "Created By   {}\n".format(obj_data['info'][5])
        line += "Narrative    {}\n".format(obj_data['info'][7])
        line += "Description  {}\n".format(obj_data['data']['description'])
        line += "Number of Elements {}\n".format(str(len(obj_data['data']['elements'])))
        rpt_list = []
        rpt_list = [["Name",obj_data['info'][1]],["Type",obj_data['info'][2]],["Created By",obj_data['info'][5]],
                    ["Narrative",obj_data['info'][7]],["Description",obj_data['data']['description']],
                    ["Number of Elements",str(len(obj_data['data']['elements']))]]
        
        for ele in obj_data['data']['elements']:
            gref = obj_data['data']['elements'][ele]['ref']
            genome = self.dfu.get_objects({'object_refs': [gref]})
            name = genome['data'][0]['info'][1]
            sci_name = genome['data'][0]['data']['scientific_name']
            
            line += "  Element:   {0:s}\t{1:s}\t{2:s}\n".format(ele, name, sci_name)
            rpt_list += [["Element:",ele, name, sci_name]]
        return line, rpt_list

    # Describe a GenomeSet
    #
    def readGenomeSet(self, obj_name, pyStr, format):
        myGS = pyStr['elements']
        line = ''
        rpt_list = []
        
        for ele in myGS:
            genome = self.dfu.get_objects({'object_refs': [myGS[ele]['ref']]})
            line += self.getGenomeSet(obj_name,myGS[ele]['ref'], genome['data'][0], format, 'line')
            rpt_list.extend(self.getGenomeSet(obj_name,myGS[ele]['ref'], genome['data'][0], format,'rpt'))

        return line, rpt_list

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
