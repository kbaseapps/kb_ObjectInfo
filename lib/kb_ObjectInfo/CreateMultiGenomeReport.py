import time
import os
from installed_clients.DataFileUtilClient import DataFileUtil

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class CreateMultiGenomeReport:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)

    # Listing of the Elements of a GenomeSet
    #
    def getGenomeSet(self, obj_id, obj_data, format):

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
        if format == 'list':
            line = name + "\n"
            line += "\tObjectID:     {0:s}\n\tScientName:   {1:s}\n\tSize:         {2}\n\tSource:       {3:s}\n\tDomain:       {4:s}\n\tAssembly Ref: {5:s}\n".format(
                obj_id, sci_name, size, source, domain, assembly)
            line += "\tFeatures:     {0}\n\tContigs:      {1}\n\tPct. GC:      {2}\n\tGenetic Code: {3}\n".format(num_feat, num_ctg, gc_cont, gen_code)

            for feat in sorted(features):
                line += "\t{:8s}      {}\n".format(feat, features[feat])

        if format == 'tab':
            lst = [name, obj_id, sci_name, size, source, domain, assembly, num_feat, num_ctg, gc_cont, gen_code]
            line = "\t".join(lst)
            for feat in sorted(features):
                line += "\t" + str(features[feat])
            line += "\n"
        if format == 'csv':
            lst = [name, obj_id, sci_name, size, source, domain, assembly, num_feat, num_ctg, gc_cont, gen_code]
            line = ",".join(lst)
            for feat in sorted(features):
                line += "," + str(features[feat])
            line += "\n"
        return line


    # Metadata for a GenomeSet
    #
    def getGenomeSetMeta(self, obj_data):
        line = ''
        line += "Name         {}\n".format(obj_data['info'][1])
        line += "Type         {}\n".format(obj_data['info'][2])
        line += "Created By   {}\n".format(obj_data['info'][5])
        line += "Narrative    {}\n".format(obj_data['info'][7])
        line += "Description  {}\n".format(obj_data['data']['description'])
        line += "Num Elements {}\n".format(str(len(obj_data['data']['elements'])))
        for ele in obj_data['data']['elements']:
            line += "  Element:   {}\n".format(ele)
        return line

    # Describe a GenomeSet
    #
    def readGenomeSet(self, obj_name, pyStr, format):
        myGS = pyStr['elements']
        line = ''
        if format == 'list':
            line = "Description for: " + obj_name + "\n"
        if format == 'tab':
            lst = ["Name", "ObjectID", "ScientName", "Size", "Source", "Domain", "Assembly Ref", "Features", "Contigs", "Pct. GC",
                   "Genetic Code", "Num CDS", "Num gene", "Num other", "Num rRNA", "Num tRNA"]
            line = "\t".join(lst) + "\n"
        if format == 'csv':
            lst = ["Name", "ObjectID", "ScientName", "Size", "Source", "Domain", "Assembly Ref", "Features", "Contigs", "Pct. GC",
                   "Genetic Code", "Num CDS", "Num gene", "Num other", "Num rRNA", "Num tRNA"]
            line = ",".join(lst) + "\n"

        for ele in myGS:
            genome = self.dfu.get_objects({'object_refs': [myGS[ele]['ref']]})
            line += self.getGenomeSet(myGS[ele]['ref'], genome['data'][0], format)

#        print "LINE:", line
        return line

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
