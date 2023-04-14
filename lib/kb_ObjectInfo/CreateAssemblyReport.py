import time
import os
from .CreateFasta_Report import CreateFasta

from installed_clients.DataFileUtilClient import DataFileUtil

class CreateAssemblyReport:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)

               
    def assembly_overview(self,obj_list):
        header = "OVERVIEW"
        this_list = [["Assembly Name","Type","Assembly Type"]]
        
        for assembly in obj_list['data']:
            name = "Assembly Data Object"
            type_assem = ''
            object_type = ''
            if 'info' in assembly:
                name = assembly['info'][1]
                object_type = assembly['info'][2]
            if 'type' in assembly['data']:
                type_assem = assembly['data']['type']

            this_list.extend([[name,object_type,type_assem]])
        
        return (header,this_list)
        
    def assembly_metadata(self,obj_list):
        header = "METADATA"
        this_list = [['Assembly Name','Assembly ID', 'DNA size', 'GC content', 'Number contigs',
                'FastA handle reference', 'MD5', 'Type', 'Taxon reference','Original filename']]
                
        list = ['assembly_id', 'dna_size', 'gc_content', 'num_contigs',
                'fasta_handle_ref', 'md5', 'type', 'taxon_ref']
                
        for assembly in obj_list['data']:
            name = "Assembly Data Object"
            if 'info' in assembly:
                name = assembly['info'][1]

            # Create the row for the one assembly
            assem_list = [name]
            for item in list:
                if item in assembly['data']:
                    assem_list.append(str(assembly['data'][item]))
                else:
                    assem_list.append(" ")

            if 'fasta_handle_info' in assembly['data'] and 'node_file_name' in assembly['data']['fasta_handle_info']:
                assem_list.append(assembly['data']['fasta_handle_info']['node_file_name'])
            else:
                assem_list.append(" ")
                
            # Add the row to the list that will be returned
            this_list.extend([assem_list])
            
        return (header,this_list)
        
    def assembly_dnabases(self,obj_list):
        header = "DNA Composition"
        this_list = [["Assembly Name","Total DNA Bases",
                        "A Count","A Percent","C Count","C Percent",
                        "G Count","G Percent","T Count","T Percent"]]
 
        for assembly in obj_list['data']:
            name = "Assembly Data Object"
            if 'info' in assembly:
                name = assembly['info'][1]
            dna_size = 1.0
            pct = 1.00
            if 'dna_size' in assembly['data']:
                dna_size = assembly['data']['dna_size']
            
            assem_list = [name,str(dna_size)]
            for base in ["A","C","G","T"]:
                pct = round(100 * assembly['data']['base_counts'][base] / dna_size,2)
                assem_list.append(str(assembly['data']['base_counts'][base]))
                assem_list.append(str(pct))
                
            this_list.extend([assem_list])
            
        return(header,this_list)
                
    def assembly_contigs(self,obj_list):
        header = "Contigs in the Assembly"
        this_list= [["Assembly Name","Contig Name","Length","GC content","Number of Ns","Contig ID","Description"]]
        
        for assembly in obj_list['data']:
            name = "Assembly Data Object"
            if 'info' in assembly:
                name = assembly['info'][1]
 
            if 'contigs' in assembly['data']:
                myContig = assembly['data']['contigs']
                for ctg in myContig:
                    list = ['length', 'gc_content', 'Ncount', 'contig_id', 'description']
                    ctg_list = [name,ctg]
                
                    for item in list:
                        if item in myContig[ctg]:
                            ctg_list.append(format(myContig[ctg][item]))
                        else:
                            ctg_list.append("")

                        this_list.append(ctg_list)
                        
        return(header,this_list)
        
                        
    def assembly_dna(self,obj_list,scratch):
        header = "Contig FastA files found in the download files."
        fasta_list = []
        dna_string = ""
        cf = CreateFasta(self.config)
        
        for assembly in obj_list['data']:
            name = "Assembly Data Object"
            input_ref = ''
            if 'info' in assembly:
                name = assembly['info'][1]
                input_ref = str(assembly['info'][6]) + "/" + str(assembly['info'][0]) + "/" + str(assembly['info'][4])

                fasta_list = cf.get_assembly_sequence(input_ref)
                report_path = os.path.join(scratch, name + '.fna')
            
#           Write the DNA string out to a Fasta file
                report_txt = open(report_path, "w")
                for dna_seq in fasta_list:
                    dna = "\n".join(dna_seq)+"\n"
                    report_txt.write(dna)
                    dna_string += dna
                report_txt.close()
            
        return(header)
   
