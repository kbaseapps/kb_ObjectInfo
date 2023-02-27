# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import csv

from Bio import SeqIO
from pprint import pprint, pformat
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from .CreateFasta_Report import CreateFasta
from .CreateFeatureLists_Report import CreateFeatureLists
from .CreateMultiGenomeReport import CreateMultiGenomeReport
from .Report_creator import Report_creator

#END_HEADER


class kb_ObjectInfo:
    '''
    Module Name:
    kb_ObjectInfo

    Module Description:
    A KBase module: kb_ObjectInfo
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_ObjectInfo"
    GIT_COMMIT_HASH = "25e0fc5060a406e7477beea677ec631031b2bdef"

    #BEGIN_CLASS_HEADER

#   rpt_list is a list of the rows that will form the output file
#   rpt_string is used for the PREVIEW for the user and for the html version of the file
#   htmltable is the html table version of rpt_list/rpt_string
#   report_path is the physical location of the output files
#   report_txt is an internal name used for referencing the report_path
#   rpt_writer is a writer object responsible for converting user data to delimited strings
#   report_format is tsv, csv, or something else selected by the user
#   rpt_delimiter is the delimiter used in the file. "\t" for tsv and ',' for csv, etc.

    def make_HTML(self,rpt_list):
        table = "<table>\n"

        # Create the table's row data
        for row in rpt_list:
            table += "  <tr>\n"
            for col in row:
                if col < '     ':
                    table += "  </tr>\n<table>\n"
                    
                table += "    <td>{0}</td>\n".format(col)
            table += "  </tr>\n"

        table += "</table>\n"
        return table

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
        self.workspaceURL = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = os.path.abspath(config['scratch'])
        self.config = config
        #END_CONSTRUCTOR
        pass

    def assembly_metadata_report(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of type "AssemblyMetadataReportParams" (A
           'typedef' can also be used to define compound or container
           objects, like lists, maps, and structures.  The standard KBase
           convention is to use structures, as shown here, to define the
           input and output of your function.  Here the input is a reference
           to the Assembly data object, a workspace to save output, and a
           length threshold for filtering.) -> structure: parameter
           "input_ref" of type "assembly_ref", parameter "workspace_name" of
           String, parameter "showContigs" of type "boolean" (A boolean. 0 =
           false, other = true.), parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN assembly_metadata_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting Assembly MetaData Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")
        
        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']
        if 'showContigs' not in params:
            raise ValueError('Parameter showContigs is not set in input arguments')
        showContigs_orig = params['showContigs']
        showContigs = None
        try:
            showContigs = int(showContigs_orig)
        except ValueError:
            raise ValueError('Cannot parse integer from showContigs parameter (' + str(showContigs_orig) + ')')
        if showContigs < 0:
            raise ValueError('showContigs parameter cannot be negative (' + str(showContigs) + ')')
        if showContigs > 1:
            raise ValueError('showContigs parameter cannot be greater than one (' + str(showContigs) + ')')

        report_format = params['report_format']
        rpt_delimiter = "\t"
        if report_format == 'csv':
            rpt_delimiter = ','

        # Step 3 - Get the data and save the output to a file.
        data_file_cli = DataFileUtil(self.callback_url)
        assembly = data_file_cli.get_objects({'object_refs': [input_ref]})
        name = "Assembly Data Object"
        object_type = ''
        if 'info' in assembly['data'][0]:
            name = assembly['data'][0]['info'][1]
            object_type = assembly['data'][0]['info'][2]
        assembly_metadata = assembly['data'][0]['data']

        rpt_list = []
        rpt_list = [["Name="+name," Type="+object_type],[""],["METADATA"]]

        dna_size = 1.0
        list = ['assembly_id', 'dna_size', 'gc_content', 'num_contigs',
                'fasta_handle_ref', 'md5', 'type', 'taxon_ref']
        for item in list:
            if item in assembly_metadata:
                rpt_list.append([item,str(assembly_metadata[item])])
                if item == 'dna_size':
                    dna_size = assembly_metadata['dna_size']

        if 'fasta_handle_info' in assembly_metadata and 'node_file_name' in assembly_metadata['fasta_handle_info']:
            rpt_list.append(["Original filename",assembly_metadata['fasta_handle_info']['node_file_name']])
            
        rpt_list.append([""])
        rpt_list.append(["DNA BASES","COUNTS","PERCENT"])
        pct = 1.00
        for base in assembly_metadata['base_counts']:
            pct = 100 * assembly_metadata['base_counts'][base] / dna_size
            rpt_list.append([base,str(assembly_metadata['base_counts'][base]),str(pct)])

        rpt_list.append([""])
        rpt_list.append(["CONTIGS IN THE ASSEMBLY"])
        rpt_list.append(["Name","Length","GC content","Number of Ns","Contig ID","Description"])

        if 'contigs' in assembly_metadata:
            myContig = assembly_metadata['contigs']
            for ctg in myContig:
                list = ['length', 'gc_content', 'Ncount', 'contig_id', 'description']
                ctg_list = [ctg]
                
                for item in list:
                    if item in myContig[ctg]:
                        ctg_list.append(format(myContig[ctg][item]))
                    else:
                        ctg_list.append("")

                rpt_list.append(ctg_list)

        rpt_delimiter = "\t"
        rpt_string = "Data Columns are tab-delimited\n"
        report_path = os.path.join(self.scratch, 'assembly_meta_tab_file.tsv')
        if report_format == 'csv':
            rpt_delimiter = ','
            rpt_string = "Data Columns are comma-delimited\n"
            report_path = os.path.join(self.scratch, 'assembly_meta_csv_file.csv')
        
        with open(report_path, mode='w') as report_txt:
            rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
            for rpt in rpt_list:
                rpt_writer.writerow(rpt)
                rpt_string += rpt_delimiter.join(rpt) + "\n"
        report_txt.close()

        fasta_list = []
        if showContigs:
            cf = CreateFasta(self.config)
            rpt_string += "\nFASTA of the DNA Sequences\n"
            fasta_list = cf.get_assembly_sequence(input_ref)
            report_path = os.path.join(self.scratch, 'assembly_metadata_file.fna')
            
#           Write the DNA string out to a Fasta file
            report_txt = open(report_path, "w")
            for dna_seq in fasta_list:
                dna = "\n".join(dna_seq)
                report_txt.write(dna)
                rpt_string += dna
            report_txt.close()

        report_path = os.path.join(self.scratch, 'assembly_metadata_file.html')
        report_txt = open(report_path, "w")
        report_txt.write("<pre>" + rpt_string + "</pre>")
        report_txt.close()

        cr = Report_creator(self.config)
        reported_output = cr.create_report(token, params['workspace_name'],
                                    rpt_string, self.scratch)

        output = {'report_name': reported_output['name'],
                           'report_ref': reported_output['ref']}
        mystr = pformat(output)
        logging.info(f"Returning: {mystr}")
        #END assembly_metadata_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method assembly_metadata_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]



    def genome_report(self, ctx, params):
        """
        :param params: instance of type "GenomeReportParams" -> structure:
           parameter "input_ref" of type "genome_ref", parameter
           "workspace_name" of String, parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN genome_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting Genome Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params=:\n{mystr}")
        
        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        data_file_cli = DataFileUtil(self.callback_url)
        genome = data_file_cli.get_objects({'object_refs': [input_ref]})
        genome_data = genome['data'][0]['data']

        report_format = params['report_format']
        rpt_list = []
        rpt_string = ''
        
        if report_format == 'tab':
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.delimitedTable(genome_data, 'tab', 'features')
            report_path = os.path.join(self.scratch, 'genome_tab_file.tsv')
        elif report_format == 'csv':
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.delimitedTable(genome_data, 'csv', 'features')
            report_path = os.path.join(self.scratch, 'genome_csv_file.csv')
        elif report_format == 'gff':
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.gff3(genome_data, 'features')
            report_path = os.path.join(self.scratch, 'genome_file.gff')
        elif report_format == 'fasta':
            cf = CreateFasta(self.config)
#           Before version 9 genomes, the cdss didn't exist
            if genome_data['cdss']:
                rpt_list = cf.create_Fasta_from_features(genome_data['cdss'])
            else:
                rpt_list = cf.create_Fasta_from_features(genome_data['features'])
            report_path = os.path.join(self.scratch, 'genome_file.faa')
        elif report_format == 'mRNA':
            cf = CreateFasta(self.config)
#           Before version 9 genomes, the cdss didn't exist
            if genome_data['cdss']:
                rpt_list = cf.create_Fasta_from_mRNA(genome_data['cdss'])
            else:
                rpt_list = cf.create_Fasta_from_mRNA(genome_data['features'])
            report_path = os.path.join(self.scratch, 'genome_mRNA_file.fna')
        elif report_format == 'DNA':
            cf = CreateFasta(self.config)
            report_path = os.path.join(self.scratch, 'genome_dna_file.fna')
            if 'assembly_ref' in genome_data:
                input_ref = genome_data['assembly_ref']
                rpt_list = (cf.get_assembly_sequence(input_ref))
            elif 'assembly_ref' in genome['data'][0]['info'][10]:
                input_ref = genome['data'][0]['info'][10]['assembly_ref']
                rpt_list = (cf.get_assembly_sequence(input_ref))
            elif 'contigset_ref' in genome_data:
                input_ref = genome_data['contigset_ref']
                rpt_list = (cf.get_assembly_sequence(input_ref))
            else:
                rpt_string += 'Did not find the Assembly Reference\n'
        else:
            raise ValueError('Invalid report option.' + str(report_format))

#       The rpt_list only exists if the output is tab or comma delimited
#       If the output is DNA or mRNA/Fasta, don't reset the string or use csv.writer
        if rpt_list:
            rpt_string = ''
            rpt_delimiter = "\t"
            if report_format == 'csv':
                rpt_delimiter = ','
                
            with open(report_path, mode='w') as report_txt:
                rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
                if report_format == 'csv':
                    rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
                for rpt in rpt_list:
                    rpt_writer.writerow(rpt)
                    rpt_string += rpt_delimiter.join(rpt) + "\n"
        else:
            report_txt = open(report_path, "w")
            report_txt.write(rpt_string)
            report_txt.close()
        
        report_path = os.path.join(self.scratch, 'text_file.html')
        report_txt = open(report_path, "w")
        htmltable = self.make_HTML(rpt_list)
        #report_txt.write("<pre>" + string + "</pre>")
        report_txt.write(htmltable)
        report_txt.close()

        cr = Report_creator(self.config)
        reported_output = cr.create_report(token, params['workspace_name'],
                                    rpt_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        mystr = pformat(output)
        logging.info(f"Returning: {mystr}")
        #END genome_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method genome_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def genomeset_report(self, ctx, params):
        """
        :param params: instance of type "GenomeSetReportParams" -> structure:
           parameter "input_ref" of type "genomeset_ref", parameter
           "workspace_name" of String, parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN genomeset_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting Genome Set Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        data_file_cli = DataFileUtil(self.callback_url)
        genomeset = data_file_cli.get_objects({'object_refs': [input_ref]})
        genome_name = genomeset['data'][0]['info'][1]
        genomeset_data = genomeset['data'][0]['data']

        report_format = params['report_format']

        rpt_list = []
        multi_fasta = []
        
        if report_format == 'tab':
            gsr = CreateMultiGenomeReport(self.config)
            rpt_list = gsr.readGenomeSet(genome_name, genomeset_data, 'tab')
            report_path = os.path.join(self.scratch, 'genomeset_tab_file.tsv')
        elif report_format == 'csv':
            gsr = CreateMultiGenomeReport(self.config)
            rpt_list = gsr.readGenomeSet(genome_name, genomeset_data, 'csv')
            report_path = os.path.join(self.scratch, 'genomeset_cvs_file.csv')
        elif report_format == 'list':
            gsr = CreateMultiGenomeReport(self.config)
            rpt_list = gsr.readGenomeSet(genome_name, genomeset_data, 'list')
            report_path = os.path.join(self.scratch, 'genomeset_list_file.tsv')
        elif report_format == 'meta':
            gsr = CreateMultiGenomeReport(self.config)
            rpt_list = gsr.getGenomeSetMeta(genomeset['data'][0])
            report_path = os.path.join(self.scratch, 'genomeset_meta_file.tsv')
        elif report_format == 'fasta':
            gsr = CreateMultiGenomeReport(self.config)
            rpt_list = [["Assembly Reference","Scientific Name","File Name"]]
            
#           Get the list of assembly IDs for the genomeSet
            assembly_list = gsr.getAssemblyRef(genomeset['data'][0])
            
#           For each assembly, get it's info and sequence

            for assembly in assembly_list:
                cf = CreateFasta(self.config)
                assembly_ref, sci_name = assembly.split(':')
                file_name ='G'+assembly_ref.replace('/', '_')+'.fna'
                rpt_list.append([assembly_ref,sci_name,file_name])
                
#               Save the Fasta sequences to an individual genome file
                fasta_list = cf.get_assembly_sequence(assembly_ref)
                multi_fasta.extend(fasta_list)
                report_path = os.path.join(self.scratch, file_name)
                
                with open(report_path, mode='w') as report_txt:
                    rpt_writer = csv.writer(report_txt, delimiter="\n", quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    for dna in fasta_list:
                        rpt_writer.writerow(dna)
                report_txt.close()
                
#           Set the report path for the summary table (Don't overwrite the last dna file)
            report_path = os.path.join(self.scratch, 'genomeset_fasta_file.txt')
        else:
            raise ValueError('Invalid report option.' + str(report_format))
    
        rpt_list.extend(multi_fasta)
    
        rpt_string = ''
        rpt_delimiter = "\t"
        if report_format == 'tsv':
            rpt_string = "Data Columns are tab-delimited\n"
        elif report_format == 'csv':
            rpt_string = "Data Columns are comma-delimited\n"
            rpt_delimiter = ','
            
        if rpt_list:
            print ("DEBUG: genomeSet ", report_format, rpt_list[0:3])
            
            with open(report_path, mode='w') as report_txt:
                rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
                if report_format == 'csv':
                    rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
                
                for rpt in rpt_list:
                    rpt_writer.writerow(rpt)
                    rpt_string += rpt_delimiter.join(rpt) + "\n"
        
        report_path = os.path.join(self.scratch, 'text_file.html')
        report_txt = open(report_path, "w")
        htmltable = self.make_HTML(rpt_list)
        report_txt.write(htmltable)
        report_txt.close()

        cr = Report_creator(self.config)
        reported_output = cr.create_report(token, params['workspace_name'],
                                           rpt_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        mystr = pformat(output)
        logging.info(f"Returning: {mystr}")
        #END genomeset_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method genomeset_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def domain_report(self, ctx, params):
        """
        :param params: instance of type "DomainReportParams" -> structure:
           parameter "input_ref" of type "domain_ref", parameter
           "evalue_cutoff" of Double, parameter "workspace_name" of String,
           parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN domain_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting Domain Annotation Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        data_file_cli = DataFileUtil(self.callback_url)
        domain_anno = data_file_cli.get_objects({'object_refs': [input_ref]})
        domain_data = domain_anno['data'][0]['data']

        evalue_cutoff = float(params['evalue_cutoff'])
        report_format = params['report_format']

        rpt_list1 = []
        rpt_list2 = []

        rpt_delimiter = "\t"
        if report_format == 'csv':
            rpt_delimiter = ','

        if report_format == 'tab':
            cf = CreateFeatureLists(self.config)
            report_path1 = os.path.join(self.scratch, 'domain_annotation_tab_list.tsv')
            report_path2 = os.path.join(self.scratch, 'domain_annotation_tab_count.tsv')
        elif report_format == 'csv':
            cf = CreateFeatureLists(self.config)
            report_path1 = os.path.join(self.scratch, 'domain_annotation_csv_list.csv')
            report_path2 = os.path.join(self.scratch, 'domain_annotation_csv_count.csv')
        else:
            raise ValueError('Invalid report option.' + str(report_format))

        rpt_list1 = cf.readDomainAnnList(domain_data, rpt_delimiter, evalue_cutoff)
        rpt_list2 = cf.readDomainAnnCount(domain_data,rpt_delimiter, evalue_cutoff)
    
        rpt_string = ''
        
        if rpt_list1 or rpt_list2:
            if rpt_list1:
                with open(report_path1, mode='w') as report_txt:
                    rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    if report_format == 'csv':
                        rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
                    for rpt in rpt_list1:
                        rpt_writer.writerow(rpt)
        
            if rpt_list2:
                with open(report_path2, mode='w') as report_txt:
                    rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    if report_format == 'csv':
                        rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
                    for rpt in rpt_list2:
                        rpt_writer.writerow(rpt)
                        rpt_string += rpt_delimiter.join(rpt) + "\n"

        report_path = os.path.join(self.scratch, 'text_file.html')
        report_txt = open(report_path, "w")
        htmltable = self.make_HTML(rpt_list2)
        report_txt.write(htmltable)
        htmltable = self.make_HTML(rpt_list1)
        report_txt.write(htmltable)
        report_txt.close()

        cr = Report_creator(self.config)

        reported_output = cr.create_report(token, params['workspace_name'],
                                           rpt_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        mystr = pformat(output)
        logging.info(f"Returning: {mystr}")
        #END domain_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method domain_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def tree_report(self, ctx, params):
        """
        :param params: instance of type "TreeReportParams" -> structure:
           parameter "input_ref" of type "tree_ref", parameter
           "workspace_name" of String, parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN tree_report
        #END tree_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method tree_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def featseq_report(self, ctx, params):
        """
        :param params: instance of type "FeatSeqReportParams" -> structure:
           parameter "input_ref" of type "featseq_ref", parameter
           "workspace_name" of String, parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN featseq_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting FeatureSeq/SequenceSet Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.
        
        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']
        
        cf = CreateFeatureLists(self.config)
        data_file_cli = DataFileUtil(self.callback_url)
        setseq = data_file_cli.get_objects({'object_refs': [input_ref]})
        setseq_data = setseq['data'][0]['data']

        report_format = params['report_format']

        rpt_list = []
        rpt_delimiter = "\t"
        
        if report_format == 'tab':
            rpt_list = cf.readFeatSeq(setseq_data, 'tab')
            report_path = os.path.join(self.scratch, 'sequence_set_tab_list.tsv')
        elif report_format == 'csv':
            rpt_list = cf.readFeatSeq(setseq_data, 'csv')
            report_path = os.path.join(self.scratch, 'sequence_set_csv_list.csv')
            rpt_delimiter = ','
        else:
            raise ValueError('Invalid report option.' + str(report_format))

        rpt_string = ''
        if rpt_list:
            with open(report_path, mode='w') as report_txt:
                rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
                if report_format == 'csv':
                    rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
                
                for rpt in rpt_list:
                    rpt_writer.writerow(rpt)
                    rpt_string += rpt_delimiter.join(rpt) + "\n"

        report_path = os.path.join(self.scratch, 'text_file.html')
        report_txt = open(report_path, "w")
        htmltable = self.make_HTML(rpt_list)
        report_txt.write(htmltable)
        report_txt.close()

        cr = Report_creator(self.config)

        reported_output = cr.create_report(token, params['workspace_name'],
                                           rpt_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        mystr = pformat(output)
        logging.info(f"Returning: {mystr}")
        #END featseq_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method featseq_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def protcomp_report(self, ctx, params):
        """
        :param params: instance of type "ProtCompReportParams" -> structure:
           parameter "input_ref" of type "protcomp_ref", parameter
           "workspace_name" of String, parameter "report_format" of String
        :returns: instance of type "ReportResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN protcomp_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting ProteomeComparison Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        data_file_cli = DataFileUtil(self.callback_url)
        protcomp = data_file_cli.get_objects({'object_refs': [input_ref]})
        protcomp_data = protcomp['data'][0]['data']

        report_format = params['report_format']

        rpt_list = []
        rpt_delimiter = "\t"
        
        if report_format == 'tab':
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.readProtComp(protcomp_data, 'tab')
            report_path = os.path.join(self.scratch, 'protcomp_tab_list.tsv')
        elif report_format == 'csv':
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.readProtComp(protcomp_data, 'csv')
            report_path = os.path.join(self.scratch, 'protcomp_csv_list.csv')
            rpt_delimiter = ','
        else:
            raise ValueError('Invalid report option.' + str(report_format))

        #Write to csv/tsv file. Create string for html and report output
        rpt_string = ''
        if rpt_list:
            with open(report_path, mode='w') as report_txt:
                rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
                if report_format == 'csv':
                    rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
                for rpt in rpt_list:
                    rpt_writer.writerow(rpt)
                    rpt_string += rpt_delimiter.join(rpt) + "\n"
        
        report_path = os.path.join(self.scratch, 'text_file.html')
        report_txt = open(report_path, "w")
        htmltable = self.make_HTML(rpt_list)
        report_txt.write(htmltable)
        report_txt.close()

        cr = Report_creator(self.config)

        reported_output = cr.create_report(token, params['workspace_name'],
                                           rpt_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        #END protcomp_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method protcomp_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
        
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
