# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import csv

from Bio import SeqIO
from pprint import pprint, pformat
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from .CreateFasta_Report import CreateFasta
from .CreateFeatureLists_Report import CreateFeatureLists
from .CreateMultiGenomeReport import CreateMultiGenomeReport
from .CreateAssemblyReport import CreateAssemblyReport
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
    VERSION = "1.2.2"
    GIT_URL = "https://github.com/kbaseapps/kb_ObjectInfo"
    GIT_COMMIT_HASH = "e5ad5301c6e0f8e6be92123d95efc3d3cdf2cd81"

    #BEGIN_CLASS_HEADER

#   report_list is a list of the rows that will form the output file
#   report_string is used for the PREVIEW for the user and for the html version of the file
#   html_table is the html table version of report_list/report_string
#   report_path is the physical location of the output files
#   report_txt is an internal name used for referencing the report_path
#   report_writer is a writer object responsible for converting user data to delimited strings
#   report_format is tsv, csv, or something else selected by the user
#   report_delimiter is the delimiter used in the file. "\t" for tsv and ',' for csv, etc.

    def make_HTML_table_from_list(self,report_list,style):
        '''Take a list of lists needed for csv writer and create an HTML table'''
        table = "<table style=\"border: 1px solid black;\">\n"

        # Create the table's row data
        for i, row in enumerate(report_list):
            table += "  <tr style=\"border: 1px solid black;\">\n"
            for j, col in enumerate(row):
                if i==0 and style == 'col_header':
                    table += "    <th style=\"border: 1px solid black;\">{0}</th>\n".format(col)
                elif j==0 and style == 'row_header':
                    table += "    <th style=\"border: 1px solid black;\">{0}</th>\n".format(col)
                else:
                    table += "    <td style=\"border: 1px solid black;\">{0}</td>\n".format(col)
                    
            table += "  </tr>\n"

        table += "</table>\n"
        return table

    def write_list_to_text_file(self,report_list,report_path,report_delimiter):
        '''Take a list of list and create tab or csv delimited file with csv writer'''
        report_string = ''
        report_path = os.path.join(self.scratch, report_path)
        with open(report_path, mode='w', encoding='utf-8') as report_txt:
            report_writer = csv.writer(report_txt, delimiter=report_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
            if report_delimiter == ',':
                    report_writer = csv.writer(report_txt, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
            for rpt in report_list:
                report_writer.writerow(rpt)
                report_string += report_delimiter.join(rpt) + "\n"
            
        return report_string
        
    def read_columns_from_text_file(self,report_path):
        '''Read a tab-delimited file into a list. Used here to import a GFF file'''
        report_list = []

        with open(report_path, mode='r', encoding="utf-8") as report_txt:
            for line in report_txt:
                cols = line.split("\t")
                report_list.append(cols)
            
        return report_list
        
    def truncate_report_string(self,report_string):
        '''Take the report string and truncate each of the strings for readability'''
        '''Based on PR review suggestion from Erik Pearson'''
        report_string_list = report_string.split("\n")
        new_report_string = ""
        max = 70
        for line in report_string_list:
            if len(line) > max:
                new_report_string += line[0:68] + '...'
            else:
                new_report_string += line
                
        return new_report_string
            
            
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
        self.gfu = GenomeFileUtil(self.callback_url)
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
           String, parameter "show_contigs" of type "boolean" (A boolean. 0 =
           false, other = true.)
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
        logging.info('Starting Assembly Object Info. ')
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")
        
        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']
        if 'show_contigs' not in params:
            raise ValueError('Parameter show_contigs is not set in input arguments')
        show_contigs_orig = params['show_contigs']
        show_contigs = None
        try:
            show_contigs = int(show_contigs_orig)
        except ValueError:
            raise ValueError('Cannot parse integer from show_contigs parameter (' + str(show_contigs_orig) + ')')
        if show_contigs < 0:
            raise ValueError('show_contigs parameter cannot be negative (' + str(show_contigs) + ')')
        if show_contigs > 1:
            raise ValueError('show_contigs parameter cannot be greater than one (' + str(show_contigs) + ')')
        
        create_assembly_report = CreateAssemblyReport(self.config)
        assembly = self.dfu.get_objects({'object_refs': [input_ref]})
        report_list = []
        this_list = []
        
        html_report_path = os.path.join(self.scratch, 'assembly_metadata_file.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")
        
        (header,this_list) = create_assembly_report.assembly_overview(assembly)
        
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[header]])
        report_list.extend(this_list)
        
        (header,this_list) = create_assembly_report.assembly_metadata(assembly)
        
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[],[header]])
        report_list.extend(this_list)
        
        (header,this_list) = create_assembly_report.assembly_dnabases(assembly)
                
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[],[header]])
        report_list.extend(this_list)
        
        (header,this_list) = create_assembly_report.assembly_contigs(assembly)
                
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[],[header]])
        report_list.extend(this_list)
        html_report_txt.close()

        report_string = self.write_list_to_text_file(report_list,'assembly_meta_tab_file.tsv',"\t")
        self.write_list_to_text_file(report_list,'assembly_meta_csv_file.csv',",")
                                          
        if show_contigs:
            (header) = create_assembly_report.assembly_dna(assembly,self.scratch)
            report_list.extend([[],[header]])
        
        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)
                                           
        output = {'report_name': reported_output['name'],
                           'report_ref': reported_output['ref']}
        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
        #END assembly_metadata_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method assembly_metadata_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def assemblyset_report(self, ctx, params):
        """
        :param params: instance of type "AssemblySetReportParams" ->
           structure: parameter "input_ref" of type "assemblyset_ref",
           parameter "workspace_name" of String, parameter "show_contigs" of
           type "boolean" (A boolean. 0 = false, other = true.)
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
        #BEGIN assemblyset_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting AssemblySet Object Info. ')
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")

        input_ref = params['input_ref']
        if 'show_contigs' not in params:
            raise ValueError('Parameter show_contigs is not set in input arguments')
        show_contigs_orig = params['show_contigs']
        show_contigs = None
        try:
            show_contigs = int(show_contigs_orig)
        except ValueError:
            raise ValueError('Cannot parse integer from show_contigs parameter (' + str(show_contigs_orig) + ')')
        if show_contigs < 0:
            raise ValueError('show_contigs parameter cannot be negative (' + str(show_contigs) + ')')
        if show_contigs > 1:
            raise ValueError('show_contigs parameter cannot be greater than one (' + str(show_contigs) + ')')

        assemblyset = self.dfu.get_objects({'object_refs': [input_ref]})
        assemblylist = assemblyset['data'][0]['data']['items']
        num_assem = len(assemblylist)
        create_assembly_report = CreateAssemblyReport(self.config)
        
        assemref = []
        
        for i in range(num_assem):
            assemref.append(assemblylist[i]['ref'])
        
        html_report_path = os.path.join(self.scratch, 'assembly_metadata_file.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")
        
        this_list = []
        report_list = []
        obj_list = self.dfu.get_objects({'object_refs': assemref})
        
        (header,this_list) = create_assembly_report.assembly_overview(obj_list)
        
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[header]])
        report_list.extend(this_list)

        (header,this_list) = create_assembly_report.assembly_metadata(obj_list)
        
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[],[header]])
        report_list.extend(this_list)
        
        (header,this_list) = create_assembly_report.assembly_dnabases(obj_list)
                
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[],[header]])
        report_list.extend(this_list)

        (header,this_list) = create_assembly_report.assembly_contigs(obj_list)
                
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(html_table)
        report_list.extend([[],[header]])
        report_list.extend(this_list)
        html_report_txt.close()
                                  
        if show_contigs:
            (header) = create_assembly_report.assembly_dna(obj_list,self.scratch)
            report_list.extend([[],[header]])

        report_string = self.write_list_to_text_file(report_list,'assembly_set_tab_file.tsv',"\t")
        self.write_list_to_text_file(report_list,'assembly_set_csv_file.csv',",")
                                          
        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)
                                           
        output = {'report_name': reported_output['name'],
                    'report_ref': reported_output['ref']}
        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
        #END assemblyset_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method assemblyset_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def genome_report(self, ctx, params):
        """
        :param params: instance of type "GenomeReportParams" -> structure:
           parameter "input_ref" of type "genome_ref", parameter
           "workspace_name" of String, parameter "list_coding" of type
           "boolean" (A boolean. 0 = false, other = true.), parameter
           "list_GFF" of type "boolean" (A boolean. 0 = false, other =
           true.), parameter "fasta_AA" of type "boolean" (A boolean. 0 =
           false, other = true.), parameter "fasta_mRNA" of type "boolean" (A
           boolean. 0 = false, other = true.), parameter "show_DNA" of type
           "boolean" (A boolean. 0 = false, other = true.)
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
        my_string = pformat(params)
        logging.info(f"Params=:\n{my_string}")
        
        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        genome = self.dfu.get_objects({'object_refs': [input_ref]})
        genome_data = genome['data'][0]['data']

        report_list = []
        report_string = ''
        report_path = ''
        
        if params['list_coding']:
            create_features = CreateFeatureLists(self.config)
            report_list = create_features.delimitedTable(genome_data, 'features')
            if report_list:
                report_string += self.write_list_to_text_file(report_list,'genome_tab_file.tsv',"\t")
                self.write_list_to_text_file(report_list,'genome_csv_file.csv',",")
                
            html_report_path = os.path.join(self.scratch, 'genome_features.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(report_list,'col_header')
            html_report_txt.write("<h1>PROTEIN CODING FEATURES</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()
                
        if params['list_GFF']:
            
            gff_return = self.gfu.genome_to_gff({'genome_ref':input_ref,'is_gtf':0,'target_dir':self.scratch})
            report_string += "Created " + gff_return['file_path'] + "\n"
            report_list += self.read_columns_from_text_file(gff_return['file_path'])
            
            html_report_path = os.path.join(self.scratch, 'genome_GFF.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(report_list,'no_header')
            html_report_txt.write("<h1>GFF OUTPUT</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()
            
        if params['fasta_AA']:
            create_fasta = CreateFasta(self.config)
#           Before version 9 genomes, the cdss didn't exist
            if genome_data['cdss']:
                report_list = create_fasta.create_fasta_from_features(genome_data['cdss'])
            else:
                report_list = create_fasta.create_fasta_from_features(genome_data['features'])
            
            dna_string = self.write_list_to_text_file(report_list,'genome_file.faa',"\n")
            report_string += dna_string[0:200] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
            html_report_path = os.path.join(self.scratch, 'genome_protein_AA.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_report_txt.write("<h1>PROTEIN CODING AMINO ACIDS</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")
            html_report_txt.close()
            
        if params['fasta_mRNA']:
            create_fasta = CreateFasta(self.config)
#           Before version 9 genomes, the cdss didn't exist
            if genome_data['cdss']:
                report_list = create_fasta.create_fasta_from_mRNA(genome_data['cdss'])
            else:
                report_list = create_fasta.create_fasta_from_mRNA(genome_data['features'])
                
            dna_string = self.write_list_to_text_file(report_list,'genome_mRNA_file.fna',"\n")
            report_string += dna_string[0:200] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
            html_report_path = os.path.join(self.scratch, 'genome_protein_mRNA.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_report_txt.write("<h1>PROTEIN CODING mRNA</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")
            html_report_txt.close()
        
        if params['show_DNA']:
            create_fasta = CreateFasta(self.config)
            report_path = os.path.join(self.scratch, 'genome_dna_file.fna')
            if 'assembly_ref' in genome_data:
                input_ref = genome_data['assembly_ref']
                report_list = (create_fasta.get_assembly_sequence(input_ref))
            elif 'assembly_ref' in genome['data'][0]['info'][10]:
                input_ref = genome['data'][0]['info'][10]['assembly_ref']
                report_list = (create_fasta.get_assembly_sequence(input_ref))
            elif 'contigset_ref' in genome_data:
                input_ref = genome_data['contigset_ref']
                report_list = (create_fasta.get_assembly_sequence(input_ref))
            else:
                report_string += 'Did not find the Assembly Reference\n'
                
            dna_string = self.write_list_to_text_file(report_list,'genome_dna_file.fna',"\n")
            report_string += dna_string[0:200] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
            html_report_path = os.path.join(self.scratch, 'genome_DNA.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_report_txt.write("<h1>DNA</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")
            html_report_txt.close()

        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                    report_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
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
           "workspace_name" of String, parameter "show_genomes" of type
           "boolean" (A boolean. 0 = false, other = true.), parameter
           "list_coding" of type "boolean" (A boolean. 0 = false, other =
           true.), parameter "list_GFF" of type "boolean" (A boolean. 0 =
           false, other = true.), parameter "list_GBK" of type "boolean" (A
           boolean. 0 = false, other = true.), parameter "fasta_AA" of type
           "boolean" (A boolean. 0 = false, other = true.), parameter
           "fasta_mRNA" of type "boolean" (A boolean. 0 = false, other =
           true.), parameter "show_DNA" of type "boolean" (A boolean. 0 =
           false, other = true.)
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
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        genome_set = self.dfu.get_objects({'object_refs': [input_ref]})
        genome_name = genome_set['data'][0]['info'][1]
        genome_set_data = genome_set['data'][0]['data']

        report_list = []
        multi_fasta = []
        report_string = ''
        
#       Always do the list of genomes - short and ensures that there is an HTML link
        genome_set_report = CreateMultiGenomeReport(self.config)
        report_list = genome_set_report.readGenomeSet(genome_name, genome_set_data)
        html_table = self.make_HTML_table_from_list(report_list,'col_header')
        html_report_path = os.path.join(self.scratch, 'genomeset_comparison.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")
        html_report_txt.write("<h1>GENOME COMPARISON</h1>")
        html_report_txt.write(html_table)
        html_report_txt.close()
        
        if params['show_genomes']:
            report_string += self.write_list_to_text_file(report_list,'genomeset_tab_file.tsv',"\t")
            self.write_list_to_text_file(report_list,'genomeset_cvs_file.csv',",")
        report_string += "\n\n"
        
        if params['list_coding']:
            create_features = CreateFeatureLists(self.config)
            my_genome_set = genome_set_data['elements']
            report_list = []
            for element in my_genome_set:
                genome = self.dfu.get_objects({'object_refs': [my_genome_set[element]['ref']]})
                genome_data = genome['data'][0]['data']

                report_list += create_features.delimitedTable(genome_data, 'features')
                
            if report_list:
                report_string += self.write_list_to_text_file(report_list,'genome_tab_file.tsv',"\t")
                report_string += "\n\n"
                self.write_list_to_text_file(report_list,'genome_csv_file.csv',",")
                
            html_report_path = os.path.join(self.scratch, 'genome_features.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(report_list,'col_header')
            html_report_txt.write("<h1>PROTEIN CODING FEATURES</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()
                
        if params['list_GFF']:
            my_genome_set = genome_set_data['elements']
    
            for element in my_genome_set:
                gff_return = self.gfu.genome_to_gff({'genome_ref':my_genome_set[element]['ref'],'is_gtf':0,'target_dir':self.scratch})
                report_string += "Created " + gff_return['file_path'] + "\n\n\n"
             
        if params['list_GBK']:
            my_genome_set = genome_set_data['elements']
    
            for element in my_genome_set:
                gbk_return = self.gfu.genome_to_genbank({'genome_ref':my_genome_set[element]['ref'],'target_dir':self.scratch})
                file_path = gbk_return['genbank_file']['file_path'].replace('/kb/module/work/tmp/','')
                report_string += "Created " + file_path + "\n\n\n"
            
        if params['fasta_AA']:
            my_genome_set = genome_set_data['elements']
    
            for element in my_genome_set:
                cds_return = self.gfu.genome_proteins_to_fasta({'genome_ref':my_genome_set[element]['ref']})
                report_string += "Created " + cds_return['file_path'] + "\n\n\n"
               
        if params['fasta_mRNA']:
            my_genome_set = genome_set_data['elements']
    
            for element in my_genome_set:
                mrna_return = self.gfu.genome_features_to_fasta({'genome_ref':my_genome_set[element]['ref']})
                report_string += "Created " + mrna_return['file_path'] + "\n\n\n"
 
        if params['show_DNA']:
            report_list = [["Assembly Reference","Scientific Name","File Name"]]
            
#           Get the list of assembly IDs for the genomeSet
            assembly_list = genome_set_report.getAssemblyRef(genome_set['data'][0])
            
#           For each assembly, get it's info and sequence

            for assembly in assembly_list:
                create_fasta = CreateFasta(self.config)
                assembly_ref, sci_name = assembly.split(':')
                file_name ='G'+assembly_ref.replace('/', '_')+'.fna'
                report_list.append([assembly_ref,sci_name,file_name])
                
#               Save the Fasta sequences to an individual genome file
                fasta_list = create_fasta.get_assembly_sequence(assembly_ref)
                multi_fasta.extend(fasta_list)
                report_path = os.path.join(self.scratch, file_name)
                
                with open(report_path, mode='w', encoding="utf-8") as report_txt:
                    report_writer = csv.writer(report_txt, delimiter="\n", quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    for dna in fasta_list:
                        report_writer.writerow(dna)
                report_txt.close()

#           Set the report path for the summary table (Don't overwrite the last dna file)
            report_string += self.write_list_to_text_file(report_list,'genomeset_DNA_meta_file.tsv',"\t")
            self.write_list_to_text_file(report_list,'genomeset_DNA_meta_file.csv',",")
            
#           The fasta not included in the HTML due to length
            html_report_path = os.path.join(self.scratch, 'genomeset_DNA_meta.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(report_list,'col_header')
            html_report_txt.write("<h1>GENOME DOWNLOAD METADATA</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()

        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
        #END genomeset_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method genomeset_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def genomecomp_report(self, ctx, params):
        """
        :param params: instance of type "GenomeCompReportParams" ->
           structure: parameter "input_ref" of type "genomecomp_ref",
           parameter "workspace_name" of String
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
        #BEGIN genomecomp_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting GenomeComparison Object Info. ')
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")
        
        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']
        
        create_features = CreateFeatureLists(self.config)
        data_file_cli = DataFileUtil(self.callback_url)
        gencomp = self.dfu.get_objects({'object_refs': [input_ref]})
        gencomp_data = gencomp['data'][0]['data']
        
        report_list1 = []
        report_list2 = []
        create_features = CreateFeatureLists(self.config)
        report_string = ''
        html_report_path = os.path.join(self.scratch, 'pangenome_genome_comparison.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")
#
#       OVERVIEW
#
        name = gencomp['data'][0]['info'][1]
        num_genomes = 'unk'
        if "Number genomes" in gencomp['data'][0]['info'][10]:
            num_genomes = gencomp['data'][0]['info'][10]["Number genomes"]

        overview_list = [["Name",name],
                        ["Number of Genomes",num_genomes],
                        ["Pangenome Reference", gencomp_data["pangenome_ref"]],
                        ["Comparison Name", gencomp_data["name"]],
                        ["Number of Core Families", str(gencomp_data["core_families"])],
                        ["Number of Core Functions", str(gencomp_data["core_functions"])] ]
        if overview_list:
            report_string += "\nOVERVIEW\n"
            report_string += self.write_list_to_text_file(overview_list,'gencomp_overview_tab.tsv',"\t")
            self.write_list_to_text_file(overview_list,'gencomp_overview_csv.csv',",")
                
            html_table = self.make_HTML_table_from_list(overview_list,'row_header')
            html_report_txt.write("<h1>OVERVIEW</h1>")
            html_report_txt.write(html_table)
#
#       GENOMES
#
        genome_data = gencomp_data["genomes"]
        genome_dict = {}
        genome_list = [["Name","ID","Taxonomy","Number of Families","Number of Features","Number of Functions"]]
        sim_list = [["Genome1","Genome2","Number Families in Common","Number Functions in Common"]]
        for gen in genome_data:
            genome_dict[gen["id"]] = gen["name"]
            genome_list.append([gen["name"],gen["id"],gen["taxonomy"],
                            str(gen["families"]),str(gen["features"]),str(gen["functions"])])
            for g2 in gen["genome_similarity"]:
                sim_list.append([gen["id"],g2,
                str(gen["genome_similarity"][g2][0]),
                str(gen["genome_similarity"][g2][1])])
                
        if genome_list:
            report_string += "\nGENOMES\n"
            report_string += self.write_list_to_text_file(genome_list,'gencomp_genomes_tab.tsv',"\t")
            self.write_list_to_text_file(genome_list,'gencomp_genomes_csv.csv',",")
   
            html_table = self.make_HTML_table_from_list(genome_list,'col_header')
            html_report_txt.write("<h1>GENOMES</h1>")
            html_report_txt.write(html_table)
                
        if sim_list:
            report_string += "\nSIMILARITY LIST\n"
            report_string += self.write_list_to_text_file(sim_list,'gencomp_genomes_similarity_tab.tsv',"\t")
            self.write_list_to_text_file(sim_list,'gencomp_genomes_similarity_csv.csv',",")
   
            html_table = self.make_HTML_table_from_list(sim_list,'col_header')
            html_report_txt.write("<h1>SIMILARITY LIST</h1>")
            html_report_txt.write(html_table)
            
        html_report_txt.close()
 #
 #      FAMILIES
 #
        family_data = gencomp_data["families"]
        family_list = [["Family","Number of Genomes","Core","Fraction consistent","Fraction Genomes","Type","Protein Translation"]]
        fam_list = [["Family","Genome Name","Genome ID","Gene","Score"]]
        for fam in family_data:
            family_list.append([fam["id"],
                str(fam["number_genomes"]),
                str(fam["core"]),
                str(fam["fraction_consistent_annotations"]),
                str(fam["fraction_genomes"]),
                fam["type"],fam["protein_translation"]])
            for f2 in fam["genome_features"]:
                for i in fam["genome_features"][f2]:
                    # Leave out i[1] because it is an unknown list
                    fam_list.append([fam["id"],genome_dict[f2],f2,i[0],str(i[2])])

        if family_list:
            report_string += "\nFAMILIES\n"
            report_string += self.write_list_to_text_file(family_list,'gencomp_familes_tab.tsv',"\t")
            self.write_list_to_text_file(family_list,'gencomp_families_csv.csv',",")
                        
            html_report_path = os.path.join(self.scratch, 'gencomp_family_list.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(family_list,'col_header')
            html_report_txt.write("<h1>FAMILIES</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()
                
        if fam_list:
            report_string += "\nGENES IN FAMILIES\n"
            report_string += self.write_list_to_text_file(fam_list,'gencomp_familes2_tab.tsv',"\t")
            self.write_list_to_text_file(fam_list,'gencomp_families2_csv.csv',",")
                                
            html_report_path = os.path.join(self.scratch, 'gencomp_family_genes.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(fam_list,'col_header')
            html_report_txt.write("<h1>GENES IN FAMILIES</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()
#
#       FUNCTIONS
#
        function_data = gencomp_data["functions"]
        function_list = [["ID","Number of Genomes","Core","Fraction Consistent Families","Fraction Genomes",
                    "Primary Class","Subclass","Subsystem"]]
        fun_list = [["Function","Genome Name","Genome ID","Gene","Unk","Score"]]
        for fun in function_data:
            function_list.append([fun["id"],
                str(fun["number_genomes"]),
                str(fun["core"]),
                str(fun["fraction_consistent_families"]),
                str(fun["fraction_genomes"]),
                fun["primclass"],fun["subclass"],fun["subsystem"]])
            for f2 in fun["genome_features"]:
                for i in fun["genome_features"][f2]:
                    fun_list.append([fun["id"],genome_dict[f2],f2,i[0],str(i[1]),str(i[2])])

        if function_list:
            report_string += "\nFUNCTIONS\n"
            report_string += self.write_list_to_text_file(function_list,'gencomp_functions_tab.tsv',"\t")
            self.write_list_to_text_file(function_list,'gencomp_functions_csv.csv',",")
                                 
            html_report_path = os.path.join(self.scratch, 'gencomp_function_list.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(function_list,'col_header')
            html_report_txt.write("<h1>FUNCTIONS</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()

        if fun_list:
            report_string += "\nFEATURES WITH A FUNCTION\n"
            report_string += self.write_list_to_text_file(fun_list,'gencomp_functions2_tab.tsv',"\t")
            self.write_list_to_text_file(fun_list,'gencomp_functions2_csv.csv',",")
                                 
            html_report_path = os.path.join(self.scratch, 'gencomp_function_features.html')
            html_report_txt = open(html_report_path, "w", encoding="utf-8")
            html_table = self.make_HTML_table_from_list(fun_list,'col_header')
            html_report_txt.write("<h1>FEATURES WITH A FUNCTION</h1>")
            html_report_txt.write(html_table)
            html_report_txt.close()

        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}
                  
        #END genomecomp_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method genomecomp_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def domain_report(self, ctx, params):
        """
        :param params: instance of type "DomainReportParams" -> structure:
           parameter "input_ref" of type "domain_ref", parameter
           "evalue_cutoff" of Double, parameter "workspace_name" of String
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
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        domain_anno = self.dfu.get_objects({'object_refs': [input_ref]})
        domain_data = domain_anno['data'][0]['data']

        evalue_cutoff = float(params['evalue_cutoff'])

        report_list1 = []
        report_list2 = []
        
        create_features = CreateFeatureLists(self.config)
        (report_list1,report_list2,report_list3) = create_features.readDomainAnnList(domain_data, evalue_cutoff)
    
        report_string = ''
        resources = "<p>The following may be helpful</p><ul>"
        resources += "<li>COGS - https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/</li>"
        resources += "<li>Pfam - http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/</li>"
        resources += "<li>TIGRFAM - ftp://ftp.jcvi.org/pub/data/TIGRFAMs/ or https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/</li>"
        resources += "<li>NCBIfams - https://ftp.ncbi.nlm.nih.gov/hmm/9.0/</li>"
        resources += "<li>CDD (COGS, NCBI-curated, SMART, PRK) - ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd</li>"
        resources += "</ul>"
        if report_list1 or report_list2:
            if report_list1:
                report_string += self.write_list_to_text_file(report_list1,'domain_annotation_byGene.tsv',"\t")
                self.write_list_to_text_file(report_list1,'domain_annotation_byGene.csv',",")
                
                html_report_path = os.path.join(self.scratch, 'domain_annotation_byGene.html')
                html_report_txt = open(html_report_path, "w", encoding="utf-8")
                html_report_txt.write("<h1>LIST OF GENES AND THEIR DOMAINS</h1>")
                html_report_txt.write(resources)
                html_table = self.make_HTML_table_from_list(report_list1,'col_header')
                html_report_txt.write(html_table)
                html_report_txt.close()
                
            if report_list2:
                report_string += self.write_list_to_text_file(report_list2,'domain_annotation_byDomain.tsv',"\t")
                self.write_list_to_text_file(report_list2,'domain_annotation_byDomain.csv',",")
                        
                html_report_path = os.path.join(self.scratch, 'domain_annotation_byDomain.html')
                html_report_txt = open(html_report_path, "w", encoding="utf-8")
                html_report_txt.write("<h1>COUNTS PER DOMAIN</h1>")
                html_report_txt.write(resources)
                html_table = self.make_HTML_table_from_list(report_list2,'col_header')
                html_report_txt.write(html_table)
                html_report_txt.close()
                    
            if report_list3:
                report_string += self.write_list_to_text_file(report_list3,'domain_annotation_byCategory.tsv',"\t")
                self.write_list_to_text_file(report_list3,'domain_annotation_byCategory.csv',",")
                        
                html_report_path = os.path.join(self.scratch, 'domain_annotation_byCategory.html')
                html_report_txt = open(html_report_path, "w", encoding="utf-8")
                html_report_txt.write("<h1>COUNTS PER CATEGORY</h1>")
                html_report_txt.write(resources)
                html_table = self.make_HTML_table_from_list(report_list3,'col_header')
                html_report_txt.write(html_table)
                html_report_txt.close()

        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
        #END domain_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method domain_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def featseq_report(self, ctx, params):
        """
        :param params: instance of type "FeatSeqReportParams" -> structure:
           parameter "input_ref" of type "featseq_ref", parameter
           "workspace_name" of String
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
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.
        
        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']
        
        create_features = CreateFeatureLists(self.config)
        setseq = self.dfu.get_objects({'object_refs': [input_ref]})
        setseq_data = setseq['data'][0]['data']

        report_list = []
        seq_list = []
        dna_string = ""
 
        html_report_path = os.path.join(self.scratch, 'sequence_file.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")
        
        (header,desc_list, report_list,seq_list) = create_features.readFeatSeq(setseq_data)
        report_string = ''
        if desc_list:
            report_string += "DESCRIPTION\n"
            html_report_txt.write("<h1>DESCRIPTION</h1>")
            report_string += self.write_list_to_text_file(desc_list,'sequence_set_tab_desc.tsv',"\t")
            self.write_list_to_text_file(desc_list,'sequence_set_csv_desc.csv',",")
            html_table = self.make_HTML_table_from_list(desc_list,'row_header')
            html_report_txt.write(html_table)
        if report_list:
            report_string += "\n"+header+"\n"
            report_string += self.write_list_to_text_file(report_list,'sequence_set_tab_list.tsv',"\t")
            self.write_list_to_text_file(report_list,'sequence_set_csv_list.csv',",")
            html_table = self.make_HTML_table_from_list(report_list,'col_header')
            html_report_txt.write("<h1>"+header+"</h1>")
            html_report_txt.write(html_table)
        if seq_list:
            report_string += "\n"+header+"\n"
            dna_string += self.write_list_to_text_file(seq_list,'sequence_set_list.fasta',"\t")
            html_report_txt.write("<h1>"+header+"</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")

        html_report_txt.close()
        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        report_string += dna_string
        
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
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
           "workspace_name" of String
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
        create_features = CreateFeatureLists(self.config)
        
        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting ProteomeComparison Object Info. ')
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")

        # Step 1 - Parse/examine the parameters and catch any errors
        # It is important to check that parameters exist and are defined, and that nice error
        # messages are returned to users.

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        protcomp = self.dfu.get_objects({'object_refs': [input_ref]})
        protcomp_data = protcomp['data'][0]['data']
        report_list = []
        
        report_list1, report_list2 = create_features.readProtComp(protcomp_data)

        #Write to csv/tsv file. Create string for html and report output
        report_string = ''
        if report_list1:
            report_string += "List of Genomes\n"
            report_string += self.write_list_to_text_file(report_list1,'protcomp_tab_genome_list.tsv',"\t")
            self.write_list_to_text_file(report_list1,'protcomp_csv_genome_list.csv',",")
        if report_list2:
            report_string += "\nPairwise Comparison of Genomes\n"
            report_string += self.write_list_to_text_file(report_list2,'protcomp_tab_list.tsv',"\t")
            self.write_list_to_text_file(report_list2,'protcomp_csv_list.csv',",")
        
        html_report_path = os.path.join(self.scratch, 'protcomp_genomes.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")
        html_table = self.make_HTML_table_from_list(report_list1,'row_header')
        html_report_txt.write("<h1>LIST OF GENOMES</h1>")
        html_report_txt.write(html_table)
        html_table = self.make_HTML_table_from_list(report_list2,'col_header')
        html_report_txt.write("<h1>PAIRWISE GENOME COMPARISON</h1>")
        html_report_txt.write(html_table)
        html_report_txt.close()

        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)

        output = {'report_name': reported_output['name'],
                  'report_ref': reported_output['ref']}

        #END protcomp_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method protcomp_report return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def msa_report(self, ctx, params):
        """
        :param params: instance of type "MSAReportParams" -> structure:
           parameter "input_ref" of type "msa_ref", parameter
           "workspace_name" of String
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
        #BEGIN msa_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting Multiple Sequence Alignment Object Info. ')
        my_string = pformat(params)
        logging.info(f"Params:\n{my_string}")

        input_ref = params['input_ref']
        create_fasta = CreateFasta(self.config)
        
        html_report_path = os.path.join(self.scratch, 'MSA_file.html')
        html_report_txt = open(html_report_path, "w", encoding="utf-8")

        report_list = []
        msa_list = self.dfu.get_objects({'object_refs': [input_ref]})['data'][0]['data']
        
        name = "Generic MSA"
        if 'name' in msa_list:
            name = msa_list['name']
        desc = 'Unknown'
        if 'description' in msa_list:
            desc = msa_list['description']
        align_length = 0
        if 'alignment_length' in msa_list:
            align_length = msa_list['alignment_length']
        seq_type = 'Unknown'
        if 'sequence_type' in msa_list:
            seq_type = msa_list['sequence_type']
        row_order = []
        if 'row_order' in msa_list:
            row_order = msa_list['row_order']
        row_labels = {}
        if 'default_row_labels' in msa_list:
            row_labels = msa_list['default_row_labels']
        alignment = {}
        if 'alignment' in msa_list:
            alignment = msa_list['alignment']
        
        report_list = [["MSA Name","Description","Alignment Length","Sequence Type"]]
        report_list.extend([[name, desc, str(align_length),seq_type]])
        html_report_txt.write("<h1>MSA Overview</h1>")
        html_table = self.make_HTML_table_from_list(report_list,'col_header')
        html_report_txt.write(html_table)
        
        report_list.extend([[],['Row Lables'],["Row","Label"]])
        this_list = [["Row","Label"]]
        longlabel = 0
        for row in row_labels.keys():
            report_list.extend([[row,row_labels[row]]])
            this_list.extend([[row,row_labels[row]]])
            if len(row) > longlabel:
                longlabel = len(row)
            
        report_string = self.write_list_to_text_file(report_list,'MSA_tab_file.tsv',"\t")
        self.write_list_to_text_file(report_list,'MSA_csv_file.csv',",")
                        
        html_report_txt.write("<h1>Row Labels</h1>")
        html_table = self.make_HTML_table_from_list(this_list,'col_header')
        html_report_txt.write(html_table)
        
##      ALIGNMENT
        msa_list = [["CLUSTAL W multiple sequence alignment"],[]]
        colsz = 50
        start = 0
        while True:
            end = start + colsz
            if end > align_length:
                end = align_length
            for row in row_order:
                tmp_str = '{0:{3}} {1:50} {2:8}'.format(row,alignment[row][start:end],str(end),longlabel)
                #msa_list.append([[row,alignment[row][start:end],str(end)]])
                msa_list.append([tmp_str])
            msa_list.append([])
            start += colsz
            if start > align_length:
#                False
                break
                
        msa_string = self.write_list_to_text_file(msa_list,'MSA_alignment_file.aln',"\n")
                                
        html_table = self.make_HTML_table_from_list(msa_list,'col_header')
        html_report_txt.write("<h1>MSA Alignment</h1>")
        html_report_txt.write("<pre>" + msa_string + "</pre>")
        html_report_txt.close()
        
##      FASTA
        fasta_list = []
        for row in row_order:
            align = alignment[row].replace('-','')
            fasta_list.extend([[">" + row ]])
            fasta_list.extend(create_fasta.splitSequence(align))
        
        fasta_file_name = 'MSA_fasta_file.fna'
        if seq_type == 'protein':
            fasta_file_name = 'MSA_protein_file.faa'
        fasta_string = self.write_list_to_text_file(fasta_list,fasta_file_name,"\n")
                                    
        create_report = Report_creator(self.config)
        report_string = self.truncate_report_string(report_string)
        reported_output = create_report.create_report(token, params['workspace_name'],
                                           report_string, self.scratch)
                                           
        output = {'report_name': reported_output['name'],
                    'report_ref': reported_output['ref']}
        my_string = pformat(output)
        logging.info(f"Returning: {my_string}")
        #END msa_report

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method msa_report return value ' +
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
