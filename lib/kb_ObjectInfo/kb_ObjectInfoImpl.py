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
    VERSION = "1.2.0"
    GIT_URL = "https://github.com/kbaseapps/kb_ObjectInfo"
    GIT_COMMIT_HASH = "0e8d1c51dd1f4d33c1803ad230d13421b3ca47fc"

    #BEGIN_CLASS_HEADER

#   rpt_list is a list of the rows that will form the output file
#   rpt_string is used for the PREVIEW for the user and for the html version of the file
#   htmltable is the html table version of rpt_list/rpt_string
#   report_path is the physical location of the output files
#   report_txt is an internal name used for referencing the report_path
#   rpt_writer is a writer object responsible for converting user data to delimited strings
#   report_format is tsv, csv, or something else selected by the user
#   rpt_delimiter is the delimiter used in the file. "\t" for tsv and ',' for csv, etc.

    def make_HTML(self,rpt_list,style):
        table = "<table style=\"border: 1px solid black;\">\n"

        # Create the table's row data
        for i, row in enumerate(rpt_list):
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

    def write_to_file(self,rpt_list,rpt_path,rpt_delimiter):
        rpt_string = ''
        report_path = os.path.join(self.scratch, rpt_path)
        with open(report_path, mode='w') as report_txt:
            rpt_writer = csv.writer(report_txt, delimiter=rpt_delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
            if rpt_delimiter == ',':
                    rpt_writer = csv.writer(report_txt, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL, dialect='excel')
            for rpt in rpt_list:
                rpt_writer.writerow(rpt)
                rpt_string += rpt_delimiter.join(rpt) + "\n"
            report_txt.close()
            
        return rpt_string
       
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
        
        car = CreateAssemblyReport(self.config)
        assembly = self.dfu.get_objects({'object_refs': [input_ref]})
        rpt_list = []
        this_list = []
        
        html_report_path = os.path.join(self.scratch, 'assembly_metadata_file.html')
        html_report_txt = open(html_report_path, "w")
        
        (header,this_list) = car.assembly_overview(assembly)
        
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)
        
        (header,this_list) = car.assembly_metadata(assembly)
        
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)
        
        (header,this_list) = car.assembly_dnabases(assembly)
                
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)
        
        (header,this_list) = car.assembly_contigs(assembly)
                
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)
        html_report_txt.close()

        rpt_string = self.write_to_file(rpt_list,'assembly_meta_tab_file.tsv',"\t")
        self.write_to_file(rpt_list,'assembly_meta_csv_file.csv',",")
                                          
        if showContigs:
            (header) = car.assembly_contigs(assembly)
            rpt_list.extend([[],[header]])
        
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

    def assemblyset_report(self, ctx, params):
        """
        :param params: instance of type "AssemblySetReportParams" ->
           structure: parameter "input_ref" of type "assemblyset_ref",
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
        #BEGIN assemblyset_report
        token = ctx['token']

        # Logging statements to stdout/stderr are captured and available as the App log
        logging.info('Starting AssemblySet Object Info. ')
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")

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

        assemblyset = self.dfu.get_objects({'object_refs': [input_ref]})
        assemblylist = assemblyset['data'][0]['data']['items']
        num_assem = len(assemblylist)
        car = CreateAssemblyReport(self.config)
        
        assemref = []
        
        for i in range(num_assem):
            assemref.append(assemblylist[i]['ref'])
        
        html_report_path = os.path.join(self.scratch, 'assembly_metadata_file.html')
        html_report_txt = open(html_report_path, "w")
        
        this_list = []
        rpt_list = []
        obj_list = self.dfu.get_objects({'object_refs': assemref})
        
        (header,this_list) = car.assembly_overview(obj_list)
        
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[header]])
        rpt_list.extend(this_list)

        (header,this_list) = car.assembly_metadata(obj_list)
        
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)
        
        (header,this_list) = car.assembly_dnabases(obj_list)
                
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)

        (header,this_list) = car.assembly_contigs(obj_list)
                
        htmltable = self.make_HTML(this_list,'col_header')
        html_report_txt.write("<h1>"+header+"</h1>")
        html_report_txt.write(htmltable)
        rpt_list.extend([[],[header]])
        rpt_list.extend(this_list)
        html_report_txt.close()
                                  
        if showContigs:
            (header) = car.assembly_contigs(obj_list)
            rpt_list.extend([[],[header]])
        
        output = {'report_name': reported_output['name'],
                    'report_ref': reported_output['ref']}
        mystr = pformat(output)
        logging.info(f"Returning: {mystr}")
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
           "workspace_name" of String, parameter "listCoding" of type
           "boolean" (A boolean. 0 = false, other = true.), parameter
           "listGFF" of type "boolean" (A boolean. 0 = false, other = true.),
           parameter "FastaAA" of type "boolean" (A boolean. 0 = false, other
           = true.), parameter "FastamRNA" of type "boolean" (A boolean. 0 =
           false, other = true.), parameter "showDNA" of type "boolean" (A
           boolean. 0 = false, other = true.)
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

        genome = self.dfu.get_objects({'object_refs': [input_ref]})
        genome_data = genome['data'][0]['data']

        rpt_list = []
        rpt_string = ''
        report_path = ''
        
        if params['listCoding']:
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.delimitedTable(genome_data, 'features')
            if rpt_list:
                rpt_string += self.write_to_file(rpt_list,'genome_tab_file.tsv',"\t")
                self.write_to_file(rpt_list,'genome_csv_file.csv',",")
                
            html_report_path = os.path.join(self.scratch, 'genome_features.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(rpt_list,'col_header')
            html_report_txt.write("<h1>PROTEIN CODING FEATURES</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()
                
        if params['listGFF']:
            cf = CreateFeatureLists(self.config)
            rpt_list = cf.gff3(genome_data, 'features')
            rpt_string += self.write_to_file(rpt_list,'genome_file.gff',"\t")
                            
            html_report_path = os.path.join(self.scratch, 'genome_GFF.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(rpt_list,'col_header')
            html_report_txt.write("<h1>GFF OUTPUT</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()
            
        if params['FastaAA']:
            cf = CreateFasta(self.config)
#           Before version 9 genomes, the cdss didn't exist
            if genome_data['cdss']:
                rpt_list = cf.create_Fasta_from_features(genome_data['cdss'])
            else:
                rpt_list = cf.create_Fasta_from_features(genome_data['features'])
            
            dna_string = self.write_to_file(rpt_list,'genome_file.faa',"\n")
            rpt_string += dna_string[0:200] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
            html_report_path = os.path.join(self.scratch, 'genome_protein_AA.html')
            html_report_txt = open(html_report_path, "w")
            html_report_txt.write("<h1>PROTEIN CODING AMINO ACIDS</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")
            html_report_txt.close()
            
        if params['FastamRNA']:
            cf = CreateFasta(self.config)
#           Before version 9 genomes, the cdss didn't exist
            if genome_data['cdss']:
                rpt_list = cf.create_Fasta_from_mRNA(genome_data['cdss'])
            else:
                rpt_list = cf.create_Fasta_from_mRNA(genome_data['features'])
                
            dna_string = self.write_to_file(rpt_list,'genome_mRNA_file.fna',"\n")
            rpt_string += dna_string[0:200] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
            html_report_path = os.path.join(self.scratch, 'genome_protein_mRNA.html')
            html_report_txt = open(html_report_path, "w")
            html_report_txt.write("<h1>PROTEIN CODING mRNA</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")
            html_report_txt.close()
        
        if params['showDNA']:
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
                
            dna_string = self.write_to_file(rpt_list,'genome_dna_file.fna',"\n")
            rpt_string += dna_string[0:200] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
            html_report_path = os.path.join(self.scratch, 'genome_protein_mRNA.html')
            html_report_txt = open(html_report_path, "w")
            html_report_txt.write("<h1>DNA</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")
            html_report_txt.close()

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
           "workspace_name" of String, parameter "showGenomes" of type
           "boolean" (A boolean. 0 = false, other = true.), parameter
           "showDNA" of type "boolean" (A boolean. 0 = false, other = true.)
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

        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']

        genomeset = self.dfu.get_objects({'object_refs': [input_ref]})
        genome_name = genomeset['data'][0]['info'][1]
        genomeset_data = genomeset['data'][0]['data']

        rpt_list = []
        multi_fasta = []
        rpt_string = ''
        
        if params['showGenomes']:
            gsr = CreateMultiGenomeReport(self.config)
            rpt_list = gsr.readGenomeSet(genome_name, genomeset_data)
            rpt_string += self.write_to_file(rpt_list,'genomeset_tab_file.tsv',"\t")
            self.write_to_file(rpt_list,'genomeset_cvs_file.csv',",")
            htmltable = self.make_HTML(rpt_list,'col_header')
            
            html_report_path = os.path.join(self.scratch, 'genomeset_comparison.html')
            html_report_txt = open(html_report_path, "w")
            html_report_txt.write("<h1>GENOME COMPARISON</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()
        
        if params['showDNA']:
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
            rpt_string += self.write_to_file(rpt_list,'genomeset_DNA_meta_file.tsv',"\t")
            self.write_to_file(rpt_list,'genomeset_DNA_meta_file.csv',",")
            
#           The fasta not included in the HTML due to length
            html_report_path = os.path.join(self.scratch, 'genomeset_DNA_meta.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(rpt_list,'col_header')
            html_report_txt.write("<h1>GENOME DOWNLOAD METADATA</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()

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
        mystr = pformat(params)
        logging.info(f"Params:\n{mystr}")
        
        if 'workspace_name' not in params:
            raise ValueError('Parameter workspace_name is not set in input arguments')
        workspace_name = params['workspace_name']
        if 'input_ref' not in params:
            raise ValueError('Parameter input_ref is not set in input arguments')
        input_ref = params['input_ref']
        
        cf = CreateFeatureLists(self.config)
        data_file_cli = DataFileUtil(self.callback_url)
        gencomp = self.dfu.get_objects({'object_refs': [input_ref]})
        gencomp_data = gencomp['data'][0]['data']
        
        rpt_list1 = []
        rpt_list2 = []
        cf = CreateFeatureLists(self.config)
        rpt_string = ''
        html_report_path = os.path.join(self.scratch, 'pangenome_genome_comparison.html')
        html_report_txt = open(html_report_path, "w")
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
            rpt_string += "\nOVERVIEW\n"
            rpt_string += self.write_to_file(overview_list,'gencomp_overview_tab.tsv',"\t")
            self.write_to_file(overview_list,'gencomp_overview_csv.csv',",")
                
            htmltable = self.make_HTML(overview_list,'row_header')
            html_report_txt.write("<h1>OVERVIEW</h1>")
            html_report_txt.write(htmltable)
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
            rpt_string += "\nGENOMES\n"
            rpt_string += self.write_to_file(genome_list,'gencomp_genomes_tab.tsv',"\t")
            self.write_to_file(genome_list,'gencomp_genomes_csv.csv',",")
   
            htmltable = self.make_HTML(genome_list,'col_header')
            html_report_txt.write("<h1>GENOMES</h1>")
            html_report_txt.write(htmltable)
                
        if sim_list:
            rpt_string += "\nSIMILARITY LIST\n"
            rpt_string += self.write_to_file(sim_list,'gencomp_genomes_sim_tab.tsv',"\t")
            self.write_to_file(sim_list,'gencomp_genomes_sim_csv.csv',",")
   
            htmltable = self.make_HTML(sim_list,'col_header')
            html_report_txt.write("<h1>SIMILARITY LIST</h1>")
            html_report_txt.write(htmltable)
            
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
            rpt_string += "\nFAMILIES\n"
            rpt_string += self.write_to_file(family_list,'gencomp_familes_tab.tsv',"\t")
            self.write_to_file(family_list,'gencomp_families_csv.csv',",")
                        
            html_report_path = os.path.join(self.scratch, 'gencomp_family_list.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(family_list,'col_header')
            html_report_txt.write("<h1>FAMILIES</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()
                
        if fam_list:
            rpt_string += "\nGENES IN FAMILIES\n"
            rpt_string += self.write_to_file(fam_list,'gencomp_familes2_tab.tsv',"\t")
            self.write_to_file(fam_list,'gencomp_families2_csv.csv',",")
                                
            html_report_path = os.path.join(self.scratch, 'gencomp_family_genes.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(fam_list,'col_header')
            html_report_txt.write("<h1>GENES IN FAMILIES</h1>")
            html_report_txt.write(htmltable)
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
            rpt_string += "\nFUNCTIONS\n"
            rpt_string += self.write_to_file(function_list,'gencomp_functions_tab.tsv',"\t")
            self.write_to_file(function_list,'gencomp_functions_csv.csv',",")
                                 
            html_report_path = os.path.join(self.scratch, 'gencomp_function_list.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(function_list,'col_header')
            html_report_txt.write("<h1>FUNCTIONS</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()

        if fun_list:
            rpt_string += "\nFEATURES WITH A FUNCTION\n"
            rpt_string += self.write_to_file(fun_list,'gencomp_functions2_tab.tsv',"\t")
            self.write_to_file(fun_list,'gencomp_functions2_csv.csv',",")
                                 
            html_report_path = os.path.join(self.scratch, 'gencomp_function_features.html')
            html_report_txt = open(html_report_path, "w")
            htmltable = self.make_HTML(fun_list,'col_header')
            html_report_txt.write("<h1>FEATURES WITH A FUNCTION</h1>")
            html_report_txt.write(htmltable)
            html_report_txt.close()

        cr = Report_creator(self.config)

        reported_output = cr.create_report(token, params['workspace_name'],
                                           rpt_string, self.scratch)

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

        domain_anno = self.dfu.get_objects({'object_refs': [input_ref]})
        domain_data = domain_anno['data'][0]['data']

        evalue_cutoff = float(params['evalue_cutoff'])

        rpt_list1 = []
        rpt_list2 = []
        
        cf = CreateFeatureLists(self.config)
        rpt_list1 = cf.readDomainAnnList(domain_data, evalue_cutoff)
        rpt_list2 = cf.readDomainAnnCount(domain_data, evalue_cutoff)
    
        rpt_string = ''
        
        if rpt_list1 or rpt_list2:
            if rpt_list2:
                rpt_string += self.write_to_file(rpt_list2,'domain_annotation_byDomain.tsv',"\t")
                self.write_to_file(rpt_list2,'domain_annotation_byDomain.csv',",")
                        
                html_report_path = os.path.join(self.scratch, 'domain_annotation_byDomain.html')
                html_report_txt = open(html_report_path, "w")
                html_report_txt.write("<h1>COUNTS PER DOMAIN</h1>")
                htmltable = self.make_HTML(rpt_list2,'col_header')
                html_report_txt.write(htmltable)
                html_report_txt.close()

            if rpt_list1:
                rpt_string += self.write_to_file(rpt_list1,'domain_annotation_byGene.tsv',"\t")
                self.write_to_file(rpt_list1,'domain_annotation_byGene.csv',",")
                
                html_report_path = os.path.join(self.scratch, 'domain_annotation_byGene.html')
                html_report_txt = open(html_report_path, "w")
                html_report_txt.write("<h1>LIST OF GENES AND THEIR DOMAINS</h1>")
                
                htmltable = self.make_HTML(rpt_list1,'col_header')
                html_report_txt.write(htmltable)
                html_report_txt.close()
                
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
        setseq = self.dfu.get_objects({'object_refs': [input_ref]})
        setseq_data = setseq['data'][0]['data']

        rpt_list = []
        seq_list = []
        dna_string = ""
 
        html_report_path = os.path.join(self.scratch, 'sequence_file.html')
        html_report_txt = open(html_report_path, "w")
        
        (header,desc_list, rpt_list,seq_list) = cf.readFeatSeq(setseq_data)
        rpt_string = ''
        if desc_list:
            rpt_string += "DESCRIPTION\n"
            html_report_txt.write("<h1>DESCRIPTION</h1>")
            rpt_string += self.write_to_file(desc_list,'sequence_set_tab_desc.tsv',"\t")
            self.write_to_file(desc_list,'sequence_set_csv_desc.csv',",")
            htmltable = self.make_HTML(desc_list,'row_header')
            html_report_txt.write(htmltable)
        if rpt_list:
            rpt_string += "\n"+header+"\n"
            rpt_string += self.write_to_file(rpt_list,'sequence_set_tab_list.tsv',"\t")
            self.write_to_file(rpt_list,'sequence_set_csv_list.csv',",")
            htmltable = self.make_HTML(rpt_list,'col_header')
            html_report_txt.write("<h1>"+header+"</h1>")
            html_report_txt.write(htmltable)
        if seq_list:
            rpt_string += "\n"+header+"\n"
            dna_string += self.write_to_file(seq_list,'sequence_set_list.fasta',"\t")
            html_report_txt.write("<h1>"+header+"</h1>")
            html_report_txt.write("<pre>" + dna_string + "</pre>")

        html_report_txt.close()
        cr = Report_creator(self.config)
        rpt_string += dna_string
        
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
        cf = CreateFeatureLists(self.config)
        
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

        protcomp = self.dfu.get_objects({'object_refs': [input_ref]})
        protcomp_data = protcomp['data'][0]['data']
        rpt_list = []
        
        rpt_list1, rpt_list2 = cf.readProtComp(protcomp_data)

        #Write to csv/tsv file. Create string for html and report output
        rpt_string = ''
        if rpt_list1:
            rpt_string += "List of Genomes\n"
            rpt_string += self.write_to_file(rpt_list1,'protcomp_tab_genome_list.tsv',"\t")
            self.write_to_file(rpt_list1,'protcomp_csv_genome_list.csv',",")
        if rpt_list2:
            rpt_string += "\nPairwise Comparison of Genomes\n"
            rpt_string += self.write_to_file(rpt_list2,'protcomp_tab_list.tsv',"\t")
            self.write_to_file(rpt_list2,'protcomp_csv_list.csv',",")
        
        html_report_path = os.path.join(self.scratch, 'protcomp_genomes.html')
        html_report_txt = open(html_report_path, "w")
        htmltable = self.make_HTML(rpt_list1,'row_header')
        html_report_txt.write("<h1>LIST OF GENOMES</h1>")
        html_report_txt.write(htmltable)
        htmltable = self.make_HTML(rpt_list2,'col_header')
        html_report_txt.write("<h1>PAIRWISE GENOME COMPARISON</h1>")
        html_report_txt.write(htmltable)
        html_report_txt.close()

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
