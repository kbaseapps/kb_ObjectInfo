import time
import os
import shutil
import logging
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    logging.info(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class Report_creator:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = os.path.abspath(config['scratch'])

    # -----------------------------------------------------------------
    #    Create a Delimited Table version of the genes in a genome
    #


    def create_report(self, token, ws, report_string, read_file_path):
        # type: (object, object, object, object) -> object
        output_html_files = list()
        output_zip_files = list()

        html_string = ""
        html_count = 0
        html_file_list = {}
        # Make HTML folder
        html_folder = os.path.join(read_file_path, 'html')
        if not os.path.exists(html_folder):
            os.mkdir(html_folder)
        listdir = os.listdir(read_file_path)
        listdir.sort()
        
        for file in listdir:
            label = ".".join(file.split(".")[1:])
            if (file.endswith(".zip")):
                desc = 'Zip file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".faa")):
                desc = 'FASTA file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".fna")):
                desc = 'FASTA file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".txt")):
                desc = 'Text file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".tab") or file.endswith(".tsv")):
                desc = 'Tab-delimited text file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".csv")):
                desc = 'Comma-delimited text file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".gff")):
                desc = 'GFF3 text file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
            elif (file.endswith(".aln")):
                desc = 'MSA text file generated for output '
                output_zip_files.append({'path': os.path.join(read_file_path, file),
                                         'name': file,
                                         'label': label,
                                         'description': desc})
                                    
            elif (file.endswith(".html")):
                # Move html into html folder
                shutil.move(os.path.join(read_file_path, file), os.path.join(html_folder, file))

                # Create an index file -  header lines
                with open('/kb/module/data/index_start.txt', 'r') as start_file:
                    html_string = start_file.read()
                    
                html_string += "        </div>    </div>    <div id=\"body\">\n"
                html_string += "        <iframe id=\"content\" "
                html_string += "style=\"width: 100%; border: none; \" src=\"" + file + "\"></iframe>\n    </div>"

                # Add the closing lines
                with open('/kb/module/data/index_end.txt', 'r') as end_file:
                    html_string += end_file.read()
                
                # Write the index html file to the directory
                file_name = "index"+str(html_count)+".html"
                with open(os.path.join(html_folder, file_name), 'w') as index_file:
                    index_file.write(html_string)
            
                html_file_list[file_name] = file
                html_count += 1

        # Get an ID for the folder in shock
        shock = self.dfu.file_to_shock({'file_path': html_folder,
                                        'make_handle': 0,
                                        'pack': 'zip'})
                
        # Create links for all of the html index files
        desc = 'Open the text Report in a new window'
        for name in html_file_list.keys():
            output_html_files.append({'shock_id': shock['shock_id'],
                                    'name': name,
                                    'label': 'HTML Link '+html_file_list[name],
                                    'description': ''})
                                        
        short_report = report_string
        
        if len(report_string) > 1000:
            short_report = report_string[0:1000] + ".... <b>Summary has been truncated. Use the links or files for full output.</b>\n"
            
        report_params = {
            'objects_created': [],
            'message': "PREVIEW of the first 1000 characters of the output. Tab and comma-delimited output may look funny on the screen because they are intended to be read by a computer. For a 'pretty' display, use the Links below to get a full output in a new browser window. Really long output may bog down some browsers. In this case, download the file and view locally. The Files below will download the output when you click on the name. \n\n" + short_report,
            'direct_html': '',
            'direct_html_link_index': None,
            'file_links': output_zip_files,
            'html_links': output_html_files,
            'workspace_name': ws,
            'report_object_name': 'kb_ObjectInfo_report'
        }
        
        kbase_report_client = KBaseReport(self.callback_url, token=token)
        output = kbase_report_client.create_extended_report(report_params)
        return output
