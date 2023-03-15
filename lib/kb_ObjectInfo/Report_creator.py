import time
import os
import shutil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport

class Report_creator:
    def __init__(self, config):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = os.path.abspath(config['scratch'])

    def create_html_links(self, label,file,index_name,html_folder):
        html_string = ""
        output_link = {}
        with open('/kb/module/data/index_start.txt', 'r') as start_file:
            html_string = start_file.read()
            
        html_string += "        </div>    </div>    <div id=\"body\">\n"
        html_string += "        <iframe id=\"content\" "
        html_string += "style=\"width: 100%; border: none; \" src=\"" + file + "\"></iframe>\n    </div>"
        
        with open('/kb/module/data/index_end.txt', 'r') as end_file:
            html_string += end_file.read()

        with open(os.path.join(html_folder, index_name), 'w') as index_name:
            index_name.write(html_string)

        shock = self.dfu.file_to_shock({'file_path': html_folder,
                                        'make_handle': 0,
                                        'pack': 'zip'})
        desc = 'Open the text Report in a new window'
        output_link = {'shock_id': shock['shock_id'],
                                  'name': index_name,
                                  'label': label,
                                  'description': ''}
        return output_link
        
    # -----------------------------------------------------------------
    #    Create a Delimited Table version of the genes in a genome
    #


    def create_report(self, token, ws, report_string, read_file_path):
        # type: (object, object, object, object) -> object
        output_html_files = list()
        output_zip_files = list()
        first_file = ""
        html_string = ""
        html_count = 0
        with open('/kb/module/data/index_start.txt', 'r') as start_file:
            html_string = start_file.read()
        html_string += "        </div>    </div>    <div id=\"body\">\n"
        
        # Make HTML folder
        html_folder = os.path.join(read_file_path, 'html')
        if not os.path.exists(html_folder):
            os.mkdir(html_folder)

        for file in os.listdir(read_file_path):
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
            elif (file.endswith(".html")):
                # Move html into html folder
                shutil.move(os.path.join(read_file_path, file), os.path.join(html_folder, file))
                index_name = "index"+str(html_count)+".html"
                output_link = self.create_html_links(file,file,index_name,html_folder)
    
                if (first_file == ""):
                    first_file = file
                    html_string += "        <iframe id=\"content\" "
                    html_string += "style=\"width: 100%; border: none; \" src=\"" + first_file + "\"></iframe>\n    </div>"
                    
                html_count += 1

        with open('/kb/module/data/index_end.txt', 'r') as end_file:
            html_string += end_file.read()

        with open(os.path.join(html_folder, "index.html"), 'w') as index_file:
            index_file.write(html_string)

        shock = self.dfu.file_to_shock({'file_path': html_folder,
                                        'make_handle': 0,
                                        'pack': 'zip'})
        desc = 'Open the text Report in a new window'
        output_html_files.append({'shock_id': shock['shock_id'],
                                  'name': 'index.html',
                                  'label': 'HTML Link',
                                  'description': ''})
                                  
        short_report = report_string[0:1000]
        report_params = {
            'objects_created': [],
            'message': "PREVIEW of the first 1000 characters of the output. Tab and comma-delimited output may look funny on the screen because they are intended to be read by a computer. For a 'pretty' output, use the data panel and the objects data viewer, either in the narrative or the landing page.\nFiles in the 'Links' section below will open in a new browser window. Links in the 'Files' section below will download the output when you click on the name.\n\n" + short_report,
            'direct_html': '',
            'direct_html_link_index': None,
            'file_links': output_zip_files,
            'html_links': output_html_files,
            'workspace_name': ws,
            'report_object_name': 'kb_ObjectInfo_report'
        }
#       Keep this code in for a little while. I might change my mind. Was > 1000000
        if len(report_string) < 1:
            report_params = {
                'objects_created': [],
                'message': '',
                'direct_html': "Links in the 'HTML' section below will open in a new window. Links in the 'Files' section below will download the output when you click on the name.",
                'direct_html_link_index': None,
                'file_links': output_zip_files,
                'html_links': output_html_files,
                'workspace_name': ws,
                'report_object_name': 'kb_ObjectInfo_report'
            }
        kbase_report_client = KBaseReport(self.callback_url, token=token)
        output = kbase_report_client.create_extended_report(report_params)
        return output
