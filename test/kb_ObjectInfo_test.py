# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import requests

from os import environ
try:
    from configparser import ConfigParser  # py3
except:
    from configparser import ConfigParser  # py2

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from kb_ObjectInfo.kb_ObjectInfoImpl import kb_ObjectInfo
from kb_ObjectInfo.kb_ObjectInfoServer import MethodContext
from kb_ObjectInfo.authclient import KBaseAuth as _KBaseAuth

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

class kb_ObjectInfoTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_ObjectInfo'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_ObjectInfo',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_ObjectInfo(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        suffix = int(time.time() * 1000)
        cls.wsName = "test_Reports_" + str(suffix)
        cls.ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})

#       Prepare the Assembly File
#       Set the name
#        cls.test_filename = '92477.assembled.fna'
        cls.test_filename = 'test.fna'
#       Set the path to file in scratch
        cls.test_path = os.path.join(cls.scratch, cls.test_filename)
#       Copy from local to scratch
        shutil.copy(os.path.join("data", cls.test_filename), cls.test_path)

#       Prepare a minimal Genome from gbff File
        cls.genbank_file_name = 'minimal.gbff'
#       Set the path to file in scratch
        cls.genbank_file_path = os.path.join(cls.scratch, cls.genbank_file_name)
#       Copy from local to scratch
        shutil.copy(os.path.join('data', cls.genbank_file_name), cls.genbank_file_path)
#       Convert gbff to Genome Object and save the genome_ref
        cls.gfu = GenomeFileUtil(cls.callback_url)
        genome_object_name = 'minimal_Genome'
        cls.minimal_genome_ref = cls.gfu.genbank_to_genome({'file': {'path': cls.genbank_file_path},
                                                    'workspace_name': cls.ws_info[1],
                                                    'genome_name': 'minimal_test_genome',
                                                    "generate_ids_if_needed": 1,
                                                    "generate_missing_genes": 1
                                                    })['genome_ref']
                                                    
#       Prepare the Genome from gbff File
        cls.genbank_file_name = 'Carsonella_ruddii_HT_isolate_Thao2000.gbff'
#       Set the path to file in scratch
        cls.genbank_file_path = os.path.join(cls.scratch, cls.genbank_file_name)
#       Copy from local to scratch
        shutil.copy(os.path.join('data', cls.genbank_file_name), cls.genbank_file_path)
#       Convert gbff to Genome Object and save the genome_ref
        cls.gfu = GenomeFileUtil(cls.callback_url)
        genome_object_name = 'test_Genome'
        genome_ref = "39201/10/1"
        genome_ref = "40843/4/1"
        cls.genome_ref = cls.gfu.genbank_to_genome({'file': {'path': cls.genbank_file_path},
                                                    'workspace_name': cls.ws_info[1],
                                                    'genome_name': 'first_test_genome',
                                                    "generate_ids_if_needed": 1,
                                                    "generate_missing_genes": 1
                                                    })['genome_ref']
#       Prepare the GenomeSet
        cls.genbank_file_name = 'Carsonella_ruddii_Thao2000.gbk_genome.gbff'
#       Set the path to file in scratch
        cls.genbank_file_path = os.path.join(cls.scratch, cls.genbank_file_name)
#       Copy from local to scratch
        shutil.copy(os.path.join('data', cls.genbank_file_name), cls.genbank_file_path)
#       Convert gbff to Genome Object and save the genome_ref
        cls.genome_ref2 = cls.gfu.genbank_to_genome({'file': {'path': cls.genbank_file_path},
                                                     'workspace_name': cls.ws_info[1],
                                                     'genome_name': 'second_test_genome',
                                                     "generate_ids_if_needed": 1,
                                                     "generate_missing_genes": 1
                                                      })['genome_ref']

        #       Prepare the Domain Annotation File
        domain_data_file = 'Carsonella.Domains.json'
        cls.domain_file = os.path.join(cls.scratch, domain_data_file)
        shutil.copy(os.path.join('data', domain_data_file), cls.domain_file)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_ObjectInfo_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def get_fasta_file(self, filename, obj_name):
        assemblyUtil = AssemblyUtil(self.callback_url)
        assembly_ref = assemblyUtil.save_assembly_from_fasta({'file': {'path': filename},
                                                              'workspace_name': self.getWsName(),
                                                              'assembly_name': obj_name
                                                              })
        return assembly_ref

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def getDomainInfo(self, domain_name):
        domain_name = "COG0486"
        lib_i = 0
        if hasattr(self.__class__, 'domainInfo_list'):
            try:
                info = self.__class__.domainInfo_list[lib_i]
                name = self.__class__.domainName_list[lib_i]
                if info != None:
                    if name != domain_name:
                        self.__class__.domainInfo_list[lib_i] = None
                        self.__class__.domainName_list[lib_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform json to kbase DomainAnnotation object and upload to ws
        # create object
        with open(self.domain_file, 'r') as domain_fh:
            domain_obj = json.load(domain_fh)

#        genome_ref = self.gfu.genbank_to_genome({'file': {'path': self.genbank_file_path},
#                                                    'workspace_name': self.getWsName(),
#                                                    'genome_name': 'minimal_test_genome'
#                                                    })['genome_ref']

        domain_obj['used_dms_ref'] = 'KBasePublicGeneDomains/All'
        domain_obj['genome_ref'] = self.genome_ref

        provenance = [{}]
        new_obj_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(),
            'objects': [
                {
                    'type': 'KBaseGeneFamilies.DomainAnnotation',
                    'data': domain_obj,
                    'name': 'test_DOMAINS',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'domainInfo_list'):
            self.__class__.domainInfo_list = []
            self.__class__.domainName_list = []
        for i in range(lib_i+1):
            try:
                assigned = self.__class__.domainInfo_list[i]
            except:
                self.__class__.domainInfo_list.append(None)
                self.__class__.domainName_list.append(None)

        self.__class__.domainInfo_list[lib_i] = new_obj_info
        self.__class__.domainName_list[lib_i] = domain_name

        return str(new_obj_info[6]) + "/" + str(new_obj_info[0]) + "/" + str(new_obj_info[4])

    def getGenomeSet(self):

        # build GenomeSet obj
        testGS = {
            'description': 'two genomes',
            'elements': dict()
        }
        testGS['elements']['Carsonella 1'] = {'ref': self.genome_ref}
        testGS['elements']['Carsonella 2'] = {'ref': self.genome_ref2}
        print ("OLD:",self.genome_ref,'  ',self.genome_ref2)
        new_obj_info = self.getWsClient().save_objects({'workspace': self.ws_info[1],
                                                        'objects': [
                                                        {
                                                            'type': 'KBaseSearch.GenomeSet',
                                                            'data': testGS,
                                                            'name': 'test_genomeset',
                                                            'meta':  {},
                                                            'provenance': [
                                                                {
                                                                    'service': 'kb_phylogenomics',
                                                                    'method': 'test_view_fxn_profile'
                                                                }
                                                            ]
                                                        }]
                                                        })[0]
#        print ("NEW:", new_obj_info)
#
        return str(new_obj_info[6]) + "/" + str(new_obj_info[0]) + "/" + str(new_obj_info[4])

    def mytest_assembly_metadata(self):

        assembly_ref = self.get_fasta_file(self.test_path,
                                           'TestAssembly3')
        ret = self.getImpl().assembly_metadata_report(self.getContext(),
                                                      {'workspace_name': self.ws_info[1],
                                                       'assembly_input_ref': assembly_ref,
                                                       'showContigs': 1
                                                       })
        # Validate the returned data
        print(("ASSEMBLY METADATA RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def test_genome_tab(self):
        genome_ref = '40843/4/1'
        ret = self.getImpl().genome_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'genome_input_ref': self.genome_ref,
                                            'report_format': 'tab'
                                            })
        print(("GENOME TAB RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genome_gff(self):
        genome_ref = '40843/4/1'
        ret = self.getImpl().genome_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'genome_input_ref': self.genome_ref,
                                            'report_format': 'gff'
                                            })
        print(("GENOME GFF RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genome_fasta(self):
        ret = self.getImpl().genome_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'genome_input_ref': self.minimal_genome_ref,
                                            'report_format': 'fasta'
                                            })
        print(("GENOME FASTA RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genome_mrna(self):
        genome_ref = '40843/4/1'
        ret = self.getImpl().genome_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'genome_input_ref': self.minimal_genome_ref,
                                            'report_format': 'mRNA'
                                            })
        print(("GENOME MRNA RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genome_DNA(self):
        genome_ref = '40843/4/1'
        ret = self.getImpl().genome_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'genome_input_ref': self.minimal_genome_ref,
                                            'report_format': 'DNA'
                                            })
        print(("GENOME DNA RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_domain_annotation(self):
        domain_ref = self.getDomainInfo('test_domain')
        ret = self.getImpl().domain_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'evalue_cutoff': '1e-20',
                                            'domain_annotation_input_ref': domain_ref,
                                            'report_format': 'tab'
                                            })
        print(("DOMAIN ANNOTATION RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genomeset_meta(self):
        genomeset_ref = self.getGenomeSet()
        ret = self.getImpl().genomeset_report(self.getContext(),
                                              {'workspace_name': self.ws_info[1],
                                               'genomeset_input_ref': genomeset_ref,
                                               'report_format': 'meta'
                                               })
        print(("GENOME SET META RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genomeset_list(self):
        genomeset_ref = self.getGenomeSet()
        ret = self.getImpl().genomeset_report(self.getContext(),
                                              {'workspace_name': self.ws_info[1],
                                               'genomeset_input_ref': genomeset_ref,
                                               'report_format': 'list'
                                               })
        print(("GENOME SET LIST RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genomeset_tab(self):
        genomeset_ref = self.getGenomeSet()
        ret = self.getImpl().genomeset_report(self.getContext(),
                                              {'workspace_name': self.ws_info[1],
                                               'genomeset_input_ref': genomeset_ref,
                                               'report_format': 'tab'
                                               })
        print(("GENOME SET TAB RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genomeset_csv(self):
        genomeset_ref = self.getGenomeSet()
        ret = self.getImpl().genomeset_report(self.getContext(),
                                              {'workspace_name': self.ws_info[1],
                                               'genomeset_input_ref': genomeset_ref,
                                               'report_format': 'csv'
                                               })
        print(("GENOME SET CSV RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_genomeset_fasta(self):
        genomeset_ref = self.getGenomeSet()
        ret = self.getImpl().genomeset_report(self.getContext(),
                                              {'workspace_name': self.ws_info[1],
                                               'genomeset_input_ref': genomeset_ref,
                                               'report_format': 'fasta'
                                               })
        print(("GENOME SET FASTA RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_featureSet_ordered(self):
        featset_ref = '27092/14/1'
        ret = self.getImpl().featseq_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'feature_sequence_input_ref': featset_ref,
                                            'report_format': 'csv'
                                            })
        print(("FEATURE SET ORDERED RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass

    def mytest_featureSet_unordered(self):
        featset_ref = '20563/36/1'
        ret = self.getImpl().featseq_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'feature_sequence_input_ref': featset_ref,
                                            'report_format': 'csv'
                                            })
        print(("FEATURE SET UNORDERED RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass
      
    def mytest_ProtComp(self):
        protcomp_ref = '29939/15/1'
        ret = self.getImpl().protcomp_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'protcomp_input_ref': protcomp_ref,
                                            'report_format': 'csv'
                                            })
        print(("PROTEOME COMPARISON RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass
    
    def mytest_sequenceSet(self):
        featset_ref = '27092/23/1'
        ret = self.getImpl().featseq_report(self.getContext(),
                                           {'workspace_name': self.ws_info[1],
                                            'feature_sequence_input_ref': featset_ref,
                                            'report_format': 'tab'
                                            })
        print(("SEQUENCE SET RETURNED", ret))
        self.assertIn('report_name', ret[0])
        self.assertIn('report_ref', ret[0])
        pass
    
    
