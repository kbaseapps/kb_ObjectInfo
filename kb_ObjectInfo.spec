/*
A KBase module: kb_ObjectInfo
*/

module kb_ObjectInfo {
    /* 
        A 'typedef' allows you to provide a more specific name for
        a type.  Built-in primitive types include 'string', 'int',
        'float'.  Here we define a type named assembly_ref to indicate
        a string that should be set to a KBase ID reference to an
        Assembly data object.
        
        KIDL Specification, KBase SDK 1.2.0 documentation
        (https://kbase.github.io/kb_sdk_docs/references/KIDL_spec.html)
    */

    /* A boolean. 0 = false, other = true. */
    typedef int boolean;

    typedef string assembly_ref;
    typedef string genome_ref;
    typedef string genomeset_ref;
    typedef string domain_ref;
    typedef string tree_ref;
    typedef string featseq_ref;
    typedef string protcomp_ref;
     /*
        A 'typedef' can also be used to define compound or container
        objects, like lists, maps, and structures.  The standard KBase
        convention is to use structures, as shown here, to define the
        input and output of your function.  Here the input is a
        reference to the Assembly data object, a workspace to save
        output, and a length threshold for filtering.
    */

    typedef structure {
        assembly_ref assembly_input_ref;
        string workspace_name;
        boolean showContigs;
    } AssemblyMetadataReportParams;

    typedef structure {
        genome_ref genome_input_ref;
        string workspace_name;
        string report_format;
    } GenomeReportParams;

    typedef structure {
        genomeset_ref genomeset_input_ref;
        string workspace_name;
        string report_format;
    } GenomeSetReportParams;

    typedef structure {
        domain_ref domain_annotation_input_ref;
        float evalue_cutoff;
        string workspace_name;
        string report_format;
    } DomainReportParams;

    typedef structure {
        tree_ref tree_input_ref;
        string workspace_name;
        string report_format;
    } TreeReportParams;

    typedef structure {
        featseq_ref feature_sequence_input_ref;
        string workspace_name;
        string report_format;
    } FeatSeqReportParams;

    typedef structure {
        protcomp_ref protcomp_input_ref;
        string workspace_name;
        string report_format;
    } ProtCompReportParams;

    /*
        Here is the definition of the output of the function.  The output
        can be used by other SDK modules which call your code, or the output
        visualizations in the Narrative.  'report_name' and 'report_ref' are
        special output fields- if defined, the Narrative can automatically
        render your Report.
    */
    
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */

    
    funcdef assembly_metadata_report(AssemblyMetadataReportParams params)
        returns (ReportResults output) authentication required;
    funcdef genome_report(GenomeReportParams params)
        returns (ReportResults output) authentication required;
    funcdef genomeset_report(GenomeSetReportParams params)
        returns (ReportResults output) authentication required;
    funcdef domain_report(DomainReportParams params)
        returns (ReportResults output) authentication required;
    funcdef tree_report(TreeReportParams params)
        returns (ReportResults output) authentication required;
    funcdef featseq_report(FeatSeqReportParams params)
        returns (ReportResults output) authentication required;
    funcdef protcomp_report(ProtCompReportParams params)
        returns (ReportResults output) authentication required;
};
