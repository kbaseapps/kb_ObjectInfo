# kb_ObjectInfo release notes
=========================================
1.2.0
-----
1. Add selection boxes for multiple formats

2. Combine tab and comma delimited options to all objects where available. Create both outputs with minimal code change.

3. Assembly - In the HTML format, make the top section into a table and make the FastA section into preformatted text.

4. Add GenomeComparison Data format in Downloaded Excel or tsv file

5. Add formatting comments from previous review 

6. Remove the tree option - It is now available through regular export

7. Remove the old_git directory 

8. Remove the test file that ends in .orig 

9. Remove the Meta and Vertical options from GenomeSet. They don’t add anything unique.

10. Remove the old image files (kaleidoscope.png) left over from previous version.

1.1.0
-----
1. Whenever csv or tsv are generated, use the Python csv module instead of a write statement. This will avoid any problems with internal commas or tab that may occur in the data.

2. Make improvements in the usefulness and readability of the output. For example, instead of just giving an ID, attempt to lookup the name as well

3. Make csv an option in cases where tab is the only current option.

4. Make the html output a table instead of plain text. It will improve readability. Although it won’t help users tell the difference between tab and csv output options, it should be easier for those that view it in the browser.

5. Verify that all the output options are complete, 1) the PREVIEW in the report section, 2) the HTML viewing options, and 3) the download option. Some of these may have been incomplete before because they were not fully tested.

6. Move get_assembly_sequence from the Impl file to CreateFasta_Report.py because it is doing similar job.

7. The move from v8 genomes to v9 genomes made several changes. The change from the string ‘function’ to the list ‘functions’ has already been added. There was also a change that added ‘cdss’ and ‘non_coding_features’ as subsets of features. The cdss is a better place to get the protein translations and mRNAs because they don’t rely on checking the type.

1.0.1
-----
Changes needed for metrics
Change the params to have "input_ref" instead of "genome_input_ref" (or "[type]_input_ref"). 

Assembly Object Info - the "Include a FASTA of the Contigs", the FASTA should be a separate 
file with a .fna extension instead of making it all one file.

Comment out the print statement in the assembly section that prints all the output to the log.
Review all the other print statements in the other apps.

Remove the print statements. Use import logging and log statements instead to add comments to the log.

1.0.0
-----
Initial version of the module copied from the third-party version Report_util_landml

- Update the apps to recognize the newer version of the genomes that have multiple functions for features.
- Update the version and icon to match current standards.
- Update the names of the apps to a better description (e.g., Object Info instead of Text Report)
- Update to current version of SDK and python.
