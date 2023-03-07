# kb_ObjectInfo release notes
=========================================
1.1.0
-----
Add selection boxes so users can select  multiple formats in one run.

Remove the Meta and Vertical options from GenomeSet. They don’t add anything unique.

1.0.1
-----
Needed for metrics
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
