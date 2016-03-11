# banzai
Scripts for analyzing Illumina-generated environmental DNA sequence data.

!!!! THESE FILES WORK ONLY FOR THE MBARI BOG DATA !!!!

Modifications to the banzai.sh script to process MBARI BOG MBON DNA sequence files.
BOG data are indexed and are already demultiplexed so the demultiplexing and cutadap calls in the original banzai.sh file have been
commented out.<br><br>  Additional modifications are:<br>
 1. Copy BOG fastq files from the network to the mbonserver, run PEAR, then delete the copied file.<br>
 2. Add MEGAN taxonomic ranks, Kingdom, Class, Order, and Species.<br>
 3. Manipulate all taxonomic ranks to create input for biom files, also create Phinch style ranking.<br>
 4. Create PEAR report summary with mean assembled and non assembled percentages<br>
 Create logfile2.txt without annoying read output<br>
 5. Remove some of the really big files after processing such as nosingle.txt, 1_demult_concat.fasta.derep, 1_demult_concat.fasta.gz, and the dup_temp directory<br>
 6. Added biom metadata file.<br><br>
 Modifications to analyses_prelim.R:<br>
 1.  Create OTU_table.txt with taxonomic ranking for biom file.
 2.  Add sample_ids to OTU table heading instead of tag ids
 3.  Add debugging print statements
 4.  Comment out graphing for community level indices as I couldn't get it to work on the mbonlive server
 <br><br>
New files are:<br>
1. Merge_taxa_ranks.py - Make observational metadata, concatenate taxa ranks<br>
2. merge_otu_taxa.py - Merge taxonomy to OTU table.  Does not include no hits, not assigned.  Output is OTU_table_taxa.txt<br>
3. merge_otu_taxa2.py - Merge taxonomy to OTU table.  Does not include no hits, not assigned.  Output is OTU_table_taxa_phinch.txt with Phinch prefixes denoting taxonomic rank.  i.e. p__Proteobacteria"<br>
4. merge_otu_taxa_all.py - Merge all OTU with taxonomy, even those without assignment<br>
5. mk_PEAR_report.py - Summary PEAR report
6. 

Modifications will continue with each target primer, 28S, COI, 16S as well as testing and analysis scripts.
