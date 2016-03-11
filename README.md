# banzai
Scripts for analyzing Illumina-generated environmental DNA sequence data.
!!!! THESE FILES WORK ONLY FOR THE MBARI BOG DATA !!!!

Modifications to the banzai.sh script to process MBARI BOG MBON DNA sequence files.
BOG data are indexed and are already demultiplexed so the demultiplexing and cutadap calls in the original banzai.sh file have been
commented out.  Additional modifications are:
 Copy BOG fastq files from the network to the mbonserver, run PEAR, then delete the copied file.
 Add MEGAN taxonomic ranks, Kingdom, Class, Order, and Species.
 Manipulate all taxonomic ranks to create input for biom files, also create Phinch style ranking.
 Create PEAR report summary with mean assembled and non assembled percentages
 Create logfile2.txt without annoying read output
 Remove some of the really big files after processing such as nosingle.txt, 1_demult_concat.fasta.derep, 1_demult_concat.fasta.gz,
  and the dup_temp directory
 Added biom metadata file.
 Modifications to analyses_prelim.R:
   Create OTU_table.txt with taxonomic ranking for biom file.
   Add sample_ids to OTU table heading instead of tag ids
   Add debugging print statements
   Comment out graphing for community level indices as I couldn't get it to work on the mbonlive server
   
New files are:
  Merge_taxa_ranks.py - Make observational metadata, concatenate taxa ranks
  merge_otu_taxa.py - Merge taxonomy to OTU table.  Does not include no hits, not assigned.  Output is OTU_table_taxa.txt
  merge_otu_taxa2.py - Merge taxonomy to OTU table.  Does not include no hits, not assigned.  Output is OTU_table_taxa_phinch.txt 
    with Phinch prefixes denoting taxonomic rank.  i.e. p__Proteobacteria"
  merge_otu_taxa_all.py - Merge all OTU with taxonomy, even those without assignment
  mk_PEAR_report.py - Summary PEAR report
