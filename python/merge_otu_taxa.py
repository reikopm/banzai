#!/usr/bin/python
# useage :
#  python merge_otu_taxa.py /home/mbonteam/MBARI/reiko/processed/Analysis_20160108_0100/all_lib
#
# Input argument is the directory path to the taxonomic levels,ie /home/mbonteam/MBARI/reiko/processed/Analysis_20160108_0100/all_lib
# Merges the taxonomic levels to the OTU table output.  
# Output looks like:
#OTU_ID	10013c01_12	10013c02_12	10013c03_12
#OTU_1010	4	15	7	Proteobacteria;Gammaproteobacteria;Methylococcales;Methylococcaceae;Methylobacter;

import sys, getopt, csv, os

def main(argv):

	ddir = sys.argv[1]
	fname = sys.argv[2]
#	print(ddir)
	os.chdir(ddir)
	#print(os.getcwd())
	mapper = dict()
#	with open('BOG18S_obs_md.txt', 'r') as f1:
	with open(fname, 'r') as f1:
		reader = csv.reader(f1, delimiter='\t')
		for row in reader:
			#print(row)
			# Column 0 contains the match; but we only want the left-most (before semi-colon)
			i = row[0]
			# Column 1 contains the target value for output
			t = row[2]
			mapper[i] = t        
	i=1
	with open('./OTUs_swarm/OTU_table.txt', 'r') as f2:
		with open('OTU_table_taxa.txt', 'wb') as fo:
			reader = csv.reader(f2, delimiter='\t')
			writer = csv.writer(fo, delimiter='\t')
			for row in reader:
				if row[0] in mapper:
					row.append( mapper[ row[0] ] )
					writer.writerow(row)
				else:
					if i:
						row.append('	taxonomy')
						writer.writerow([row])
						i=0
				
						
if __name__ == "__main__":
    main(sys.argv[1:])
   #main(**vars(args))