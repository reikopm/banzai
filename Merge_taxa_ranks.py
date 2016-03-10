#!/usr/bin/python
# The banzai pipeline invokes MEGAN, metagenomic analysis program.  The output from MEGAN are files with the Kingdom, Phylum, Class, Order, Genus, Species.  These are matched by DUP number and a taxonomic ranking is created.  This ranking is attached to an OTU table.
# Called from the main banzai.sh file, or i.e. BOG28S_banzai.sh
# Input is the directory path to the MEGAN output, usually the ../all_lib subdirectory
# i.e. python "/home/mbonteam/MBARI/reiko/scripts/Make_taxa_ranks.py" "${ANALYSIS_DIR}"/all_lib
# rpm Jan 2016
import sys, getopt, csv, os, glob

def main(argv):

	ddir = sys.argv[1]
#	print(ddir)
	os.chdir(ddir)
	print(os.getcwd())
	mapper = dict()
	filename=glob.glob('*meganout_Species_mod.csv')
	i=0
	if filename[0]=='N_meganout_Species_mod.csv':
		i=1
	print('filename', filename[i])
	with open(filename[i], 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]
			t = row[2]
			mapper[i] = t        


	filename=glob.glob('*meganout_Genus_mod.csv')
	i=0
	if filename[0]=='N_meganout_Genus_mod.csv':
		i=1
	print('filename', filename[i])
	with open(filename[i], 'r') as f2:			
		reader = csv.reader(f2)
		with open('Genus_CSV.csv', 'wb') as fo:
			writer = csv.writer(fo)
			for row in reader:
				row.pop(1)
				if row[0] in mapper:
					row.append( mapper[ row[0] ] )
				writer.writerow(row)
				
	mapper = dict()
	filename=  glob.glob('*meganout_Family_mod.csv')	
	i=0
	if filename[0]=='N_meganout_Family_mod.csv':
		i=1
	print('filename', filename[i])
	with open(filename[i], 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]
			t = row[2]
			mapper[i] = t        

	with open('Genus_CSV.csv', 'r') as f2:
		with open('Family_CSV.csv', 'wb') as fo:
			reader = csv.reader(f2)
			writer = csv.writer(fo)
			for row in reader:
				if row[0] in mapper:
					row.insert(1,mapper[ row[0] ] )
				writer.writerow(row)
				
	mapper = dict()
	filename=  glob.glob('*meganout_Order_mod.csv')
	i=0
	if filename[0]=='N_meganout_Order_mod.csv':
		i=1
	print('filename', filename[i])	
	with open(filename[i], 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]
			t = row[2]
			mapper[i] = t        

	with open('Family_CSV.csv', 'r') as f2:
		with open('Order_CSV.csv', 'wb') as fo:
			reader = csv.reader(f2)
			writer = csv.writer(fo)
			for row in reader:
				if row[0] in mapper:
					row.insert(1,mapper[ row[0] ] )
				writer.writerow(row)
	mapper = dict()
	filename=  glob.glob('*meganout_Class_mod.csv')
	i=0
	if filename[0]=='N_meganout_Class_mod.csv':
		i=1
	print('filename', filename[i])	
	with open(filename[i], 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]
			t = row[2]
			mapper[i] = t        

	with open('Order_CSV.csv', 'r') as f2:
		with open('Class_CSV.csv', 'wb') as fo:
			reader = csv.reader(f2)
			writer = csv.writer(fo)
			for row in reader:
				if row[0] in mapper:
					row.insert(1,mapper[ row[0] ] )
				writer.writerow(row)
				
	mapper = dict()
	filename=  glob.glob('*meganout_Phylum_mod.csv')
	i=0
	if filename[0]=='N_meganout_Phylum_mod.csv':
		i=1
	print('filename', filename[i])
	with open(filename[i], 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]
			t = row[2]

			mapper[i] = t        

	with open('Class_CSV.csv', 'r') as f2:
		with open('Phylum_CSV.csv', 'wb') as fo:
			reader = csv.reader(f2)
			writer = csv.writer(fo)
			for row in reader:
				if row[0] in mapper:
					row.insert(1,mapper[ row[0] ] )
				writer.writerow(row)

	mapper = dict()
	filename=  glob.glob('*meganout_Kingdom_mod.csv')
	i=0
	if filename[0]=='N_meganout_Kingdom_mod.csv':
		i=1
	print('filename', filename[i])	
	with open(filename[i], 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]			
			t = row[2]
			mapper[i] = t        

	with open('Phylum_CSV.csv', 'r') as f2:
		with open('phinch_obs_md.csv', 'wb') as fo:
			reader = csv.reader(f2)
			writer = csv.writer(fo)
			for row in reader:
				if row[0] in mapper:
					row.insert(1,mapper[ row[0] ] )
				writer.writerow(row)
				
	# with open('Phylum.csv', 'r') as f2:
		# with open('phinch_obs_md.csv', 'wb') as fo:
			# reader = csv.reader(f2)
			# writer = csv.writer(fo)
			# for row in reader:
				# row.extend('"')
				# writer.writerow(row)
						
if __name__ == "__main__":
   main(sys.argv[1:])