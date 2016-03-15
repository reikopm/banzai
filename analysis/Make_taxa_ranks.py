#!/usr/bin/python
# Merge MEGAN taxonomic ranks for phinch style biom file
import sys, getopt, csv, os

def main(argv):

	ddir = sys.argv[1]
#	print(ddir)
	os.chdir(ddir)
	print(os.getcwd())
	mapper = dict()
	with open('meganout_Species_mod.csv', 'r') as f1:
		reader = csv.reader(f1)
		for row in reader:
			i = row[0]
			t = row[2]
			mapper[i] = t        

	with open('meganout_Genus_mod.csv', 'r') as f2:
		with open('Genus_CSV.csv', 'wb') as fo:
			reader = csv.reader(f2)
			writer = csv.writer(fo)
			for row in reader:
				row.pop(1)
				if row[0] in mapper:
					row.append( mapper[ row[0] ] )
				writer.writerow(row)
				
	mapper = dict()
	with open('meganout_Family_mod.csv', 'r') as f1:
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

	with open('meganout_Order_mod.csv', 'r') as f1:
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

	with open('meganout_Class_mod.csv', 'r') as f1:
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
	with open('meganout_Phylum_mod.csv', 'r') as f1:
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
	with open('meganout_Kingdom_mod.csv', 'r') as f1:
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
				
						
if __name__ == "__main__":
   main(sys.argv[1:])