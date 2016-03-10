#!/bin/bash
#cd /usr/bin/
echo $(date +%H:%M) "blastn /MBON/mbonteam/MBARI/reiko/processed/Analysis_20160224_1433/all_lib/OTUs_swarm/OTUs.fasta"
	blastn \
		-query "/MBON/mbonteam/MBARI/reiko/processed/Analysis_20160224_1433/all_lib/OTUs_swarm/OTUs.fasta" \
		-db "/MBON/blastdb/Silva/SILVA123_SSURef.ncbi.db" \
		-num_threads "8" \
		-perc_identity "97" \
		-word_size "20" \
		-evalue "1e-20" \
		-max_target_seqs "400" \
		-culling_limit "20" \
		-outfmt "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle"	\
		-out "/MBON/mbonteam/MBARI/reiko/processed/Analysis_20160224_1433/all_lib/OTUs_Kelp18s_silva_blast.txt"
echo $(date +%H:%M) "Output is /MBON/mbonteam/MBARI/reiko/processed/Analysis_20160224_1433/all_lib/OTUs_Kelp18s_silva_blast.txt"		