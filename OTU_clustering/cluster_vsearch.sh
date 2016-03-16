#!/usr/bin/env bash

# This is a shell script to cluster a fasta file of sequences into OTUs using vsearch.
# It requires one and only one argument: The absolute path to a fasta file.
# This file should contain only unique sequences - no duplicates.

# example usage: bash cluster_vsearch.sh "/Users/threeprime/Desktop/Analysis_20150727_0004/all_lib/no_duplicates.fasta"

# "${1}" is the first argument
infile="${1}"

# 0 to 1, inclusive
cluster_percent="0.99"

# define output files (these will be in the same directory as the infile)
dir_out="${infile%/*}"/OTUs_vsearch
mkdir "${dir_out}"
out_fasta="${dir_out}"/OTUs.fasta
logfile="${dir_out}"/OTUs.log
out_uc="${dir_out}"/OTUs.uc

echo $(date +%H:%M) "Clustering OTUs..."

vsearch \
	--cluster_fast "${infile}" \
	--id "${cluster_percent}" \
	--centroids "${out_fasta}" \
	--fasta_width 0 \
	--sizein \
	--sizeout \
	--uc "${out_uc}" \
	--log "${logfile}"

# vsearch --cluster_fast "${infile}" --id "${cluster_percent}" --centroids "${out_fasta}" --fasta_width 0 --sizein --sizeout --uc "${out_uc}" --log "${logfile}"

# NOT NECESSARY WITH VSEARCH
# remove the annoying line breaks
# echo $(date +%H:%M) "Removing line breaks within fasta sequences generated by usearch (awk)..."
# awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${DEREP_INPUT%/*}"/9_OTUs_linebreaks.fasta > "${DEREP_INPUT%/*}"/9_OTUs.fasta
# usearch option -notmatched disappeared with version 8
# awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${DEREP_INPUT%/*}"/9_notmatched_linebreaks.fasta > "${DEREP_INPUT%/*}"/9_notmatched.fasta

# rm "${DEREP_INPUT%/*}"/9_OTUs_linebreaks.fasta # "${DEREP_INPUT%/*}"/9_notmatched_linebreaks.fasta

exit

################################################################################
# RESOLVE OTUS AND DUPLICATES
################################################################################
# Note that this drops sequences determined to be chimeras by usearch
DUPS_TO_OTUS="${DEREP_INPUT%/*}"/dups_to_otus.csv
awk -F'[\t;]' 'BEGIN{ print "Query,Match" } { if ($4 == "otu") {print $1 "," $1} else if ($4 == "match") { print $1 "," $7 } else if ($4 == "chimera") { print $1 "," "chimera"} }' "${UPARSE_OUT}" > "${DUPS_TO_OTUS}"

# Assign the path for the OTU table
OTU_table="${DEREP_INPUT%/*}"/OTU_table.csv

# Convert duplicate table to OTU table using R script (arguments: duplicate table, dup to otu table, otu table path, concatenated directory (obsolete?))
Rscript "$SCRIPT_DIR/dup_to_OTU_table.R" "${duplicate_table}" "${DUPS_TO_OTUS}" "${OTU_table}" "${CONCAT_DIR}"






# TAG_COLUMN_NAME="tag_sequence"
# LIBRARY_COLUMN_NAME="library"
#
# SEQUENCING_METADATA="/Users/threeprime/Desktop/20150717/libraries/kelly_lab/SEQUENCING_POOL_20150618.csv"
#
# LIB_COL=$(awk -F',' -v LIB_COL_NAME=$LIBRARY_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == LIB_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
# LIBS=$(awk -F',' -v LIBCOL=$LIB_COL 'NR>1 {print $LIBCOL}' $SEQUENCING_METADATA | sort | uniq)
#
# echo $LIBS
#
# TAG_COL=$(awk -F',' -v TAG_COL_NAME=$TAG_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == TAG_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
# TAGS=$(awk -F',' -v TAGCOL=$TAG_COL 'NR>1 {print $TAGCOL}' $SEQUENCING_METADATA | sort | uniq)
#
# echo $TAGS
