#!/usr/bin/env bash
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 
# FOR MBARI BOG PROCESSING ONLY !!!!!
#
# This script is customized to run the BOG18S sequences and will not work with other sequence data.
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Pipeline for analysis of MULTIPLEXED Illumina data, a la Jimmy
# Modified for the MBARI Biological Oceanography Group 18S sequences.  Data are collected from CTD cruises in Monterey Bay, CA.  eDNA extraction and PCR processing are carried out at MBARI by Kris Walz.  Stanford University processes the Illumina sequencing and the raw sequences, as *.fastq files, are run through the Banzai Pipeline.  Our raw sequences are demulitplexed and so the demultiplexing and the cutadapt processing has been commented out in this script file.  Additionally our raw sequences do not have the barcode tag attached to the fastq sequence headers.  Barcodes are attached after the quality filtering.  Modifications were made by Reiko Michisaki, Dec 2015
# Runs all data files in specified PARENT_DIR. Modifications made to unzip *.fastq.gz files to /home/mbonteam/MBARI/reiko/raw/BOG_Data then delete the unzipped file after PEAR processing.  This will keep the processing from filling up the disk.  Reiko Michisaki Jan 2016

echo
echo
echo -e '\t' "\x20\xf0\x9f\x8f\x84" " "  "\xc2\xa1" BANZAI !
echo
echo


################################################################################
# CHECK FOR RAW DATA 
################################################################################

# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)

# MAKE MODIFICATIONS HERE 
# For SCRIPT_DIR, DATA_DIR and MY_SCRIPT_DIR.  SCRIPT_DIR is the location of the banzai repository.  For mbonteam this should point to /MBON/mbonteam/dev DATA_DIR setup is specific to MBARI BOG data which is stored on a linked network drive.  Data are copied to DATA_DIR for processing and then delete.  MY_SCRIPT_DIR contains custom banzai.sh, banzai_params.sh, and metadata for target loci.
CURR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT_DIR=$(dirname $CURR_DIR)

DATA_DIR="/home/reiko/MBARI/reiko/raw/BOG_data"
MY_SCRIPT_DIR="/home/reiko/MBARI/reiko/scripts"

# Read in the parameter file (was source "$SCRIPT_DIR/banzai_params.sh"; now argument 1)
param_file="${1}"
source "${param_file}"

# check if param file exists:
if [[ -s "${param_file}" ]] ; then
	echo "Reading analysis parameters from:"
	echo "${param_file}"
else
	echo
	echo 'ERROR! Could not find analysis parameter file. You specified the file path:'
	echo
	echo "${param_file}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi


# check if sequencing metadata exists
if [[ -s "${SEQUENCING_METADATA}" ]] ; then
	echo "Reading sequencing metadata from:"
	echo "${SEQUENCING_METADATA}"
else
	echo
	echo 'ERROR! Could not find sequencing metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi

################################################################################
# CHECK FOR DEPENDENCIES
################################################################################
dependencies=($( echo pear cutadapt vsearch swarm seqtk python blastn R ))
echo 'Checking for dependencies:' "${dependencies[@]}"
for i in "${dependencies[@]}"; do
	if hash "${i}" 2>/dev/null; then
	# if command -v "${i}" >/dev/null 2>&1; then
		echo 'Found program' "${i}" 'in' $( which "${i}" )
	else
		echo 'ERROR: A program on which this script depends was not found:' "${i}"
		echo 'Aborting script.'
		exit
	fi
done

# Specify compression utility
if hash pigz 2>/dev/null; then
	ZIPPER="pigz"
	echo "pigz installation found"
else
	ZIPPER="gzip"
	echo "pigz installation not found; using gzip"
fi

# Detect number of cores on machine; set variable
n_cores=$(getconf _NPROCESSORS_ONLN)
if [ $n_cores -gt 1 ]; then
	echo "$n_cores cores detected."
else
	n_cores=1
	echo "Multiple cores not detected."
fi

# make an analysis directory with starting time timestamp
ANALYSIS_DIR="${ANALYSIS_DIRECTORY}"/Analysis_"${START_TIME}"
mkdir "${ANALYSIS_DIR}"

# Write a log file of output from this script (everything that prints to terminal)
LOGFILE="${ANALYSIS_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1

echo "Analysis started at ""${START_TIME}" " and is located in ""${ANALYSIS_DIR}"

# Copy these files into that directory as a verifiable log you can refer back to.
#cp /home/mbonteam/MBARI/reiko/scripts/banzai_BOG18S.sh "${ANALYSIS_DIR}"/analysis_script.txt
cp $BASH_SOURCE "${ANALYSIS_DIR}"/analysis_script.txt
cp "${param_file}" "${ANALYSIS_DIR}"/analysis_parameters.txt



################################################################################
# LOAD MULTIPLEX TAGS
################################################################################
TAG_COL=$(awk -F',' -v TAG_COL_NAME=$TAG_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == TAG_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
TAGS=$(awk -F',' -v TAGCOL=$TAG_COL 'NR>1 {print $TAGCOL}' $SEQUENCING_METADATA | uniq)
ii=0
for T in $TAGS; do
	TAGS_SEQ[ii]=$T
	ii=$ii+1
done
#echo "${TAGS}"
R1_COL=$(awk -F',' -v R1_COL_NAME=$R1_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == R1_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
R1=$(awk -F',' -v R1COL=$R1_COL 'NR>1 {print $R1COL}' $SEQUENCING_METADATA | uniq)
R2_COL=$(awk -F',' -v R2_COL_NAME=$R2_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == R2_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
R2=$(awk -F',' -v R2COL=$R2_COL 'NR>1 {print $R2COL}' $SEQUENCING_METADATA | uniq)
echo "${R1}"
#TAGS=$(awk -F',' -v TAGCOL=$TAG_COL 'NR>1 {print $TAGCOL}' $SEQUENCING_METADATA | sort | uniq)
N_index_sequences=$(echo $TAGS | awk '{print NF}')

# check if number of tags is greater than one:
# if [[ "${N_index_sequences}" -gt 1 ]]; then
	 echo "Multiplex tags read from sequencing metadata (""${N_index_sequences}"") total"
# else
  # echo
  # echo 'ERROR:' "${N_index_sequences}" 'index sequences found. There should probably be more than 1.'
  # echo
  # echo 'Aborting script.'
	# exit
# fi

declare -a TAGS_ARRAY=($TAGS)


################################################################################
# Read in primers and create reverse complements. R2_COLUMN_NAME="R2"
################################################################################
PRIMER1_COLNUM=$(awk -F',' -v PRIMER1_COL=$PRIMER_1_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == PRIMER1_COL) print i; exit}' $SEQUENCING_METADATA)
PRIMER2_COLNUM=$(awk -F',' -v PRIMER2_COL=$PRIMER_2_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == PRIMER2_COL) print i; exit}' $SEQUENCING_METADATA)
PRIMER1=$(awk -F',' -v PRIMER1_COL=$PRIMER1_COLNUM 'NR==2 {print $PRIMER1_COL}' $SEQUENCING_METADATA)
PRIMER2=$(awk -F',' -v PRIMER2_COL=$PRIMER2_COLNUM 'NR==2 {print $PRIMER2_COL}' $SEQUENCING_METADATA)
echo "Primers read from sequencing metadata:" "${PRIMER1}" "${PRIMER2}"

if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
  echo 'Primers read from metadata columns' "${PRIMER1_COLNUM}" 'and' "${PRIMER2_COLNUM}"
  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
else
  echo 'ERROR:' 'At least one primer is not valid'
  echo 'Looked in metadata columns' "${PRIMER1_COLNUM}" 'and' "${PRIMER2_COLNUM}"
  echo 'Aborting script'
  exit
fi

# make primer array
read -a primers_arr <<< $( echo $PRIMER1 $PRIMER2 )

# Reverse complement primers
PRIMER1RC=$( echo ${PRIMER1} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )
PRIMER2RC=$( echo ${PRIMER2} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )

# make primer array
read -a primersRC_arr <<< $( echo $PRIMER1RC $PRIMER2RC )


################################################################################
# Calculate the expected size of the region of interest, given the total size of fragments, and the length of primers and tags
################################################################################
EXTRA_SEQ=${TAGS_ARRAY[0]}${TAGS_ARRAY[0]}$PRIMER1$PRIMER2
LENGTH_ROI=$(( $LENGTH_FRAG - ${#EXTRA_SEQ} ))
LENGTH_ROI_HALF=$(( $LENGTH_ROI / 2 ))
echo "**** Extra_seq=", "${EXTRA_SEQ}", " Length_ROI=", "${LENGTH_ROI}"

################################################################################
# Find raw sequence files
################################################################################
# Look for any file with '.fastq' in the name in the parent directory
# note that this will include ANY file with fastq -- including QC reports!
# TODO make LIBRARY_DIRECTORIES an array by wrapping it in ()
LIBRARY_DIRECTORIES=$( find "$PARENT_DIR" -name '*.fastq*' -print0 | xargs -0 -n1 dirname | sort --unique )

# PEAR v0.9.6 does not correctly merge .gz files.
# Look through files and decompress if necessary.
# raw_files=($( find "${PARENT_DIR}" -name '*.fastq*' ))
# for myfile in "${raw_files[@]}"; do
	# if [[ "${myfile}" =~ \.gz$ ]]; then
		# echo $(date +%H:%M) "decompressing "${myfile}""
		# "${ZIPPER}" -d "${myfile}"
	# fi
# done

# Count library directories and print the number found
# TODO if LIBRARY_DIRECTORIES is an array, its length is "${#LIBRARY_DIRECTORIES[@]}"
N_library_dir=$(echo $LIBRARY_DIRECTORIES | awk '{print NF}')
echo "${N_library_dir}"" library directories found:"

# Show the libraries that were found:
# TODO for i in "${LIBRARY_DIRECTORIES[@]}"; do echo "${i##*/}" ; done
#for i in $LIBRARY_DIRECTORIES; do echo "${i##*/}" ; done

# Assign it to a variable for comparison
LIBS_FROM_DIRECTORIES=$(for i in $LIBRARY_DIRECTORIES; do echo "${i##*/}" ; done)

# Read library names from file or sequencing metadata
if [ "${READ_LIB_FROM_SEQUENCING_METADATA}" = "YES" ]; then
	LIB_COL=$(awk -F',' -v LIB_COL_NAME=$LIBRARY_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == LIB_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
	LIBS=$(awk -F',' -v LIBCOL=$LIB_COL 'NR>1 {print $LIBCOL}' $SEQUENCING_METADATA)
	N_libs=$(echo $LIBS | awk '{print NF}')
	echo "Library names read from sequencing metadata (""${N_libs}"") total"
	echo "${LIBS}"
else
	LIBS=$(tr '\n' ' ' < "${LIB_FILE}" )
	N_libs=$(echo $LIBS | awk '{print NF}')
	echo "Library names read from lib file (""${LIBS}"") total"
fi
# make library names into an array
# TODO LIBS_ARRAY is never used
# declare -a LIBS_ARRAY=($LIBS)

# Check that library names are the same in the metadata and file system
if [ "$LIBS_FROM_DIRECTORIES" != "$LIBS" ]; then
	echo "Warning: Library directories and library names in metadata are NOT the same. Something will probably go wrong later..."
else
	echo "Library directories and library names in metadata are the same - great jorb."
fi


# Unique samples are given by combining the library and tags
# TODO originally contained sort | uniq; this is unnecessary I think
#LIB_TAG_MOD=$( awk -F',' -v LIBCOL=$LIB_COL -v TAGCOL=$TAG_COL 'NR>1 { print "lib_" $LIBCOL "_tag_" $TAGCOL }' $SEQUENCING_METADATA | sort | uniq )
LIB_TAG_MOD_COL=$(awk -F',' -v LIB_TAG_MOD_COL_NAME=$LIBRARY_TAG_COMBO_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == LIB_TAG_MOD_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
LIB_TAG_MOD=$(awk -F',' -v LIBTAGMODCOL=$LIB_TAG_MOD_COL 'NR>1 {print $LIBTAGMODCOL}' $SEQUENCING_METADATA)
echo "LIB_TAG_MOD"
LIB_TAG_MOD
# create a file to store tag efficiency data
TAG_COUNT="${ANALYSIS_DIR}"/tag_count.txt
echo "library   tag   left_tagged   right_tagged" >> "${TAG_COUNT}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
# BOG data are on a network disk.  Copy file to mbonserver, process, then remove
################################################################################
ii=0
for C_LIB in $LIBS; do
	echo "Library " "${PARENT_DIR}"/"${C_LIB}"
	raw_files=($( find "${PARENT_DIR}"/"${C_LIB}" -name '*.fastq*' ))
 	for myfile in "${raw_files[@]}"; do
	 echo $(date +%H:%M) "copying "${myfile}""
	 cp "${myfile}"  "${DATA_DIR}"/"${myfile##*/}"
	 if [[ "${myfile}" =~ \.gz$ ]]; then
		 echo $(date +%H:%M) "decompressing "${myfile}""
		 "${ZIPPER}" -d "${DATA_DIR}"/"${myfile##*/}"
	 fi
	done	
   
   CURRENT_LIB="${DATA_DIR}"
   
 # Identify the forward and reverse fastq files.
	READS=($(find "${CURRENT_LIB}" -name '*.fastq*' | sort))
	READ1="${READS[0]}"
	READ2="${READS[1]}"

#	LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/} #old code
	LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${C_LIB}
	mkdir "${LIB_OUTPUT_DIR}"

	##############################################################################
	# MERGE PAIRED-END READS AND QUALITY FILTER (PEAR)
	##############################################################################

	##############################################################################
	# CALCULATE EXPECTED AND MINIMUM OVERLAP OF PAIRED END SEQUENCES
	##############################################################################
	LENGTH_READ=$( head -n 100000 "${READ1}" | awk '{print length($0);}' | sort -nr | uniq | head -n 1 )
	#echo " **** DEBUG " "${LENGTH_READ}"
	OVERLAP_EXPECTED=$(($LENGTH_FRAG - (2 * ($LENGTH_FRAG - $LENGTH_READ) ) ))
	MINOVERLAP=$(( $OVERLAP_EXPECTED / 2 ))

	##############################################################################
	# CALCULATE MAXIMUM AND MINIMUM LENGTH OF MERGED READS
	##############################################################################
	ASSMAX=$(( $LENGTH_FRAG + 50 ))
	ASSMIN=$(( $LENGTH_FRAG - 50 ))
	####  DEBUG  ###
	# COMMENT OUT SECTION FOR PRODUCTION.  These are the values Kevan used when manually running PEAR
	ASSMAX=300  #180
	ASSMIN=75   #120
	MINOVERLAP=50
	PHRED=33
	ET=0


	if [ "$ALREADY_PEARED" = "YES" ]; then
		MERGED_READS="$PEAR_OUTPUT"
		echo "Paired reads have already been merged."
	else
		echo $(date +%H:%M) "Merging reads in library" "${C_LIB}""..."
		MERGED_READS_PREFIX="${LIB_OUTPUT_DIR}"/1_merged
		MERGED_READS="${LIB_OUTPUT_DIR}"/1_merged.assembled.fastq
		pear \
			--forward-fastq "${READ1}" \
			--reverse-fastq "${READ2}" \
			--output "${MERGED_READS_PREFIX}" \
			-v $MINOVERLAP \
			-m $ASSMAX \
			-n $ASSMIN \
			-t $min_seq_length \
			-q $Quality_Threshold \
			-u $UNCALLEDMAX \
			-g $TEST \
			-p $PVALUE \
			-s $SCORING \
			-j $n_cores

		# check pear output:
		if [[ ! -s "${MERGED_READS}" ]] ; then
		    echo 'ERROR: No reads were merged.'
		    echo 'Maybe the script should exit, but this could be a library-specific problem.'
		fi



	fi

	################################################################################
	# EXPECTED ERROR FILTERING (usearch)
	################################################################################
	# FILTER READS (This is the last step that uses quality scores, so convert to fasta)
	if [ "${Perform_Expected_Error_Filter}" = "YES" ]; then
		echo $(date +%H:%M) "Filtering merged reads..."
		FILTERED_OUTPUT="${LIB_OUTPUT_DIR}"/2_filtered.fasta
		vsearch \
			--fastq_filter "${MERGED_READS}" \
			--fastq_maxee "${Max_Expected_Errors}" \
			--fastaout "${FILTERED_OUTPUT}" \
			--fasta_width 0
	else
		# Convert merged reads fastq to fasta
		echo  $(date +%H:%M) "converting fastq to fasta..."
		FILTERED_OUTPUT="${MERGED_READS%.*}".fasta
		seqtk seq -A "${MERGED_READS}" > "${FILTERED_OUTPUT}"
	fi

	# Compress merged reads
  echo $(date +%H:%M) "Compressing PEAR output..."
  #find "${LIB_OUTPUT_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
  rm "${MERGED_READS}"
  echo $(date +%H:%M) "PEAR output compressed."





	if [ "${RENAME_READS}" = "YES" ]; then
		echo $(date +%H:%M) "Renaming reads in library" "${C_LIB}""..."
		# TODO remove whitespace from sequence labels?
		# sed 's/ /_/'

		# original (usearch7):	sed -E "s/ (1|2):N:0:[0-9]/_"${CURRENT_LIB##*/}"_/" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/tmp.fasta
		# update for usearch8, which without warning removes any part of the sequence ID following a space.
		# holy shit this ads the _"${CURRENT_LIB##*/}"_ to EVERY line
		# sed -E "s/$/_"${CURRENT_LIB##*/}"_/" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/tmp.fasta
		# sed -E "s/>([a-zA-Z0-9-]*:){4}/>/" "${CURRENT_LIB}"/tmp.fasta > "${FILTERED_OUTPUT%.*}"_renamed.fasta
		# rm "${CURRENT_LIB}"/tmp.fasta

		# updated 20150521; one step solution using awk; removes anything after the first space!
		FILTERED_RENAMED="${FILTERED_OUTPUT%.*}"_renamed.fasta
		awk -F'[: ]' '{
				if ( /^>/ )
					print ">"$4":"$5":"$6":"$7"'${C_LIB}'_tag_";
				else
					print $0
				}' "${FILTERED_OUTPUT}" > "${FILTERED_RENAMED}"

		FILTERED_OUTPUT="${FILTERED_RENAMED}"

		# rm "${FILTERED_OUTPUT}"

	else
		echo "Reads not renamed"
	fi


	################################################################################
	# HOMOPOLYMERS (grep, awk)
	################################################################################
	if [ "${REMOVE_HOMOPOLYMERS}" = "YES" ]; then
		echo $(date +%H:%M) "Removing homopolymers..."
		HomoLineNo="${CURRENT_LIB}"/homopolymer_line_numbers.txt
		grep -E -i -B 1 -n "(A|T|C|G)\1{$HOMOPOLYMER_MAX,}" "${FILTERED_OUTPUT}" | \
			cut -f1 -d: | \
			cut -f1 -d- | \
			sed '/^$/d' > "${HomoLineNo}"
		if [ -s "${HomoLineNo}" ]; then
			DEMULTIPLEX_INPUT="${CURRENT_LIB}"/3_no_homopolymers.fasta
			awk 'NR==FNR{l[$0];next;} !(FNR in l)' "${HomoLineNo}" "${FILTERED_OUTPUT}" > "${DEMULTIPLEX_INPUT}"
			awk 'NR==FNR{l[$0];next;} (FNR in l)' "${HomoLineNo}" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/homopolymeric_reads.fasta
		else
			echo "No homopolymers found" > "${CURRENT_LIB}"/3_no_homopolymers.fasta
			DEMULTIPLEX_INPUT="${FILTERED_OUTPUT}"
		fi
	else
		echo "Homopolymers not removed."
		DEMULTIPLEX_INPUT="${FILTERED_OUTPUT}"
	fi

	################################################################################
	# DEMULTIPLEXING (awk)
	# BOG 18S modification
	# Raw sequences are already demultiplexed
	################################################################################
	# make a directory to put all the demultiplexed files in
	DEMULTIPLEXED_DIR="${LIB_OUTPUT_DIR}"/demultiplexed
	mkdir "${DEMULTIPLEXED_DIR}"

	# Copy sequences to fasta files into separate directories based on tag sequence on left side of read
	# TODO test for speed against removing the tag while finding it: wrap first tag regex in gsub(/pattern/,""):  awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,""){if . . .
	# 20150522 changed {0,9} to {3} to eliminate flexibility (that could result in a read being assigned to >1 sample)
	# awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &

	echo $(date +%H:%M) "Demultiplexing: removing tags and adding to sequence ID in library" "${C_LIB}""..."
	#for TAG_SEQ in $TAGS; do
	TAG_SEQ="${TAGS_SEQ[$ii]}"
	(	TAG_DIR="${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
		mkdir "${TAG_DIR}"
		demult_file_L="${TAG_DIR}"/1_tagL_removed.fasta
	    demult_file_R="${TAG_DIR}"/2_notags.fasta

	# Raw sequences are already demultiplexed so skip this step
	# NOTE BOG samples/individual barcoded raw seq files MUST be run individually !!!!  re: TAGS
		# Left side tag
		# awk 'gsub(/^.{3}'"$TAG_SEQ"'/,"") {
			# if (a && a !~ /^.{3}'"$TAG_SEQ"'/)
				# print a;
			# print
		# } {a=$0}' "${DEMULTIPLEX_INPUT}" > "${demult_file_L}"

		# # Right side tag
		# TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
		# awk 'gsub(/'"$TAG_RC"'.{3}$/,"") {
			# if (a && a !~ /'"$TAG_RC"'.{3}$/)
				# print a "tag_""'"$TAG_SEQ"'";
			# print
		# } {a = $0}' "${TAG_DIR}"/1_tagL_removed.fasta > "${demult_file_R}"

		# echo "${CURRENT_LIB##*/}" "${TAG_SEQ}" $(wc -l "${demult_file_L}" | \
			# awk '{ print ($1/2) }') $(wc -l "${demult_file_R}" | \
			# awk '{ print ($1/2)}') >> "${TAG_COUNT}" ) &
			echo " filter output sed tag " "$TAG_SEQ"
		cat "${FILTERED_OUTPUT}" | sed 's/_tag_/_tag_'"${TAG_SEQ}"'/'> "${demult_file_R}"
		
    )
	#done



	wait

	# echo $(date +%H:%M) "Demultiplexing: removing right tag and adding tag sequence to sequence ID in library" "${CURRENT_LIB##*/}""..."
	# for TAG_SEQ in $TAGS; do
	# (	TAG_DIR="${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	# 	demult_file_R="${TAG_DIR}"/2_notags.fasta
	# 	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
	# 	# 20150522 changed {0,9} to {3} to eliminate flexibility (that could result in a read being assigned to >1 sample)
	# 	awk 'gsub(/'"$TAG_RC"'.{3}$/,"") {if (a && a !~ /'"$TAG_RC"'.{3}$/) print a "tag_""'"$TAG_SEQ"'"; print } {a = $0}' "${TAG_DIR}"/1_tagL_removed.fasta > "${demult_file_R}" ) &
	# done
	#
	# wait
	#remove the unzipped fastq sequence data files
    rm "${DATA_DIR}"/*.fastq
	ii=$ii+1
done
################################################################################
# END LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################










# TODO add single if/else for CONCATENATE_SAMPLES: assign directory as appropriate, correct references within loop to be extendable
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	# do the stuff here
	WORKING_DIR="${CONCAT_DIR}"
else
	WORKING_DIR="${LIBRARY_DIRECTORIES}"
fi

################################################################################
# CONCATENATE SAMPLES
################################################################################
# TODO could move this first step up above any loops (no else)
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then

	# TODO MOVE THE VARIABLE ASSIGNMENT TO TOP; MOVE MKDIR TO TOP OF CONCAT IF LOOP
	echo $(date +%H:%M) "Concatenating fasta files..."
	CONCAT_DIR="$ANALYSIS_DIR"/all_lib
	mkdir "${CONCAT_DIR}"
	CONCAT_FILE="${CONCAT_DIR}"/1_demult_concat.fasta

	# TODO could move this into above loop after demultiplexing?
	ii=0
	for CURRENT_LIB in $LIBS; do

		LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB}

		TAG_SEQ="${TAGS_SEQ[$ii]}"
		cat "${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"/2_notags.fasta >> "${CONCAT_FILE}"
		echo "${TAG_SEQ}"
		

		echo $(date +%H:%M) "Compressing fasta files..." "${LIB_OUTPUT_DIR}"
		find "${LIB_OUTPUT_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
		echo $(date +%H:%M) "fasta files compressed."
		ii=$ii+1
	done

	################################################################################
	# Count the occurrences of '_tag_' + the 12 character bar code tag following it in the concatenated file
	################################################################################
  # TODO !!! This will fail if there are underscores in the library names !!!
	# an attempt at making this robust to underscores
	# grep -E -o '_lib_.+?(?=_tag)_tag_.{6}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags

	echo $(date +%H:%M) "Counting reads associated with each sample index (primer tag)..."
	grep -E -o '_tag_.{12}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags
	#grep -E -o '_lib_[^_]*_tag_.{6}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags
	sed -i 's/_tag_/_/' "${CONCAT_DIR}"/1_demult_concat.fasta
	echo $(date +%H:%M) "Summary of sequences belonging to each sample index found in ""${CONCAT_DIR}""/1_demult_concat.fasta.tags"
	################################################################################




	################################################################################
	# PRIMER REMOVAL
	################################################################################
	# (moot for concatenated file): echo $(date +%H:%M) "Removing primers in library" "${CURRENT_LIB##*/}""..."
	# Remove PRIMER1 and PRIMER2 from the BEGINNING of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.

	# count lines in primer removal input
	echo $(date +%H:%M) "Counting sequences in primer removal input..."
	seq_N_demult_concat=$( grep -e '^>' --count "${CONCAT_FILE}" )
	echo $(date +%H:%M) "${seq_N_demult_concat}" "sequences found in file" "${CONCAT_FILE}"

	# TODO wrap in '( ) &' to force into background and allow parallel processing
	# i.e.
	# for primer in "${primers_arr[@]}"; do
	# 	( cutadapt -g ^"${primer}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_FILE}" > "${CONCAT_DIR}"/5_L"${primer}"_removed.fasta ) &
	# done
	# wait

	# # remove primer 1 from left side of sequences
	# primerL1_removed="${CONCAT_DIR}"/5_primerL1_removed.fasta
	# ( cutadapt \
		# -g ^"${PRIMER1}" \
		# -e "${PRIMER_MISMATCH_PROPORTION}" \
		# # -m "${LENGTH_ROI_HALF}" \
# #		--discard-untrimmed \
		# "${CONCAT_FILE}" > "${primerL1_removed}" ) &

	# # remove primer 2 from left side of sequences
	# primerL2_removed="${CONCAT_DIR}"/5_primerL2_removed.fasta
	# ( cutadapt \
		# -g ^"${PRIMER2}" \
		# -e "${PRIMER_MISMATCH_PROPORTION}" \
		# # -m "${LENGTH_ROI_HALF}" \
# #		--discard-untrimmed \
		# "${CONCAT_FILE}" > "${primerL2_removed}" ) &

	# wait

	# # compress left primer removal input  DEBUG
	# echo $(date +%H:%M) "Compressing left primer removal input..."
	# # "${ZIPPER}" "${CONCAT_DIR}"/1_demult_concat.fasta
	# echo $(date +%H:%M) "Left primer removal input compressed."

	# # check for cutadapt/primer removal success.
	# if [[ ! -s "${primerL1_removed}" ]]; then
	  # echo 'ERROR: cutadapt did not process primerL1 reads correctly. This file is empty or absent:'
		# echo "${primerL1_removed}"
	  # # echo 'Aborting script'
	  # # exit
	# fi
	# # check for cutadapt/primer removal success.
	# if [[ ! -s "${primerL2_removed}" ]]; then
	  # echo 'ERROR: cutadapt did not process primerL2 reads correctly. This file is empty or absent:'
		# echo "${primerL2_removed}"
	  # # echo 'Aborting script'
	  # # exit
	# fi

	# # Remove the reverse complement of primer 1 from the right side of sequences
	# primerR1_removed="${CONCAT_DIR}"/6_primerR1_removed.fasta
	# ( cutadapt \
		# -a "${PRIMER2RC}"$ \
		# -e "${PRIMER_MISMATCH_PROPORTION}" \
		# # -m "${LENGTH_ROI_HALF}" \
# #		--discard-untrimmed \
		# "${primerL1_removed}" > "${primerR1_removed}" ) &

	# # Remove the reverse complement of primer 2 from the right side of sequences
	# primerR2_removed="${CONCAT_DIR}"/6_primerR2_removed.fasta
	# ( cutadapt \
		# -a "${PRIMER1RC}"$ \
		# -e "${PRIMER_MISMATCH_PROPORTION}" \
		# # -m "${LENGTH_ROI_HALF}" \
# #		--discard-untrimmed \
		# "${primerL2_removed}" > "${primerR2_removed}" ) &

	# wait

	# # check for cutadapt/primer removal success.
	# if [[ ! -s "${primerR1_removed}" ]]; then
		# echo 'ERROR: cutadapt did not process primerR1 reads correctly. This file is empty or absent:'
		# echo "${primerR1_removed}"
		# # echo 'Aborting script'
		# # exit
	# fi
	# # check for cutadapt/primer removal success.
	# if [[ ! -s "${primerR2_removed}" ]]; then
		# echo 'ERROR: cutadapt did not process primerR2 reads correctly. This file is empty or absent:'
		# echo "${primerR2_removed}"
		# # echo 'Aborting script'
		# # exit
	# fi

	# Remove reverse-complement the sequences in which the RC of primer 1 was found on the right side
	# seqtk seq -r "${CONCAT_DIR}"/6_primerR1_removed.fasta > "${CONCAT_DIR}"/6_primerR1_removedRC.fasta

	# paste together the contents of the files that primers were removed from.  DEBUG
	# DEREP_INPUT="${CONCAT_DIR}"/7_no_primers.fasta

	# cat "${CONCAT_DIR}"/6_primerR1_removedRC.fasta "${CONCAT_DIR}"/6_primerR2_removed.fasta > "${DEREP_INPUT}"
    DEREP_INPUT="${CONCAT_DIR}"/1_demult_concat.fasta
	# check that it worked (derep input / no primers)
	if [[ ! -s "${DEREP_INPUT}" ]] ; then
	    echo 'ERROR: Input file for dereplication is empty or absent.'
	    echo 'This will cause problems for all remaining steps, so script will exit.'
	    exit
	fi



	################################################################################
	# CONSOLIDATE IDENTICAL SEQUENCES (DEREPLICATION)
	################################################################################
	echo $(date +%H:%M) "Identifying identical sequences... python /home/reiko/banzai/dereplication/dereplicate_fasta.py" "${DEREP_INPUT}"
	derep_output="${DEREP_INPUT}".derep
	python "$SCRIPT_DIR"/dereplication/dereplicate_fasta.py "${DEREP_INPUT}"

	# check for derep output
	if [[ ! -s "${derep_output}" ]] ; then
	    echo 'ERROR: python dereplication output is empty or absent.'
	    echo 'This will cause problems for all remaining steps, so script will exit.'
	    exit
	fi


	# Exclude singleton sequences (if NF > 2), count the number of sequences per duplicate (print NF-1), sort them by the number of sequences per duplicate (sort -nr), and precede with a name ("DUP_X", where X is the line number,sorted)
	echo $(date +%H:%M) "Counting duplicates per identical sequence and excluding singletons... (awk)"
	no_singletons="${DEREP_INPUT%/*}"/nosingle.txt

	awk -F';' '{
		if (NF > 2)
			print NF-1 ";" $0
		}' "${derep_output}" | \
	sort -nr | \
	awk -F';' '{
		print ">DUP_" NR ";" $0
	}' > "${no_singletons}"

	# check output
	if [[ ! -s "${no_singletons}" ]] ; then
	    echo 'There was a problem generating the nosingletons file. It is empty or absent.'
	    echo 'This will cause problems counting sequences for dereplication.'
	fi



	# COUNT OCCURRENCES PER SAMPLE (LIBRARY + TAG) PER DUPLICATE

	# assign a path for the output (a table of counts of each duplicate sequence in each unique combination of library and primer ("tag") indexes, and a fasta of all the duplicate sequences.
	duplicate_table="${DEREP_INPUT%/*}"/duplicate_table.csv
	duplicate_fasta="${DEREP_INPUT%/*}"/duplicates.fasta

	# make a directory to store the temporary duplicate files
	temp_dir="${DEREP_INPUT%/*}"/dup_temp
	mkdir "${temp_dir}"

	# set a file prefix for the batches of samples
	sample_batch_prefix="${temp_dir}"/sample_batch_

	# split the sample identifiers (lib + tag combination) into batches of no more than the number of available cores
	echo $LIB_TAG_MOD | tr ' ' '\n' | split -l "${n_cores}" - "${sample_batch_prefix}"


	# for each of the batches of files
	for batch in "${sample_batch_prefix}"* ; do

		echo processing batch "${batch##*/}"

		# 	current_batch=$( cat "${batch}" ) # this reads whitespace rather than newline

		for sample in $( cat "$batch" ) ; do

			# say that it's being processed
			echo $(date +%H:%M) "Processing" "${sample}""..."

			# Isolate this process to be put in the background
			(

			# write an output file called *.dup, start by printing the lib/tag being processed, then print a count the occurrences of the current lib/tag on each line of the input file
			awk 'BEGIN {print "'$sample'" ; FS ="'${sample}'" } { print NF -1 }' "${no_singletons}" > "${temp_dir}"/"${sample}".dup

			) &

		done

		wait

	done

	# write a file of names of each of the duplicates:
	dupnames="${temp_dir}"/dupnames
	awk -F';' 'BEGIN {print "sample"} {print $1}' "$no_singletons" | sed 's/>//' > "${dupnames}"

	# first, count the number of duplicate files:
	n_files=$(find "${temp_dir}" -type f -name '*.dup*' | wc -l)

	# I think this was only relevant for a different approach
	max_files=$(ulimit -n)

	# this will paste row by row... takes 48s on a set of 300 files (samples) each containing 630023 lines (duplicates)
	paste -s -d, "${dupnames}" "${temp_dir}"/*.dup > "${duplicate_table}"

	# this will do columns; it takes a very long time.
	# for file in "${temp_dir}"/*; do cat final.dup | paste - $file >temp; cp temp final.dup; done; rm temp

	# cleanup
	# rm "${infile%/*}"/*.dup
	# rm "${sample_batch_prefix}"*

	# say that you're finished.
	echo $(date +%H:%M) "Identical sequences consolidated in file ""${duplicate_table}"

	# OLD/UNSTABLE/BAD -- but generated CSV in different orientation
	# This will start as many processes as you have libraries... be careful!
	# echo $(date +%H:%M) "Consolidating identical sequences per sample... (awk)"
	# for CURRENT_LIB in $LIBRARY_DIRECTORIES; do
	# 	( for TAG_SEQ in $TAGS; do
	# 		LIB_TAG="lib_${CURRENT_LIB##*/}_tag_${TAG_SEQ}"
	# 		echo $(date +%H:%M) "Processing" "${LIB_TAG}""..."
	# 		# the output of the awk function gsub is the number of replacements, so you could use this instead... however, it appears slower?
	# 		# ( awk 'BEGIN {print "'$LIB_TAG'" } { print gsub(/"'$LIB_TAG'"/,"") }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/"${LIB_TAG}".dup ) &
	# 		awk 'BEGIN {print "'$LIB_TAG'" ; FS ="'${LIB_TAG}'" } { print NF -1 }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/"${LIB_TAG}".dup
	# 	done ) &
	# done
	# wait
	# # Write a csv file of the number of occurrences of each duplicate sequence per tag. (rows = sequences, cols = samples)
	# duplicate_table="${DEREP_INPUT%/*}"/dups.csv
	# find "${DEREP_INPUT%/*}" -type f -name '*.dup' -exec paste -d, {} \+ | awk '{ print "DUP_" NR-1 "," $0 }' > "${duplicate_table}"
	# # delete all of the '.dup' files
	# rm ${DEREP_INPUT%/*}/*.dup

	# Write fasta file in order to blast sequences
	echo $(date +%H:%M) "Writing fasta file of duplicate sequences"
	awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' "${no_singletons}" > "${duplicate_fasta}"

	# check if duplicate fasta and duplicate table exist. (Might need to check size)
	if [[ ! -s "${duplicate_fasta}" ]] ; then
	    echo 'There was a problem generating the duplicate fasta file. It is empty or absent.'
	    echo 'The remainder of the script, including OTU clustering, depends on this file.'
	    echo 'Aborting script.'
	    exit
	fi
	if [[ ! -s "${duplicate_table}" ]] ; then
	    echo 'There was a problem generating the duplicate table. It is empty or absent.'
	    echo 'Aborting script.'
	    exit
	fi

################################################################################


	################################################################################
	# CLUSTER OTUS
	################################################################################
	# Note that identical (duplicate) sequences were consolidated earlier;
	# This step outputs a file (*.uc) that lists, for every sequence, which sequence it clusters with
	if [ "$CLUSTER_OTUS" = "NO" ]; then
		BLAST_INPUT="${duplicate_fasta}"
	else
		echo $(date +%H:%M) "Clustering OTUs..."

		case "${cluster_method}" in

		    "swarm" )

		        echo $(date +%H:%M) 'Clustering sequences into OTUs using swarm'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh "${duplicate_fasta}"

		    ;;

		    "vsearch" )

		        # echo $(date +%H:%M) 'Clustering sequences into OTUs using vsearch'
		        # source "${SCRIPT_DIR}"/OTU_clustering/cluster_vsearch.sh "${duplicate_fasta}"
						echo "Sorry, OTU clustering with vsearch has not been implemented yet."
						echo $(date +%H:%M) 'Clustering sequences into OTUs using swarm'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh "${duplicate_fasta}"

		    ;;

		    "usearch" )

		        echo $(date +%H:%M) 'Clustering sequences into OTUs using usearch'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_usearch.sh "${duplicate_fasta}"

		    ;;

		    * )

		        echo "${cluster_method}" 'is an invalid clustering method.'
		        echo 'Must be one of swarm, vsearch, usearch, or none.'
		        echo $(date +%H:%M) 'Clustering sequences into OTUs using swarm'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh "${duplicate_fasta}"

		    ;;

		esac

		# check that dup to otu map is greater than 12 bytes
		minsize=12
		size_dup_otu_map=$(wc -c <"${dup_otu_map}")
		if [ $size_dup_otu_map -lt $minsize ]; then
		    echo 'There was an error generating the dup-to-otu map.'
		fi


		# Assign the path for the OTU table
		# OTU_table="${dir_out}"/OTU_table.csv

		# Convert duplicate table to OTU table using R script (arguments: (1) duplicate table, (2) dup to otu table, (3) otu table path
		echo 'dup_to_OTU_table.R processing ...'
		Rscript "$SCRIPT_DIR/dup_to_OTU_table.R" "${duplicate_table}" "${dup_otu_map}" "${OTU_table}"

		# check if OTU table and OTU fasta exist (and/or are of size gt 1?)
		if [[ ! -s "${OTU_fasta}" ]] ; then
		    echo 'There was a problem generating the OTU fasta file. It is empty or absent.'
		    echo 'Aborting script.'
		    exit
		fi
		if [[ ! -s "${OTU_table}" ]] ; then
		    echo 'There was a problem generating the OTU table. It is empty or absent.'
		    echo 'Aborting script.'
		    exit
		fi

	fi

	##############################################################################
	# CHECK FOR CHIMERAS
	##############################################################################
	if [[ "${remove_chimeras}" = "YES" ]] ; then
		echo $(date +%H:%M) 'Looking for chimeras in OTU fasta file using vsearch'
		source "${SCRIPT_DIR}"/chimera_check.sh "${OTU_fasta}"
		BLAST_INPUT="${chimera_free_fasta}"
	else
		BLAST_INPUT="${OTU_fasta}"
	fi

	################################################################################
	# BLAST CLUSTERS
	################################################################################
	blast_output="${DEREP_INPUT%/*}"/BOG18S.xml
	
	echo ":::::::::::::::::::::: BLAST ::::::::::::::::::::::::::::::::"
	echo $(date +%H:%M) " BLASTing..." 
	echo " Database: " "${BLAST_DB}"
	echo " Percent identity: " "${PERCENT_IDENTITY}"
	echo " Word size: " "${WORD_SIZE}"
	echo " E value: " "${EVALUE}"
	echo " Maximum target matches: " "${MAXIMUM_MATCHES}"
	echo " Culling Limit: " "${culling_limit}"
	echo " Output format: 5"
	
	blastn \
		-query "${BLAST_INPUT}" \
		-db "$BLAST_DB" \
		-num_threads "$n_cores" \
		-perc_identity "${PERCENT_IDENTITY}" \
		-word_size "${WORD_SIZE}" \
		-evalue "${EVALUE}" \
		-max_target_seqs "${MAXIMUM_MATCHES}" \
		-culling_limit "${culling_limit}" \
		-outfmt 5 \
		-out "${blast_output}"

	# check for blast output
	if [[ ! -s "${blast_output}"  ]]; then
		echo
		echo 'BLAST failed: the output file is empty or absent.'
	    echo 'File should be:' "${blast_output}"
		echo
	fi


else
################################################################################
# DON'T CONCATENATE SAMPLES
################################################################################

	################################################################################
	# PRIMER REMOVAL
	################################################################################
	echo $(date +%H:%M) "Removing primers..."
	for TAG_SEQ in $TAGS; do
		TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
		# Remove PRIMER1 from the beginning of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.
		cutadapt -g ^"${PRIMER1}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${TAG_DIR}"/2_notags.fasta > "${TAG_DIR}"/5_primerL1_removed.fasta
		cutadapt -g ^"${PRIMER2}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${TAG_DIR}"/2_notags.fasta > "${TAG_DIR}"/5_primerL2_removed.fasta
		# Remove the primer on the other end of the reads by reverse-complementing the files and then trimming PRIMER1 and PRIMER2 from the left side.
		# NOTE cutadapt1.7 will account for anchoring these to the end of the read with $
		seqtk seq -r "${TAG_DIR}"/5_primerL1_removed.fasta | cutadapt -g ^"${PRIMER2}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${TAG_DIR}"/6_primerR1_removed.fasta
		seqtk seq -r "${TAG_DIR}"/5_primerL2_removed.fasta | cutadapt -g ^"${PRIMER1}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${TAG_DIR}"/6_primerR2_removed.fasta
		seqtk seq -r "${TAG_DIR}"/6_primerR1_removed.fasta > "${TAG_DIR}"/6_primerR1_removedRC.fasta
		cat "${TAG_DIR}"/6_primerR1_removedRC.fasta "${TAG_DIR}"/6_primerR2_removed.fasta > "${TAG_DIR}"/7_no_primers.fasta
	done

	################################################################################
	# CONSOLIDATE IDENTICAL SEQUENCES
	################################################################################
	for TAG_SEQ in $TAGS; do
		TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"

		DEREP_INPUT="${TAG_DIR}"/7_no_primers.fasta

		# usearch -derep_fulllength "${DEREP_INPUT}" -sizeout -strand both -uc "${TAG_DIR}"/derep.uc -output "${TAG_DIR}"/7_derep.fasta
		echo $(date +%H:%M) "Consolidating identical sequences..."
		python "$SCRIPT_DIR/dereplication/dereplicate_fasta.py" "${DEREP_INPUT}"

		# REMOVE SINGLETONS
		# usearch -sortbysize "${TAG_DIR}"/7_derep.fasta -minsize 2 -sizein -sizeout -output "${TAG_DIR}"/8_nosingle.fasta
		# COUNT DUPLICATES PER READ, REMOVE SINGLETONS
		awk -F';' '{ if (NF > 2) print NF-1 ";" $0 }' "${DEREP_INPUT}".derep | sort -nr | awk -F';' '{ print ">DUP_" NR ";" $0}' > ${DEREP_INPUT%/*}/nosingle.txt

		# count the duplicates
		awk 'BEGIN { FS ="_tag_'${TAG_SEQ}'" } { print NF -1 }' "${DEREP_INPUT%/*}"/nosingle.txt > ${DEREP_INPUT%/*}/"${TAG_SEQ}".dup
#		awk 'BEGIN { FS ="_tag_'${TAG_SEQ}'" } { print NF -1 }' "${DEREP_INPUT%/*}"/nosingle.txt > ${DEREP_INPUT%/*}/"${TAG_SEQ}".dup

		# Write fasta file in order to blast sequences
		awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/no_duplicates.fasta

		# CLUSTER SEQUENCES
		if [ "$CLUSTER_OTUS" = "NO" ]; then
			BLAST_INPUT=${DEREP_INPUT%/*}/no_duplicates.fasta
		else
			CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
			UPARSE_OUT="${DEREP_INPUT%/*}"/OTU_uparse.txt
			usearch -cluster_otus "${DEREP_INPUT%/*}"/nosingle.txt -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${TAG_DIR}"/9_OTUs.fasta -uparseout "${UPARSE_OUT}"
			BLAST_INPUT="${TAG_DIR}"/9_OTUs.fasta
		fi

		# BLAST CLUSTERS
		blastn -query "${BLAST_INPUT}" -db "$BLAST_DB" -num_threads "$n_cores" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${TAG_DIR}"/10_BLASTed.xml
	done
fi



################################################################################
# TAXONOMIC ANNOTATION
################################################################################
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	DIRECTORIES="${DEREP_INPUT%/*}"
else
	DIRECTORIES=$( find "${ANALYSIS_DIR}"/demultiplexed -type d -d 1 )
fi

for DIR in "$DIRECTORIES"; do

	# Some POTENTIAL OPTIONS FOR MEGAN EXPORT:
	# {readname_taxonname|readname_taxonid|readname_taxonpath|readname_matches|taxonname_count|taxonpath_count|taxonid_count|taxonname_readname|taxonpath_readname|taxonid_readname}
	# PERFORM COMMON ANCESTOR GROUPING IN MEGAN

		# check for blast output
		if [[ -s "${blast_output}"  ]]; then

		echo $(date +%H:%M) "${blast_output}" ',BLAST output found; proceeding to MEGAN.  Taxonomic ranks exported to files.'
		# Specify paths to megan-related files
		BLAST_XML="${DIR}"/BOG18S.xml
		MEGAN_COMMAND_FILE="${DIR}"/megan_commands.txt
		MEGAN_RMA_FILE="${DIR}"/BOG18S.rma
		MEGAN_SHELL_SCRIPT="${DIR}"/megan_script.sh

		echo "import blastfile='${BLAST_XML}' meganFile='${MEGAN_RMA_FILE}' \
minScore=${MINIMUM_SCORE} \
maxExpected=${MAX_EXPECTED} \
topPercent=${TOP_PERCENT} \
minSupportPercent=${MINIMUM_SUPPORT_PERCENT} \
minSupport=${MINIMUM_SUPPORT} \
minComplexity=${MINIMUM_COMPLEXITY} \
lcapercent=${LCA_PERCENT};" > "${MEGAN_COMMAND_FILE}"
		echo "update;" >> "${MEGAN_COMMAND_FILE}"
		echo "collapse rank='$COLLAPSE_RANK1';" >> "${MEGAN_COMMAND_FILE}"
		echo "update;" >> "${MEGAN_COMMAND_FILE}"
		echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
		echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_${COLLAPSE_RANK1}.csv';" >> "${MEGAN_COMMAND_FILE}"
		if [ "$PERFORM_SECONDARY_MEGAN" = "YES" ]; then
			echo "collapse rank='$COLLAPSE_RANK2';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_${COLLAPSE_RANK2}.csv';" >> "${MEGAN_COMMAND_FILE}"
			echo "collapse rank='$COLLAPSE_RANK3';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_${COLLAPSE_RANK3}.csv';" >> "${MEGAN_COMMAND_FILE}"
			echo "collapse rank='SuperKingdom';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_kingdom.csv';" >> "${MEGAN_COMMAND_FILE}"
			echo "collapse rank='Order';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_order.csv';" >> "${MEGAN_COMMAND_FILE}"
			echo "collapse rank='class';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_class.csv';" >> "${MEGAN_COMMAND_FILE}"
			echo "collapse rank='Species';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file='${DIR}/BOG18S_meganout_species.csv';" >> "${MEGAN_COMMAND_FILE}"
			# echo "collapse rank='Phylum';" >> "${MEGAN_COMMAND_FILE}"
			# echo "update;" >> "${MEGAN_COMMAND_FILE}"
			# echo "select nodes=leaves;" >> "${MEGAN_COMMAND_FILE}"
			# echo "export what=biome data=Taxonomy file='${DIR}/BOG18S.biom';" >> "${MEGAN_COMMAND_FILE}"
			# echo "export what=CSV format=readname_taxonpath separator=comma file='${DIR}/BOG18S_tpath_species.csv';" >> "${MEGAN_COMMAND_FILE}"
		fi
		echo "quit;" >> "${MEGAN_COMMAND_FILE}"

		echo "#!/bin/bash" > "$MEGAN_SHELL_SCRIPT"
		echo "cd "${megan_exec%/*}"" >> "$MEGAN_SHELL_SCRIPT"
		echo "./"${megan_exec##*/}" -g -E -c ${DIR}/megan_commands.txt" >> "$MEGAN_SHELL_SCRIPT"

		# Run MEGAN
		sh "${DIR}"/megan_script.sh
#export what=CSV format=taxonname_count separator=comma counts=summarized file='C:\Users\reiko\Documents\reiko_BLASTed-5-ex.txt';
		# Modify the MEGAN output so that it is a standard CSV file with clusterID, N_reads, and Taxon
		echo '**sed process on megan, remove size= text'
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_${COLLAPSE_RANK1}.csv >"${DIR}"/BOG18S_meganout_${COLLAPSE_RANK1}_mod.csv
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_${COLLAPSE_RANK2}.csv >"${DIR}"/BOG18S_meganout_${COLLAPSE_RANK2}_mod.csv
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_${COLLAPSE_RANK3}.csv >"${DIR}"/BOG18S_meganout_${COLLAPSE_RANK3}_mod.csv
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_kingdom.csv >"${DIR}"/BOG18S_meganout_Kingdom_mod.csv
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_species.csv >"${DIR}"/BOG18S_meganout_Species_mod.csv
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_class.csv >"${DIR}"/BOG18S_meganout_Class_mod.csv
		sed 's|;size=|,|' <"${DIR}"/BOG18S_meganout_order.csv >"${DIR}"/BOG18S_meganout_Order_mod.csv
		rm "${DIR}"/BOG18S_meganout_${COLLAPSE_RANK1}.csv 
		# Run the R script, passing the current tag directory as the directory to which R will "setwd()"
		echo '**Rscript megan_plotter_reiko.R Based on MEGAN output.' "${DIR}"
		#Rscript "$SCRIPT_DIR/megan_plotter.R" "${DIR}"
		Rscript "$MY_SCRIPT_DIR"/megan_plotter_reiko.R "${DIR}"
	else
		echo
		echo 'BLAST failed: the output file is empty or absent.'
		echo 'File should be:' "${blast_output}"
		echo
	fi

done


################################################################################
# PRELIMINARY ANALYSES
################################################################################
# Once you have a final CSV file of the number of occurences of each OTU in each sample, run some preliminary analyses in R
# TODO rename preliminary to OTU analyses; move analysis script to OTU analysis directory
OUTPUT_PDF="${ANALYSIS_DIR}"/analysis_results_"${START_TIME}".pdf

echo $(date +%H:%M) "passing args to analyses_prelim.R for preliminary analysis..."
echo Rscript "$MY_SCRIPT_DIR"/analyses_prelim.R  "${OUTPUT_PDF}" "${OTU_table}" "${SEQUENCING_METADATA}" "${LIBRARY_COLUMN_NAME}" "${TAG_COLUMN_NAME}" "${ColumnName_SampleName}" "${ColumnName_SampleType}"
Rscript "$MY_SCRIPT_DIR"/analyses_prelim.R "${OUTPUT_PDF}" "${OTU_table}" "${SEQUENCING_METADATA}" "${LIBRARY_COLUMN_NAME}" "${TAG_COLUMN_NAME}" "${ColumnName_SampleName}" "${ColumnName_SampleType}"
echo $(date +%H:%M) "completed R for preliminary analysis..."

# EMPTY PDFs are 3829 bytes
minimumsize=4000
size_PDF=$(wc -c <"${OUTPUT_PDF}")
if [ "${size_PDF}" -lt "${minimumsize}" ]; then
    echo 'There was a problem generating the PDF.'
else
	REMOTE_PDF="${OUTPUT_PDF_DIR}"/analysis_results_"${START_TIME}".pdf
	cp "${OUTPUT_PDF}" "${REMOTE_PDF}"
	echo $(date +%H:%M) "analysis_results pdf written; R for preliminary analysis..."

fi


if [ "$PERFORM_CLEANUP" = "YES" ]; then
	echo $(date +%H:%M) "Compressing fasta, fastq, and xml files..."
	find "${ANALYSIS_DIR}" -type f -name '1_merged.assembled.fastq' -exec rm "{}" \; 
	find "${ANALYSIS_DIR}" -type f -name '2_filtered.fasta' -exec rm "{}" \; 
	find "${ANALYSIS_DIR}" -type f -name '2_filtered_renamed.fasta' -exec rm "{}" \; 
	find "${ANALYSIS_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
	find "${ANALYSIS_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
	find "${ANALYSIS_DIR}" -type f -name '*.xml' -exec ${ZIPPER} "{}" \;
	if [ "$DEBUG1" = "NO" ]; then
		rm "${ANALYSIS_DIR}"/all_lib/nosingle.txt
		rm "${ANALYSIS_DIR}"/all_lib/1_demult_concat.fasta.derep
		rm "${ANALYSIS_DIR}"/all_lib/1_demult_concat.fasta.gz
		rm -r "${ANALYSIS_DIR}"/all_lib/dup_temp
	fi
	echo $(date +%H:%M) "Cleanup performed."
else
	echo $(date +%H:%M) "Cleanup not performed."
fi

# ::::::::::::::::::: CREATE BIOM FILE ::::::::::::::::::::::::

# create observational metadata for biom file created from OTU_table.txt, from analyses_prelim.R
# The following code used the taxonomic path from MEGAN but for some reason it is not getting the complete path for each OTU_ID
		# sed -i -e 's|root;cellular organisms;Eukaryota;||' "${DIR}"/BOG18S_tpath_species.csv 
		# sed -i -e 's|root;cellular organisms;Bacteria;||' "${DIR}"/BOG18S_tpath_species.csv 
		# sed -i -e 's|;size=|,|' "${DIR}"/BOG18S_tpath_species.csv 
		# sed -i -e 's|DUP|OTU|' "${DIR}"/BOG18S_tpath_species.csv 
		# echo "#OTU_ID	size	taxonomy" > "${DIR}"/BOG18S_obs_md.txt
		# sed  's|,|\t|g' < "${DIR}"/BOG18S_tpath_species.csv >> "${DIR}"/BOG18S_obs_md.txt
		echo "Make observational metadata, concatenate taxa ranks. python $MY_SCRIPT_DIR/Merge_taxa_ranks.py" "${ANALYSIS_DIR}"/all_lib
		python "$MY_SCRIPT_DIR"/Merge_taxa_ranks.py "${ANALYSIS_DIR}"/all_lib
		sed -i -e 's/DUP/OTU/' -e 's/,/\t/' "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.csv
		sed -i -e 's/,/;/g' "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.csv
		echo "#OTU_ID	taxonomy" > "${ANALYSIS_DIR}"/all_lib/BOG18S_obs_md.txt
		cat "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.csv >> "${ANALYSIS_DIR}"/all_lib/BOG18S_obs_md.txt
		echo $(date +%H:%M) "Create new otu tables with taxonomic ranks"
		echo " Merge taxonomy to OTU table.  Does not include no hits, not assigned.  Output is OTU_table_taxa.txt
		python merge_otu_taxa.py " "${ANALYSIS_DIR}"/all_lib
		
# Stuff for phinch.  Add prefixes for the taxonomic ranks.  Then merge by OTU id and concatenate the ranks
		echo $(date +%H:%M) "Create Phinch style taxonomic names, ie. p__ prefix for phylum"
		sed -i -e 's|,"|,"f__|' "${DIR}"/BOG18S_meganout_${COLLAPSE_RANK1}_mod.csv
		sed -i -e 's|,"|,"g__|' "${DIR}"/BOG18S_meganout_${COLLAPSE_RANK2}_mod.csv
		sed -i -e 's|,"|,"p__|' "${DIR}"/BOG18S_meganout_${COLLAPSE_RANK3}_mod.csv
		sed -i -e 's|,"|,"s__|' "${DIR}"/BOG18S_meganout_Species_mod.csv
		sed -i -e 's|,"|,"c__|' "${DIR}"/BOG18S_meganout_Class_mod.csv
		sed -i -e 's|,"|,"o__|' "${DIR}"/BOG18S_meganout_Order_mod.csv
		sed -i -e 's|,"|,"k__|' "${DIR}"/BOG18S_meganout_Kingdom_mod.csv
		echo 'Make phinch observational metadata, python "$MY_SCRIPT_DIR"/Merge_taxa_ranks.py' "${ANALYSIS_DIR}"/all_lib
		python "$MY_SCRIPT_DIR"/Merge_taxa_ranks.py "${ANALYSIS_DIR}"/all_lib
		sed -i -e 's/DUP/OTU/' -e 's/,/\t/' "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.csv
		sed -i -e 's/,/;/g' "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.csv
		echo "#OTU_ID	taxonomy" > "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.txt
		cat "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.csv >> "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.txt
		
echo $(date +%H:%M) "Create new otu tables with Phinch style ranks"
echo " Merge taxonomy to OTU table.  Does not include no hits, not assigned.  Output is OTU_table_taxa_phinch.txt with Phinch prefixes denoting taxonomic rank.  i.e. p__Proteobacteria"
python merge_otu_taxa2.py "${ANALYSIS_DIR}"/all_lib
# python adds quotes even though the string is not quoted.  To get it to not put quotes in is a pain.  This is totally stupid
sed -i 's/"//g' "${ANALYSIS_DIR}"/all_lib/OTU_table_taxa_phinch.txt
python merge_otu_taxa_all.py "${ANALYSIS_DIR}"/all_lib

echo $(date +%H:%M) "Create biom file BOG18S.biom" 
biom convert -i "${ANALYSIS_DIR}"/all_lib/OTUs_swarm/OTU_table.txt -o "${ANALYSIS_DIR}"/all_lib/BOG18S_json.biom --table-type="OTU table" --to-json
biom convert -i "${ANALYSIS_DIR}"/all_lib/OTUs_swarm/OTU_table.txt -o "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5.biom --table-type="OTU table" --to-hdf5
biom add-metadata -i "${ANALYSIS_DIR}"/all_lib/BOG18S_json.biom -o "${ANALYSIS_DIR}"/all_lib/BOG18S_json_md.biom --sample-metadata-fp "$MY_SCRIPT_DIR"/BOG_kelp18S/BOG18S_biom_metadata_all.txt
biom add-metadata -i "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5.biom -o "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5_md.biom --sample-metadata-fp "$MY_SCRIPT_DIR"/BOG_kelp18S/BOG18S_biom_metadata_all.txt

biom add-metadata -i "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5_md.biom -o "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5_md_obs.biom --observation-metadata-fp "${DIR}"/BOG18S_obs_md.txt --int-fields size --sc-separated taxonomy
biom add-metadata -i "${ANALYSIS_DIR}"/all_lib/BOG18S_json_md.biom -o "${ANALYSIS_DIR}"/all_lib/BOG18S_json_md_obs.biom --observation-metadata-fp "${ANALYSIS_DIR}"/all_lib/phinch_obs_md.txt  
# ::::::::::::::::::::::: BIOM FILE :::::::::::::::::::::::::::

echo "Create a PEAR report with the names of the read files and the assembled, discarded, and not assembled read numbers.  The report also includes the PEAR parameters"
echo "mk_PEAR_report.py " "${ANALYSIS_DIR}"
python mk_PEAR_report.py "${ANALYSIS_DIR}"

FINISH_TIME=$(date +%Y%m%d_%H%M)
if [ "$NOTIFY_EMAIL" = "YES" ]; then
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME | mail -s "banzai is finished" "${EMAIL_ADDRESS}"
else
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME
fi
echo 'Data are in ' "${ANALYSIS_DIR}"
# remove some of the repetitive lines, read percents from kmers, fastq etc.  makes it difficult to read logfile.
# grep -v flag means output everything except pattern, -f means pattern is in text file, i.e. grep_patterns.txt
grep -v -f "$MY_SCRIPT_DIR"/grep_patterns.txt "${ANALYSIS_DIR}"/logfile.txt > "${ANALYSIS_DIR}"/logfile2.txt
# Quick scan of logfile2.txt for errors, error line and 3 lines before, -B 3
grep -n -B 3 "error" "${ANALYSIS_DIR}"/logfile2.txt >"${ANALYSIS_DIR}"/errors.txt
# echo -e '\n'$(date +%H:%M)'\tAll finished! Why not treat yourself to a...\n'
# echo
# echo -e '\t~~~ MAI TAI ~~~'
# echo -e '\t2 oz\taged rum'
# echo -e '\t0.75 oz\tfresh squeezed lime juice'
# echo -e '\t0.5 oz\torgeat'
# echo -e '\t0.5 oz\ttriple sec'
# echo -e '\t0.25 oz\tsimple syrup'
# echo -e '\tShake, strain, and enjoy!' '\xf0\x9f\x8d\xb9\x0a''\n'
