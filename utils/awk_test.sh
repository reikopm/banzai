
TAG_DIR="/home/mbonteam/MBARI/reiko/processed/Analysis_20160201_1602/lib2/demultiplexed/tag_AACAAC"
demult_file_L="${TAG_DIR}"/1_tagL_removed_AWK.fasta
demult_file_R="${TAG_DIR}"/2_notags_AWK.fasta
DEMULTIPLEX_INPUT="/home/mbonteam/MBARI/reiko/processed/Analysis_20160201_1602/lib2/2_filtered_renamed.fasta"
TAG_SEQ="AACAAC"
		 awk 'gsub(/^.{3}'"$TAG_SEQ"'/,"") {
			 if (a && a !~ /^.{3}'"$TAG_SEQ"'/)
				# print a;
			 print
		 } {a=$0}' "${DEMULTIPLEX_INPUT}" > "${demult_file_L}"
		
				# Right side tag
		TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
		echo 'TAG_RC' $TAG_RC
		 awk 'gsub(/'"$TAG_RC"'.{3}$/,"") {
			if (a && a !~ /'"$TAG_RC"'.{3}$/)
				print a "tag_""'"$TAG_SEQ"'";
			print
		} {a = $0}' "${TAG_DIR}"/1_tagL_removed_AWK.fasta > "${demult_file_R}"
