# Some samples had the incorrect name and need to be fixed before moving forward, Run these two commands in FASTQ directory to find and fix the names 
for f in *_Cxt_R1_trimmed.fastq.gz; do mv "$f" "${f/_Cxt/-Cxt}"; done
for f in *_Cxt_R2_trimmed.fastq.gz; do mv "$f" "${f/_Cxt/-Cxt}"; done


