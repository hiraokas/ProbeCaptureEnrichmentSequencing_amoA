#!/bin/bash
#$ -m be
#$ -cwd
#$ -pe threads 6

function usage() {
    cat <<'EOF'
========================================================================================================================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    By: Satoshi Hiraoka
    Created: <2017
    History:  20210825
    History:  20221001  (Now .gz file is acceptable for mmseq2)
    History:  20221031
	- This script is for sequence clustering (e.g., CDS, OTU, etc...)
	- CD-HIT: http://weizhongli-lab.org/cd-hit/wiki/doku.php?id=cd-hit_user_guide#hierarchically_clustering
Usage:
    geneclustering.sh [-a] [-o] [-i] [-t]
Required:
    Require: cd-hit, perl or MMseq2
Mandatory:
    -s      using tool [Nucleotide: cd-hit, vsearch, MMseq2, MMseq2LC]
                       [Protein:    cd-hit         , MMseq2, MMseq2LC]
    -i      Input filename (faa:protein, fna:nucleotide. .gz file is acceptable)
    -t      clustering threshold (98, 97, 95, 90, 85, 90-80, 90-75, 90-70, 90-60, and 90-60-30 for CD-hit)
                                 (0.0-1.0 for vsearch)
                                 (0.0-1.0 for MMseq2)
Options:
    -o      Output directory (default: ../clustering)
    -c      coverage cutoff  (=0.8) [MMseq2]
    -T      Pararel threds   (=6)
    -M      Memory (=5000M) (CD-HIT)
    -h      Print this (help print)
    -f      Force running
    -F      Force running (strong)
TIPS:
    ./clustering_seq.sh -a vsearch  -i ../genecall/.fna -t 98
    ./qsub_DDBJ.sh medium 20 20 48 ./clustering_seq.sh -s MMseq2LC -i ../gene/PublicDatabase/DSSMv0.2_with_nrandGMGC10.95.faa.gz -t 0.95 -c 0.90 -T 16
========================================================================================================================================
EOF
	return 0
}


WORKNAME="Sequence clustering"
PSICDHIT="psi-cd-hit.pl"
CDHIT="    ${HOME}/software/cd-hit-v4.8.1-2019-0228/cd-hit"
CDHIT_EST="${HOME}/software/cd-hit-v4.8.1-2019-0228/cd-hit"
CDHIT_REV="${HOME}/software/cd-hit-v4.8.1-2019-0228/clstr_rev.pl"
mmseqs="   ${HOME}/software/mmseqs/bin/mmseqs"
VSEARCH="  ${HOME}/software/vsearch/bin/vsearch"

usage_exit() {
	usage
	exit 1
}

##init
THREAD=6
MEMORY=5000
FORCE=0
COVER=0.8
start_time=`date +%s`
INPUT_FILE=""
TOOL=""
OUT_DIR="../clustering"

#==========================================================================================
#http://qiita.com/b4b4r07/items/dcd6be0bb9c9185475bb
while getopts T:M:i:o:t:s:hfFc: OPT
do
	case $OPT in
		T)	THREAD=$OPTARG
			;;
		M)	MEMORY=$OPTARG
			;;
		i)	INPUT_FILE=$OPTARG
			FILENAME=${INPUT_FILE##*/}
			BASE_FILENAME=${FILENAME%.*} 
			EXT_FILENAME=` echo ${FILENAME}| sed -e "s/.gz$//g" | rev | cut -f1 -d "." | rev`
            ;;
        s)  TOOL=$OPTARG
			;; 
		o)	OUT_DIR=$OPTARG
			;;
		t)	CLUSTERING_THRESHOLD=$OPTARG
			;;
		c)	COVER=$OPTARG
			;;
		h)	usage_exit
			;;
		v)	usage_exit
			;;
		f)  FORCE=1
			;;
		F)  FORCE=2
			;;
		\?)	usage_exit
			;;
	esac
done

TH1=0
OT1=0
ON1=0
TH2=0
OT2=0
ON2=0
TH3=0
OT3=0
ON3=0

flag=0
if [ -z ${OUT_DIR} ]; then
	echo "[**Error**]  No output directory provided."
	usage_exit
fi
if [ -z ${INPUT_FILE} ]; then
	echo "[**Error**]  No input file."
	usage_exit
fi
if [ ! -e ${INPUT_FILE} ]; then
	echo "[**Error**]  input file not exist."
	usage_exit
fi
if [ -z ${CLUSTERING_THRESHOLD} ]; then
	echo "[**Error**]  No clustering threshold input."
	usage_exit
fi

#echo ${EXT_FILENAME}
if [ ${EXT_FILENAME} = "faa" ]; then
	echo "  faa file input: treated as [protein sequence]."
	BASE_FILENAME=${BASE_FILENAME}_P
elif [ ${EXT_FILENAME} = "fna" -o  ${EXT_FILENAME} = "fa" -o  ${EXT_FILENAME} = "fasta" -o  ${EXT_FILENAME} = "fastq" ]; then
	echo "  fna file input: treated as [nucleotide sequence]."
	BASE_FILENAME=${BASE_FILENAME}_N
else
	echo "[**Error**]  unidentified sequence type: ${EXT_FILENAME}"
	usage_exit
fi

case ${CLUSTERING_THRESHOLD} in
	"98")	TH1=98
			OT1=0.98
			ON1=5
			;;
	"97")	TH1=97
			OT1=0.97
			ON1=5
			;;
	"95")	TH1=95
			OT1=0.95
			ON1=5
			;;
	"90")	TH1=90
			OT1=0.9
			ON1=5
			;;
	"85")	TH1=85
			OT1=0.85
			ON1=5
			;;
	#"80")		TH1=80
	#		OT1=0.8
	#		ON1=5
	#		;;
	#"75")		TH1=75
	#		OT1=0.75
	#		ON1=5
	#		;;
	"90-80")TH1=90
			OT1=0.9
			ON1=5
			TH2=80
			OT2=0.8
			ON2=5
			;;
	"80")	TH1=80
			OT1=0.8
			ON1=5
			;;
	"90-75")TH1=90
			OT1=0.9
			ON1=5
			TH2=75
			OT2=0.75
			ON2=5
			;;
	"90-70")TH1=90
			OT1=0.9
			ON1=5
			TH2=70
			OT2=0.7
			ON2=5
			;;
	"70")	TH1=70
			OT1=0.7
			ON1=5
			;;
	"90-60")TH1=90
			OT1=0.9
			ON1=5
			TH2=60
			OT2=0.6
			ON2=4
			;;
	"60")	TH1=60
			OT1=0.6
			ON1=5
			;;
	"90-60-30")	TH1=90
			OT1=0.9
			ON1=5
			TH2=60
			OT2=0.6
			ON2=4
			TH3=30
			OT3=0.3
			ON3=3
			;;
	\?)		#echo "[**Error**]  unidentified threshold of sequence simirarity for clustering."
			#usage_exit
			OT1=${CLUSTERING_THRESHOLD}
			;;
esac

#make directory
if [ ! -e ${OUT_DIR} ]; then
    mkdir ${OUT_DIR}
fi

#run
echo "${WORKNAME} start..." 

#==========================================================================================
#cd-hit
#==========================================================================================
if [ ${TOOL} = "cd-hit" ]; then
	echo "  cluster using cd-hit."
	output_filename1="${OUT_DIR}/${BASE_FILENAME}_"${TH1}".cdhit"
	output_filename2="${OUT_DIR}/${BASE_FILENAME}_temp"${TH2}".cdhit"
	output_filename2_5="${OUT_DIR}/${BASE_FILENAME}_"${TH2}"-"${TH1}".cdhit"
	output_filename3="${OUT_DIR}/${BASE_FILENAME}_temp"${TH3}".cdhit"
	output_filename3_5="${OUT_DIR}/${BASE_FILENAME}_"${TH3}"-"${TH2}"-"${TH1}".cdhit"

	echo "  make: ${output_filename1}.fna"
	if [ ${FORCE} -lt 2 -a -e ${output_filename1}.clstr ]; then
		echo "**Skipped 1st step.**"
	else
		echo "1st step..."
		command="${CDHIT} -i ${INPUT_FILE} -o ${output_filename1} -c ${OT1} -M ${MEMORY} -T ${THREAD} -n ${ON1} -d 0"
		echo "  ${command}"
		${command}
		mv ${output_filename1} ${output_filename1}.fna
	fi

	if [ ${TH2} -ne 0 ]; then
		echo "  make: ${output_filename2_5}.fna"
		if [ ${FORCE} -ne 1 -a -e ${output_filename2_5}.clstr ]; then
			echo "**Skipped 2nd step."
		else
			echo "2nd step..."
			command="${CDHIT} -i ${output_filename1}.fna -o ${output_filename2} -c ${OT2} -M ${MEMORY} -T ${THREAD} -n ${ON2} -d 0"
			echo "  ${command}"
			${command}
			echo "marge step (1 & 2)"
			perl ${CDHIT_REV} ${output_filename1}.clstr ${output_filename2}.clstr > ${output_filename2_5}.clstr

			mv ${output_filename2} ${output_filename2_5}.fna
		fi
	fi

	if [ ${TH3} -ne  0 ]; then
		if [ ${FORCE} -ne 1 -a -e ${output_filename3_5}.clstr ]; then
			echo "**Skipped 3rd step."
		else
			echo "3rd step..."
			#command="${PSICDHIT} -i ${output_filename2} -o ${output_filename3} -c 0.3 -core 3"
			command="${CDHIT} -i ${output_filename2_5}.fna -o ${output_filename3} -c ${OT3} -M ${MEMORY} -T ${THREAD} -n ${ON3} -d 0"
			echo "  ${command}"
			${command}
			
			echo "merge step (12 & 3)"
			perl ${CDHIT_REV} ${output_filename2_5}.clstr ${output_filename3}.clstr > ${output_filename3_5}.clstr
			mv ${output_filename3} ${output_filename3_5}.fna
		fi
	fi

#==========================================================================================
#vsearch
#==========================================================================================
elif [ ${TOOL} = "vsearch" ]; then

	if [ ${EXT_FILENAME} = "faa" ]; then
		echo "vsearch does not support protein sequence. Exit."
		exit
	fi

	echo "  cluster using vsearch."
	output_filename1="         ${OUT_DIR}/${BASE_FILENAME}_${OT1}.vsearch"
	output_filename_singleton="${OUT_DIR}/nosinglton_${BASE_FILENAME}_${OT1}.fna"

	if [ ${FORCE} -ne 1 -a -e ${output_filename1}.vserch_out ]; then
		echo "**Skipped vsearch"
	else
		echo "vsearch"
		command="${VSEARCH} --cluster_fast ${INPUT_FILE} --id ${OT1} --centroids ${output_filename1}.fna --threads ${THREAD} --msaout ${output_filename1}.msaout --sizeout"
		echo "  ${command}"
		${command}

		cat ${output_filename1}.msaout |grep ">"|cut -b 2-|sed -e "s/consensus//g" >> ${output_filename1}.vserch_out
		rm  ${output_filename1}.msaout
	fi

	if [ ${FORCE} -ne 1 -a -e ${output_filename_singleton} ]; then
		echo "**Skipped vsearch"
	else
		echo "remove singltons"
		seqkit grep -nrp "size=1\$" ${output_filename1}.fna -v > ${output_filename_singleton}
	fi


#==========================================================================================
#MMseq2 easy
#==========================================================================================
elif [ ${TOOL} = "MMseq2" ]; then
    echo "  cluster using MMseq2."
    DBname="      ${OUT_DIR}/${BASE_FILENAME}.MMseq2DB"
    outputPrefix="${OUT_DIR}/${BASE_FILENAME}_${CLUSTERING_THRESHOLD}-${COVER}.MMseq2"
    tmpdir="      ${OUT_DIR}/tmp"

	${mmseqs} easy-cluster ${INPUT_FILE} ${DBname} ${tmpdir} --min-seq-id ${CLUSTERING_THRESHOLD} -c ${COVER} --cov-mode 0 --threads ${THREAD}

    #cleanup
    mv ${DBname}_rep_seq.fasta ${outputPrefix}.fasta
    mv ${DBname}_cluster.tsv   ${outputPrefix}.tsv
    rm ${DBname}_all_seqs.fasta
    rm ${tmpdir} -r


#==========================================================================================
# #MMseq2 linclust
#==========================================================================================
 elif [ ${TOOL} = "MMseq2LC" ]; then
     echo "  cluster using MMseq2. (linclust)"
     DBname="      ${OUT_DIR}/${BASE_FILENAME}.MMseq2DB"
     outputPrefix="${OUT_DIR}/${BASE_FILENAME}_${CLUSTERING_THRESHOLD}.MMseq2LC"
     tmpdir="      ${OUT_DIR}/tmp"

	${mmseqs} easy-linclust  ${INPUT_FILE} ${DBname} ${tmpdir} --min-seq-id ${CLUSTERING_THRESHOLD} -c ${COVER} --cov-mode 0 --threads ${THREAD} --split-memory-limit 320G

    #cleanup
    mv ${DBname}_all_seq.fasta ${outputPrefix}.fasta.all.fasta
    mv ${DBname}_rep_seq.fasta ${outputPrefix}.fasta
    mv ${DBname}_cluster.tsv   ${outputPrefix}.cluster.tsv

    echo "Making ClusterID-count table. (This will takes a few minutes.) : ${outputPrefix}.clusterCount.tsv"
    cat ${outputPrefix}.cluster.tsv |  cut -f1 |sort | uniq -c | sed 's/^ *\| *$//' | awk -F " " '{print $2" "$1}' > ${outputPrefix}.clusterCount.tsv

    #rm ${DBname}_all_seqs.fasta
    rm ${tmpdir} -r
    rm ${DBname}_all_seqs.fasta  


#==========================================================================================
# else
#==========================================================================================
else
	echo "[**Error**]  undefined clustering tool: ${TOOL}"
	usage_exit
fi


#==========================================================================================
# log
#==========================================================================================
echo "${WORKNAME} done."

if [ $? -eq  0 ]; then
	echo "Success."
else
	echo "Failed."
fi

end_time=`date +%s`
PT=$((end_time - start_time))
H=`expr ${PT} / 3600`
PT=`expr ${PT} % 3600`
M=`expr ${PT} / 60`
S=`expr ${PT} % 60`
echo "Run Time: ${H}h ${M}m ${S}s"
echo ""

exit 1
