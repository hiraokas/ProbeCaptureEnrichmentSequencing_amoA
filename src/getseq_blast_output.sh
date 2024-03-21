
function usage() {
    cat <<EOF
========================================================================================================================================
Description:
    $(basename ${0}) is a tool for extracting specific sequences using blast output files (format=6)
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Created: 20190401
    History: 20190409
    History: 20210924
    History: 20220922 (accept for hmmer output)
    History: 20230317 (change output dir)
    History: 20230317
    - blast-like (blast and DIAMOND) and hmmer outputs are acceptable.
    - output dir is ../CandidateSeq_blast or ../CandidateSeq_domtblout
Usage:
    $(basename ${0}) sequence.faa output.blast 
TIPS:
    ./getseq_blast_output.sh ../genecall/megahit_all_20181117.faa  ../blast/blast_megahit_all_20181117.tsv
========================================================================================================================================
EOF
}

usage_exit() {
    usage
    exit 1
}

if [ $# -ne 2 ]; then
    usage
    exit 1
fi

input_seq_file=${1}
blast_output=${2}
echo ${blast_output}

if [ ! -e ${OutputDir} ]; then
    mkdir ${OutputDir}
fi

SeqFilename=${input_seq_file##*/}
SeqFilenameBase=${SeqFilename%.*}
SeqFileExt=${SeqFilename##*.}

BlastFilename=${blast_output##*/}
BlastFilenameBase=${BlastFilename%.*}
BlastFileExt=${SeqFilename##*.}

OutputDir="../CandidateSeq_${BlastFileExt}"  #CandidateSeq_blast / CandidateSeq_domtblout
mkdir ${OutputDir} -p


output_genelist=" ${OutputDir}/${BlastFilenameBase}_${SeqFilenameBase}.tsv"
output_geneseq="  ${OutputDir}/${BlastFilenameBase}_${SeqFilenameBase}.${SeqFileExt}"


echo "get blasted sequence IDs"
cat ${blast_output} | grep -v "^#" | sed -e 's/  */\t/g' | cut -f1  >${output_genelist}
echo "Output: ${output_genelist}"

echo "extract all blasted sequences"
seqkit grep -f ${output_genelist} ${input_seq_file} -o ${output_geneseq} #  --quiet
echo "Output: ${output_geneseq}"
echo ""

rm ${output_genelist}

#echo "filtered complete sequence"
#grep "partial=00" ${output_geneseq} | cut -f1 -d " " | cut -c2- > ${output_genelist2}

#echo "extract complete blasted sequences"
#seqkit grep -f ${output_genelist2} ${output_geneseq} -o ${output_geneseq2}
#echo "output: ${output_geneseq2}"

echo "Done."
