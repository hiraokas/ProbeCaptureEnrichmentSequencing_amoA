#!/bin/sh
#$ -m be
#$ -cwd
#$ -pe threads 6

function usage() {
    cat <<'EOF'
============================================================================================================================================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Created: ?????
    History: 20211201
    - using bioawk
Usage:
    conda activate py38
    this fastA/Q

Preparation:
    conda install bioawk
    - install clistats
    - https://github.com/dpmcmlxxvi/clistats
============================================================================================================================================================
EOF
    return 0
}   


if [ $# -lt 1 ]; then
    usage
    exit 1
fi


FILENAME=${1##*/}
BASE_FILENAME=${FILENAME%.*}

seqlen_tsv=../SummaryInformation/seqlen_${BASE_FILENAME}.tsv

bioawk -c fastx '{ print $name, length($seq) }' < ${1} > ${seqlen_tsv}
cat ${seqlen_tsv} |cut -f2|clistats
           