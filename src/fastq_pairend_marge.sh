#!/bin/sh
#$ -m be
#$ -cwd
#$ -pe threads 1

function usage() {
    cat <<'EOF'
===================================================================================================================================================================================================
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Created:  <20181103
    History:   20220418 (support .gz)
    History:   20230119
    - Merging illumina paired-end sequence into one singlend reads.
    - threads=6
    - .fastq.gz will be acceptable
Usage:
    this.sh mode[dirNoComp, dir, mate-pair]
        paired-end mode:   this.sh  paired-end  paired-end-1                                [threads=6]
        dirNoComp mode:    this.sh  dir         target_dir(after low-complex-read removal)  [threads=6]
        dir mode:          this.sh  dir         target_dir(after adapter trim)              [threads=6] 
Output:
    outputdir="../data/16_merged/"
Ext:
    ./fastq_pairend_marge.sh  dirNoComp ../data/15_nocomplex_upper100
    (./fastq_pairend_marge.sh dir       ../data/12_QC)
===================================================================================================================================================================================================
EOF
    return 0
}   

outputdir="../data/16_merged/"
flash="${HOME}/software/FLASH-1.2.11-Linux-x86_64/flash"
threads=6

if [ $# -lt 2 ]; then
    echo "Insufficient parameters"
    usage
    exit 1
fi

mode=${1}


if [ $# -gt 2 ]; then
    threads=${3}
    echo "threads: ${threads}"
fi

if [ ! -e ${outputdir} ]; then
    mkdir ${outputdir}
fi

#===================================================
if   [ ${mode} == "dirNoComp" ]; then
    targetdir=${2}

    for f in `ls ${targetdir}/*_1.fastq.gz|rev|cut -f2- -d "_"|rev|sort|uniq`; 
        do echo ${f}; 
        bn=`basename ${f}`
        option=""   
        echo "Prefix: ${bn}"
        if [ ! -e ${outputdir}/${bn}.fastq ]; then
            ${flash} -o ${bn} -d  ${outputdir} ${f}_1.fastq.gz ${f}_2.fastq.gz -t ${threads}
            rm ${outputdir}/${bn}.hist
            rm ${outputdir}/${bn}.histogram
            rm ${outputdir}/${bn}.notCombined_1.fastq        
            rm ${outputdir}/${bn}.notCombined_2.fastq  

            mv ${outputdir}/${bn}.extendedFrags.fastq ${outputdir}/${bn}.fastq
        else
            echo "@@@Skipped: "${bn}
        fi
    done


#===================================================
elif   [ ${mode} == "dir" ]; then
    targetdir=${2}

    for f in `ls ${targetdir}/*.fq|rev|cut -f2- -d "_"|rev|sort|uniq`; 
        do echo ${f}; 
        bn=`basename ${f}`
        option=""   
        echo "Prefix: ${bn}"
        if [ ! -e ${outputdir}/${bn}.extendedFrags.fastq ]; then
            ${flash} -o ${bn} -d  ${outputdir} ${f}_R1.fq ${f}_R2.fq -t ${threads}
            rm ${outputdir}/${bn}.hist
            rm ${outputdir}/${bn}.histogram
            rm ${outputdir}/${bn}.notCombined_1.fastq        
            rm ${outputdir}/${bn}.notCombined_2.fastq  
        else
            echo "@@@Skipped: "${bn}
        fi
    done


#===================================================
elif [ ${mode} == "paired-end" ]; then
    file_1=${2}
    file_2=`echo ${2} | sed -e "s/_1.fastq/_2.fastq/g"`

    echo "file1: ${file_1}"
    echo "file2: ${file_2}"

    bn=`basename ${file_1}`
    echo "Prefix: ${bn}"

    if [ ! -e ${outputdir}/${bn}.extendedFrags.fastq ]; then
        ${flash} -o ${bn} -d  ${outputdir} ${file_1} ${file_2} -t ${threads} -O
        rm ${outputdir}/${bn}.hist
        rm ${outputdir}/${bn}.histogram
        rm ${outputdir}/${bn}.notCombined_1.fastq        
        rm ${outputdir}/${bn}.notCombined_2.fastq   

        rename "_1.fastq.gz.extendedFrags" "" *
    else
        echo "@@@Skipped: ${bn}"
    fi
else
    echo "undefined mode: ${mode}"
    exit 0  
fi

echo "@@@Done."
