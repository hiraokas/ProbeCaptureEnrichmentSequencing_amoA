#!/bin/sh

function usage() {
    cat <<'EOF'
==================================================================================================================================================================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Created: 20180225
    History: 20210402
    History: 20220723 (add -Y option for retrieve unmapped reads)
    History: 20220924 (Return sorted BAM file as the last standerd output)
    History: 20230215 (add -I option for making index only)
    History: 20230215 (samtools update 1.16.1)
    History: 20240112
Require:
    -t        mapping tools [bowtie2, bwa, bwa_matepair, bowtie2strict
                             pbmm2, pbmm2CCS, pbmm2HiFi, pbmm2CCSHiFi, blasr, minimapONT, minimapHiFi]
                            (Abolish: pbalign, pbalignCCS)
    -b        base file     [fasta]
    -i        input file    [fasta/q (bowtie2)
                             fasta   (bwa, minimapONT, minimapHiFi)
                             bam     (pbmm, pbmm2CCS, pbmm2HiFi, blasr, pbalign, pbalignCCS)]
    -1, -2    pair-end file [fasta/q (bowtie2)]
Option
    -T        threads [=6] (fqrequently fails when with numbers of threads (e.g., >31))
    -X        Retrive mapped     reads accompany with mapping results.
    -Y        Retrive non-mapped reads accompany with mapping results.  (HiFi reads)
    -R        Retain unsorted bam file.
    -r        Reorder ('--reorder' for botwie2)
    -I        Only making index (bowtie2, pbmm)
Usage:
    conda activate py38
    this.sh -t bowtie2 -b assembled.fasta -i rawreads.fasta
    this.sh -t bowtie2 -b assembled.fasta -1 rawreads_R1.fasta -2 rawreads_R2.fasta

    - return: output sorted BAM file path. This will get using following command:
    - $ tail -n 1
TIPS:
    for f in ../assembly/SPAdes_*.fasta; do sample=`basename ${f} | cut -f2- -d "_"| rev | cut -f2- -d "."|rev`;./mapping.sh -t bowtie2 -b ${f} -i ../data/18_fasta/${sample}.fa -T 38; done

------------------------------------------------------------
 INSTALL TOOLS
------------------------------------------------------------
    - "samtool" can be installed by conda, but need library linking as follow:
    cd ~/anaconda3/envs/py38/lib
    ln -s libcrypto.so.1.1 libcrypto.so.1.0.0

    #pbmm2
    conda activatepy38
    conda install pbmm2 -c bioconda
==================================================================================================================================================================================
EOF
    return 0
}   

#================================================================================================
# Definition
#================================================================================================
THREADS=6
bowtie2="      ${HOME}/software/bowtie2-2.4.4-linux-x86_64/bowtie2"
bowtie2_build="${HOME}/software/bowtie2-2.4.4-linux-x86_64/bowtie2-build"
samtools="     ${HOME}/software/samtools-1.16.1/samtools"
bedtools="     ${HOME}/software/bedtools2/bin/bedtools"
bwa="bwa"
blasr="        ${HOME}/anaconda3/envs/py38/bin/blasr"
pitchfork_env="${HOME}/software/pitchfork/deployment/setup-env.sh"
pbmm2="        ${HOME}/anaconda3/envs/py38/bin/pbmm2"
#~/workspace/software/smrtlink_11.1.0.166339/smrtcmds/bin/pbmm2

minimap2="     ${HOME}/software/minimap2-2.23/minimap2"


#================================================================================================
if [ $# -le 0 ]; then
    usage
    exit 1
fi

FLAG_GETSEQ=false
FLAG_REMOVEBAM=true
FLAG_ONLYMAKINGINDEX=false
INPUT_FILE_EXT=""

while getopts t:T:XYRb:i:h1:2:I OPT
do
    case $OPT in
        T)  THREADS=$OPTARG
            ;;
        X)  FLAG_GETSEQ=true
            ;;
        Y)  FLAG_GETNONMAPPEDSEQ=true
            ;;
        R)  FLAG_REMOVEBAM=false
            ;;
        t)  MAPPING_TOOL=$OPTARG
            ;;    
        I)  FLAG_ONLYMAKINGINDEX=true
            ;;   
        b)  BASE_FILE=$OPTARG
            FILENAME=${BASE_FILE##*/}
            BASE_FILE_BASENAME=${FILENAME%.*}
            BASE_FILE_EXT=${FILENAME##*.}
            BASE_FILE_DIR=$(dirname ${BASE_FILE})
            ;;   
        i)  INPUT_FILE=$OPTARG
            INPUT_FILE_NAME=${INPUT_FILE##*/}
            INPUT_FILE_BASENAME=`echo ${INPUT_FILE_NAME}| sed -e "s/.gz//g"| rev | cut -f2- -d "."| rev`

            INPUT_FILE_EXT=${INPUT_FILE_NAME##*.}
            INPUT_FILE_BASE=${INPUT_FILE_BASENAME}
            assembly_type="single"
            ;;
        1)  INPUT_FILE_P1=$OPTARG
            InputFileName_P1=${INPUT_FILE_P1##*/}
            InputFilePrefix_P1=`echo ${InputFileName_P1}| sed -e "s/.gz//g"| rev | cut -f2- -d "."| rev`

            INPUT_FILE_EXT=${InputFileName_P1##*.}
            INPUT_FILE_BASE=${InputFilePrefix_P1}  #remain "_R1"
            assembly_type="double"
            ;;
        2)  INPUT_FILE_P2=$OPTARG
            INPUT_FILE_NAME_P2=${INPUT_FILE_P2##*/}
            assembly_type="double"
            ;;
        h)  usage_exit
            exit 0
            ;;
        \?) usage_exit
            ;;
    esac
done


if [ -z ${BASE_FILE} ] || [ ! -e ${BASE_FILE} ]; then
    echo "No base file provided."
    usage;     exit 0
fi

if [ -z ${INPUT_FILE} ] && ! "${FLAG_ONLYMAKINGINDEX}"  && [ -z ${INPUT_FILE_P1} ]; then
    echo "No input file provided."
    usage;     exit 0
fi

if [ ! -e ${INPUT_FILE} ] && ! "${FLAG_ONLYMAKINGINDEX}"  && [ ! -e ${INPUT_FILE_P1} ]; then
    echo "No input file exist."
    usage;     exit 0
fi

if [ -z ${MAPPING_TOOL} ] ; then
    echo "No mapping tool setting."
    usage;   exit 0
fi

output_dir="../mapping/${BASE_FILE_BASENAME}"
if [ ! -e ${output_dir} ]; then
    echo "@@@Make new directory: ${output_dir}"
    mkdir ${output_dir} -p
fi

INDEX_NAME="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}.ref"
INDEX_NAME_mmi_pac="  ${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}.pac.mmi"
INDEX_NAME_mmi_ont="  ${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}.ont.mmi"
#INDEX_NAME_ref_run="  ${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}.ref.running"
INDEX_NAME_ref_check="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}.ref.check"
output_sam="          ${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}.sam"
output_bam="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}.bam"
output_bam_check="    ${output_bam}.check"
output_bam_sort="     ${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_sorted.bam"
output_bam_calmd="    ${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_calmd.bam"
output_depth="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_depth.tsv"
output_summary="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_summary.tsv"
output_idxstats="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_idxstats.tsv"
output_counts="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_counts.tsv"
output_fastq="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}.fastq"
output_mappedReadSeqID="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_mapped.tsv"
output_unmapped_fastq="${output_dir}/${MAPPING_TOOL}_${BASE_FILE_BASENAME}_${INPUT_FILE_BASE}_unmapped.fastq"
output_sort_check=" ${output_bam_sort}.check"

trap 'exit 0' 1 2 3 15 

echo "================================================================================================"
echo "Query file:   ${INPUT_FILE_BASE}"
echo "Base file:"   ${BASE_FILE}
echo "Mapping tool: ${MAPPING_TOOL}"
echo "Output file:  ${output_bam_sort}"  
echo "================================================================================================"




#=====================================================================================================================================
if [ ${MAPPING_TOOL} = "bowtie2" ] || [ ${MAPPING_TOOL} = "bowtie2strict" ]  ; then

    echo ${INDEX_NAME}
    if [ ! -e ${INDEX_NAME_ref_check} ]; then
        echo "@@@Make bowtie2 ref: ${INDEX_NAME}" 
        #touch ${INDEX_NAME_ref_run}

        ${bowtie2_build} -f ${BASE_FILE} ${INDEX_NAME} --quiet --threads ${THREADS}

        #check
        if [ -e ${INDEX_NAME}.1.bt2 ] || [ -e ${INDEX_NAME}.1.bt2l ]; then
            touch ${INDEX_NAME_ref_check}
        fi
    else
        echo "@@@skip bowtie2-build: ${INDEX_NAME} is already exist."
    fi

    if  "${FLAG_ONLYMAKINGINDEX}" ;  then 
        echo "Only index mode. exit"
        exit 0
    fi

    read_format=""
    echo "@@@mapped query file format: ${INPUT_FILE_EXT}"
    if [ "${INPUT_FILE_EXT}" = "fa" ] || [ "${INPUT_FILE_EXT}" = "fasta" ]; then
        read_format="-f"
    else
        read_format="-q"
    fi

    if [ ! -e ${output_bam_check} ]; then
        echo "@@@mode: ${assembly_type}"
        if [ ${assembly_type} = "single"  ]; then
            INPUTFILE_OPTION="-U ${INPUT_FILE}"
        elif [ ${assembly_type} = "double"  ]; then
            INPUTFILE_OPTION="-1 ${INPUT_FILE_P1} -2 ${INPUT_FILE_P2}"
        fi

        if [ ${MAPPING_TOOL} = "bowtie2" ]; then
            command="${bowtie2} -N 1 -x ${INDEX_NAME} ${INPUTFILE_OPTION} --threads ${THREADS} ${read_format} --quiet | ${samtools} view -bS -@ ${THREADS} -o ${output_bam}"
        elif [ ${MAPPING_TOOL} = "bowtie2strict" ]; then
            command="${bowtie2} -N 0 -x ${INDEX_NAME} ${INPUTFILE_OPTION} --threads ${THREADS} ${read_format} --quiet | ${samtools} view -bS -@ ${THREADS} -o ${output_bam}"
        fi

        echo "${command}"
        eval "${command}"

    else
        echo "@@@Skipped mapping: ${output_bam} is already exist."
    fi


#--------------------------------------------------------
elif [ "${MAPPING_TOOL}" = "bwa" ]; then
    if [ ! -e ${INDEX_NAME}.bwt ]; then
        echo "@@@Make bwa ref: ${INDEX_NAME}"
        ${bwa} index -p ${INDEX_NAME} ${BASE_FILE} 
    else
        echo "@@@skip bwa index: ${INDEX_NAME} is already exist."
    fi

    if [ ! -e ${output_bam_check} ]; then
        #${bwa} aln   index_name SRR1163656.fastq > SRR1163656.sam
        #${bwa} bwasw index_name SRR1163656.fastq > SRR1163656.sam
        #${bwa} mem   index_name SRR1163656.fastq > SRR1163656.sam
        echo "!!!!!!nanopore setting!!!!!!!!!!"
        ${bwa} mem   ${INDEX_NAME} ${INPUT_FILE} -t ${THREADS} -x ont2d | ${samtools} view -bS > ${output_bam}  #nanopore setting

        #touch ${output_bam_check}

    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi   


#--------------------------------------------------------
elif [ "${MAPPING_TOOL}" = "bwa_matepair" ]; then
    if [ ! -e ${INDEX_NAME}.bwt ]; then
        echo "@@@Make bwa ref: ${INDEX_NAME}"
        ${bwa} index -p ${INDEX_NAME} ${BASE_FILE} 
    else
        echo "@@@skip bwa index: ${INDEX_NAME} is already exist."
    fi

    if [ ! -e ${output_bam_check} ]; then

        ${bwa} mem   ${INDEX_NAME} ${INPUT_FILE} -t ${THREADS} -w 0 -O 99 | ${samtools} view -bS > ${output_bam} 
        
        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi   


#--------------------------------------------------------
elif [ ${MAPPING_TOOL} = "pbmm2" ] || [ ${MAPPING_TOOL} = "pbmm2CCS" ] || [ ${MAPPING_TOOL} = "pbmm2HiFi" ]; then    
    if [ ${INPUT_FILE_EXT} != "bam" ] ; then
        echo "@@@ Error: please set bam file as input"
        usage
        exit 0
    fi

    if [ ! -e ${output_bam_check} ]; then

        source ${HOME}/anaconda3/etc/profile.d/conda.sh
        conda activate py38

        if [ ! -e ${INDEX_NAME_ref_check} ]; then
            echo "@@@Make pbmm2 ref: ${INDEX_NAME_mmi_pac}"
            command="pbmm2 index ${BASE_FILE} ${INDEX_NAME_mmi_pac} -j ${THREADS}"
            echo ${command}
            ${command}
            touch ${INDEX_NAME_ref_check}
        else
            echo "@@@skip pbmm2 index: ${INDEX_NAME_ref_check} is already exist."
        fi
        
        if  "${FLAG_ONLYMAKINGINDEX}" ;  then 
            echo "Only index mode. exit"
            exit 0
        fi

        if   [ ${MAPPING_TOOL} = "pbmm2" ]; then
            option="--preset SUBREAD"
        elif [ ${MAPPING_TOOL} = "pbmm2CCS" ]; then
            option="--preset CCS"
        elif [ ${MAPPING_TOOL} = "pbmm2HiFi" ]; then
            option="--preset HiFi"
        fi
        THREADS_half=`expr ${THREADS} / 2`

        command="pbmm2 align ${INDEX_NAME_mmi_pac} ${INPUT_FILE} ${output_bam} -j ${THREADS_half} --sort  -J ${THREADS_half} ${option}"
        echo ${command}
        ${command}
        
        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi  


#--------------------------------------------------------
elif [ ${MAPPING_TOOL} = "blasr" ]; then    
    if [ ${INPUT_FILE_EXT} != "bam" ] ; then
        echo "@@@error: please set bam file as input"
        usage
        exit
    fi

    if [ ! -e ${output_bam_check} ]; then   
        source ${HOME}/anaconda3/etc/profile.d/conda.sh
        conda activate py38

        ${blasr} ${INPUT_FILE} ${BASE_FILE} --nproc ${THREADS} --out ${output_bam} --bam
        
        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi    


#--------------------------------------------------------
elif [ ${MAPPING_TOOL} = "pbalign" ]; then    
    if [ ${INPUT_FILE_EXT} != "bam" ] ; then
        echo "@@@error: please set bam file as input"
        usage
        exit
    fi

    if [ ! -e ${output_bam_check} ]; then
        source ${HOME}/anaconda3/etc/profile.d/conda.sh
        conda activate py38

        if [ ! -e ${INPUT_FILE}.pbi ]; then
            pbindex ${INPUT_FILE}
        fi

        pbalign  ${INPUT_FILE} ${BASE_FILE} ${output_bam} --nproc ${THREADS} --concordant #--algorithmOptions="-minMatch 12 -bestn 10 -minPctIdentity 70.0"
        
        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi 


#--------------------------------------------------------
elif [ ${MAPPING_TOOL} = "pbalignCCS" ]; then
    if [ ${INPUT_FILE_EXT} != "bam" ] ; then
        echo "@@@error: please set bam file as input"
        usage
        exit
    fi

    if [ ! -e ${output_bam_check} ]; then
        echo "@@@Mapping ${MAPPING_TOOL}"
        
        source ${HOME}/anaconda3/etc/profile.d/conda.sh
        conda activate py38

        if [ ! -e ${INPUT_FILE}.pbi ]; then
            pbindex ${INPUT_FILE}
        fi

        pbalign ${INPUT_FILE} ${BASE_FILE} ${output_bam} --useccs useccsdenovo --nproc ${THREADS}
                
        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi   


#--------------------------------------------------------
elif [ ${MAPPING_TOOL} = "minimapONT" ]; then

    # create an index
    if [ ! -e ${INDEX_NAME_ref_check} ]; then
        echo "@@@Making minimap MMI file: ${mmi}"
        ${minimap2} -x map-ont -d ${INDEX_NAME_mmi_ont} ${BASE_FILE} -t ${THREADS} -I32g
        touch ${INDEX_NAME_ref_check}
    else
        echo "@@@Skip making minimap MMI file: ${INDEX_NAME_mmi_ont}"
    fi

    if [ ! -e ${output_bam_check} ]; then

        #Aligning reads to the reference genome
        #${minimap2} -ax map-ont ${mmi} ${INPUT_FILE} -t ${THREADS} -o ${output_sam}
        #${samtools} view  -@ ${THREADS} -Sb ${output_sam}     > ${output_bam}
        ${minimap2} -ax map-ont ${INDEX_NAME_mmi_ont} ${INPUT_FILE} -t ${THREADS}  -I32g| ${samtools} sort -@ ${THREADS} -O BAM - > ${output_bam}

        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi  


#--------------------------------------------------------
elif [ ${MAPPING_TOOL} = "minimapHiFi" ]; then

    # create an index
    if [ ! -e ${INDEX_NAME_ref_check} ]; then
        echo "@@@Making minimap MMI file: ${mmi}"
        ${minimap2} -x map-ont -d ${INDEX_NAME_mmi_ont} ${BASE_FILE} -t ${THREADS} -I32g
        touch ${INDEX_NAME_ref_check}
    else
        echo "@@@Skip making minimap MMI file: ${INDEX_NAME_mmi_ont}"
    fi

    if [ ! -e ${output_bam_check} ]; then

        ${minimap2} -ax map-hifi ${INDEX_NAME_mmi_ont} ${INPUT_FILE} -t ${THREADS} -I32g | ${samtools} sort -@ ${THREADS} -O BAM - > ${output_bam}
        #touch ${output_bam_check}
    else
        echo "@@@skip Mapping: ${output_bam} is already exist."
    fi  

#--------------------------------------------------------
else
    echo "@@@Error: Undefiend mapping tool"
    usage
    exit
fi


#=====================================================================================================================================
#file check
#=====================================================================================================================================
if [ ! -e ${output_bam_check} ]; then
    if [ -e ${output_bam} ] && [ -s ${output_bam} ]; then
        echo "@@@Output: ${output_bam}"
        touch ${output_bam_check}
    else
        echo "@@@ ERROR!!! no ${output_bam}"
        exit 0
    fi
fi
    
#=====================================================================================================================================
#make accessory files
#=====================================================================================================================================
if [ ! -e ${output_sort_check} ]; then
    echo "@@@convert bam format to various formats."
    
    #${samtools} view  -@ ${THREADS} -Sb ${output_sam}     > ${output_bam}
    echo "@@@commands----------------------------------------------------------------------------------------------------------------------"
    echo "@@@${samtools} sort  -@ ${THREADS}     ${output_bam}    -o ${output_bam_sort}    "
    echo "@@@${samtools} index -@ ${THREADS}     ${output_bam_sort}"
    echo "@@@${samtools} idxstats ${output_bam_sort}               > ${output_idxstats}"
    echo "@@@cut -f1,3 ${output_idxstats}                          > ${output_counts}"
    echo "@@@------------------------------------------------------------------------------------------------------------------------------"
    
    if [ -e ${output_bam}  ]; then
        ${samtools} sort  -@ ${THREADS}     ${output_bam}    -o ${output_bam_sort}    

    fi

    ${samtools} index -@ ${THREADS}     ${output_bam_sort}
    ${samtools} idxstats ${output_bam_sort}  > ${output_idxstats}
    cut -f1,3 ${output_idxstats}             > ${output_counts}
    #${samtools} depth    ${output_bam_sort} > ${output_depth}  # commented out to reduce strage usage

    #summary
    ${samtools} flagstat -@ ${THREADS}  ${output_bam_sort} > ${output_summary}
    #bamtools stats -in  ${output_bam_sort} > ${output_summary}
    echo "@@@Output: ${output_summary}"
    cat ${output_summary}
    
    #${bedtools} bamtobed -i mapped_biwa_5m_${sample}_bowtie2.bam > mapped_biwa_5m_${sample}_bowtie2.bed
    #rm mapped_biwa_5m_${sample}_bowtie2.sam mapped_biwa_5m_${sample}_bowtie2.bam

    #get aligned reads
    if   "${FLAG_GETSEQ}"; then
        #${bedtools} bamtofastq -i ${output_bam_sort} -fq ${output_fastq}
        samtools fastq ${output_bam_sort} --threads ${THREADS} -0 ${output_fastq}  #-o option is not work... 20220811
    fi

    #get unmapped reads
    if [ ! -z "${FLAG_GETNONMAPPEDSEQ}" ] ; then
        #get mapped read SeqID
        ${bedtools} bamtobed -i ${output_bam_sort} | cut -f4 > ${output_mappedReadSeqID}

        seqkit grep -f ${output_mappedReadSeqID} -v ${INPUT_FILE} -o ${output_unmapped_fastq}
        #${bedtools} bamtofastq -i ${output_bam_sort} -fq ${output_fastq}
    fi

    #remove bam file after check file existence
    if [ -e ${output_bam_sort}  ]; then
        echo "remove ${output_bam}"
        touch ${output_sort_check}
        
        if "${FLAG_REMOVEBAM}"; then 
            rm ${output_bam}
        fi
    fi

else
    echo "@@@Skip converting process"
fi

echo "@@@mapping All done."

echo ${output_bam_sort}

exit 0

