#!/bin/sh
#$ -m be
#$ -cwd
#$ -pe threads 1

function usage() {
    cat <<EOF
========================================================================================================================================
Description:
    Created: <20180910
    History:  20200224
    History:  20220221
    History:  20230315 (update mafft3.5)
    History:  20230315
    - Constructing phylogenetic tree from protein-coding/16S rRNA sequences stored in single multi-fasta file.
    - At least >=5 seqs were required.
Install:
    MAFFT (Also see https://mafft.cbrc.jp/alignment/software/installation_without_root.html )
        - MakeFile
            PREFIX = /usr/local    -> PREFIX = /home/hiraoka-s/workspace/software/mafft-7.490-with-extensions/ 
            BINDIR = $(PREFIX)/bin -> BINDIR = /home/hiraoka-s/workspace/software/mafft-7.490-with-extensions//bin
Usage:
    $(basename ${0}) [command] [<options>]
Required:
    -i        input multi-fasta file [.faa (amino acid) / .fasta (16S rRNA)]
    -S        tools [raxml, fasttree, fasttree_16S, fasttree_amoA]
                - fasttree      : JTT+CAT
                - fasttree_16S  : gtr, cat
                - fasttree_amoA : gamma, gtr, cat
Options:
    -F        Force running (=False)
    -M        Only perform Multiple alignment  (=False)
    -A        Only perform Tree reconstraction (=False)
    -T        Threads (=6)
    -h        print this
========================================================================================================================================
EOF
}

function usage() {
    cat <<EOF
========================================================================================================================================
INSTALL:
    wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz
    cd mafft-7.505-with-extensions/core
    mkdir ../bin
    cp Makefile Makefile_original
    cat Makefile_original | sed -e "s/\/usr\/local/\${HOME}\/software\/mafft-7.505-with-extensions/g"  > Makefile
    make clean
    make -j 4
    make install
========================================================================================================================================
EOF
}

usage_exit() {
    usage
    exit 1
}

SEQ_THRESHOLD=5
QUERY_FILE=""
THREADS=6
FORCE=0
ONLY_MULTIPLE_ALIGNMENT=0
ONLY_TREE_RECONST=0
TREE_TOOL=""

MAFFT="        ${HOME}/software/mafft-7.505-with-extensions/bin/mafft"
FastTree="     ${HOME}/software/FastTreeMP"
RAxML="        ${HOME}/software/standard-RAxML-8.2.12/raxmlHPC-PTHREADS"
ElConcatenero="${HOME}/software/ElConcatenero/ElConcatenero.py"
output_dir_alignment="../alignment/"
output_dir_prefix="../phylogeny/"
current_OMP_NUM_THREADS=${OMP_NUM_THREADS}


while getopts i:vhFMAT:S: OPT
do
    case $OPT in
        i)  FLAG_A=0
            QUERY_FILE=$OPTARG
            ;;
        h)  usage_exit
            exit 0
            ;;
        S)  TREE_TOOL=$OPTARG
            ;;
        F)  FORCE=1
            ;;
        M)  ONLY_MULTIPLE_ALIGNMENT=1
            ;;
        A)  ONLY_TREE_RECONST=1
            ;;
        T)  THREADS=$OPTARG
            ;;
        \?) usage_exit
            ;;
    esac
done

if [ -z ${QUERY_FILE} ]; then
    echo "No input file name."
    usage_exit
    exit 1
fi

#set filenames
QUERY_FILENAME="${QUERY_FILE##*/}"
QUERY_BASE_FILENAME="${QUERY_FILENAME%.*}"
output_dir="${output_dir_prefix}${QUERY_BASE_FILENAME}"
output_align_file="${output_dir}/aln_${QUERY_BASE_FILENAME}.ali"
output_align_phy="${output_dir}/aln_${QUERY_BASE_FILENAME}"  #automatically added ".phy"
output_tree_file_base="${output_dir}/${TREE_TOOL}_${QUERY_BASE_FILENAME}"
output_tree_file="${output_tree_file_base}.newick"
output_tree_file_log="${output_tree_file_base}.log"

#check input directory
if [ -z ${QUERY_FILE} ]; then
    echo "ERROR: inccorect file name: ${QUERY_FILE}"
    exit 1
fi

if [ -z ${TREE_TOOL} ]; then
    echo "ERROR: inccorect TREE_TOOL name: ${TREE_TOOL}"
    exit 1
fi

if [ ! -e ${output_dir} ]; then
    echo "make new directory: ${output_dir}"
    mkdir ${output_dir} -p
elif [ ${FORCE} = 1 ]; then
    echo "force proceeding..."
else
    echo "skiped this blast search. output dir: ${output_dir}"
fi

echo "========================================================================="
echo "MAFFT START"
echo "========================================================================="

#multiple alignment
if [ ! ${ONLY_TREE_RECONST} = 1 ]; then
    if [ ! -e ${output_align_file} ]; then
        echo "multiple alignment using MAFFT..."
        echo "  Input filename: ${QUERY_FILE}"

        #MAFFT
        command="${MAFFT} --anysymbol --auto --thread ${THREADS} --quiet ${QUERY_FILE} > ${output_align_file}"  #fasta
        echo ${command}; eval ${command}

        echo "convert fasta to phylip"
        python3 ${ElConcatenero} -if fasta -of phylip -in ${output_align_file} -o ${output_align_phy}
        
        echo "  Output: ${output_align_file}"
        echo "MAFFT END"
    else
        echo "Skiped MAFFT alignment"
    fi
else
    echo "Skiped multiple alignment."
fi
echo "========================================================================="
echo ""

if [ ONLY_MULTIPLE_ALIGNMENT = 1 ]; then
    echo "Done."
    exit 0
fi


echo "========================================================================="
echo "Phylogenetic tree" 
echo "========================================================================="

if [ ${TREE_TOOL} = "fasttree" ]; then
    echo "FastTree START"
    echo "Phylogenetic tree reconstruction using FastTree... (Jones-Taylor-Thorton + CAT model)" 
        
    #check seq count > 5
    #seq_count=`cat ${input_filename}|grep ">"| wc -l`
    #if [ ${seq_count} < ${SEQ_THRESHOLD} ];then
    #   continue
    #   echo "skiped ${input_filename}"
    #else
    #   echo "${input_filename}"
    #fi

    export OMP_NUM_THREADS=${THREADS}
    command="${FastTree} -log ${output_tree_file_log} ${output_align_file} > ${output_tree_file}"
    echo ${command}; eval ${command}

    echo "  Output: ${output_tree_file}"
    export OMP_NUM_THREADS=${current_OMP_NUM_THREADS}
    echo "FastTree END"


elif [ ${TREE_TOOL} = "fasttree_16S" ]; then
    echo "FastTree START"
    echo "Phylogenetic tree reconstruction using FastTree... (16S rRNA mode; nucleotide alignment with the GTR+CAT model, which selected in original paper)"  

    #fasttree
    export OMP_NUM_THREADS=${THREADS}
    command="${FastTree} -log ${output_tree_file_log} -gtr -nt ${output_align_file} > ${output_tree_file}"
    echo ${command}; eval ${command}
    
    echo "  Output: ${output_tree_file}"
    export OMP_NUM_THREADS=${current_OMP_NUM_THREADS}
    echo "FastTree END"

elif [ ${TREE_TOOL} = "fasttree_amoA" ]; then
    echo "FastTree START"
    echo "Phylogenetic tree reconstruction using FastTree... (16S rRNA mode; nucleotide alignment with the GTR+CAT model with gamma option)"  

    #fasttree
    export OMP_NUM_THREADS=${THREADS}
    command="${FastTree} -log ${output_tree_file_log} -gamma -gtr -nt ${output_align_file} > ${output_tree_file}"
    echo ${command}; eval ${command}
    
    echo "  Output: ${output_tree_file}"
    export OMP_NUM_THREADS=${current_OMP_NUM_THREADS}
    echo "FastTree END"


elif [ ${TREE_TOOL} = "raxml" ]; then
    echo "RAxML START"

    #raxml
    command="${RAxML} -p 1989 -m PROTGAMMAAUTO -s ${output_align_phy}.phy -w ${output_dir} -n ${TREE_TOOL}_${QUERY_BASE_FILENAME} -T ${THREADS}"
    echo ${command}; eval ${command}
    
    echo "  Output: ${output_tree_file}"  
    echo "RAxML END"
else
    echo "Undefined tool: ${TREE_TOOL}"
fi
echo "========================================================================="
echo ""

echo "All done"
exit 0

