#!/bin/sh
#$ -M hiraoka@cb.k.u-tokyo.ac.jp
#$ -m be
#$ -cwd
#$ -pe threads 6

function usage() {
    cat <<EOF
##============================================================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
	History: 20181103
    HIstory: 20230829
    - This is a script for removing chemic sequence from fasta file. 
    - You can chose two tools: vsearch and usearch2. Vsearch is recommended.
    - vsearch:
        - Threads = 1. This is a restriction by current vsearch version. 
        - See http://manpages.ubuntu.com/manpages/bionic/man1/vsearch.1.html
Usage:
    this.sh vsearch  file.fastA/Q
    this.sh usearch2 file.fastA/Q
##============================================================================
EOF
    return 0
}

if [ $# -ne 2 ]; then
    usage
    exit 1
fi

vsearch="${HOME}/software/vsearch-2.21.1-linux-x86_64/bin/vsearch"
usearch="${HOME}/software/usearch11.0.667_i86linux32"
filepath=${1}
outputdir="../removeChimera/"

threads=6
tool=${1}
fasta=${2}
str_filename=${fasta##*/}
str_ext=${str_filename##*.}
bn=`basename ${str_filename}`
prefix=` echo ${bn} | rev | cut -f2- -d "." | rev`

if [ ! -e ${outputdir} ]; then
    mkdir ${outputdir}
fi


#for f in ${filepath}/*.fastq; do \
#	echo ${f}; 
#	string_filename=${f##*/}; 
#	string_filename_we=${string_filename%.*}; 
#	out_chime_fasta=../Data/21_chimera_filtered/chimeric_${string_filename_we}.fna; 
#	out_nonchime_fasta=../Data/21_chimera_filtered/nonchimeric_${string_filename_we}.fna; 
#	out_nonchime_list=../Data/21_chimera_filtered/nonchimeric_${string_filename_we}.list; 
#	../module/usearch8.1.1861_i86linux32 -uchime_ref ${f} -db ../rdp_gold.fa -chimeras ${out_chime_fasta} -nonchimeras ${out_nonchime_fasta} -strand plus; 
#	cat ${out_nonchime_fasta} |grep ">" |cut -c 2- > ${out_nonchime_list};
#	../module/seqtk-master/seqtk subseq ${filepath}/${string_filename_we}.fastq ${out_nonchime_list} > ../Data/21_chimera_filtered/filtered_${string_filename_we}.fastq;  
#	wc ../Data/21_chimera_filtered/filtered_${string_filename_we}.fastq
#	echo "output: ../Data/21_chimera_filtered/filtered_${string_filename_we}.fastq"
#done

echo "Tool: ${tool}"
if [ ${tool} == "vsearch" ]; then
    # command details: http://manpages.ubuntu.com/manpages/bionic/man1/vsearch.1.html
    # formula: (1) in https://academic.oup.com/bioinformatics/article/27/16/2194/255262
    ${vsearch} --uchime_denovo  ${fasta} --nonchimeras ${outputdir}/${prefix}.fa --chimeras ${outputdir}/${prefix}.chimeras --sizeout --abskew 1 --threads ${threads} --fasta_score --borderline ${outputdir}/${prefix}.borderline --uchimealns ${outputdir}/${prefix}.uchimealns --dn 0.5 #--minh 0.2
elif [ ${tool} == "vsearch_uchime2" ]; then
    ${vsearch} --uchime2_denovo ${fasta} --nonchimeras ${outputdir}/${prefix}.fa --chimeras ${outputdir}/${prefix}.chimeras --sizeout --abskew 1 --threads ${threads} --fasta_score --borderline ${outputdir}/${prefix}.borderline --uchimealns ${outputdir}/${prefix}.uchimealns 
elif [ ${tool} == "usearch" ]; then
    echo "Usearch algolithm: https://www.drive5.com/usearch/manual/chimeras.html"

     ${usearch} -uchime3_denovo ${fasta} -uchimeout ${outputdir}/${prefix}.out -nonchimeras ${outputdir}/${prefix}.${str_ext}

else
    echo "Illigal tool set. Exit."
fi


echo "Output: ${outputdir}/${bn}"
echo "Done."
