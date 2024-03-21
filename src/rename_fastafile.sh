#!/bin/sh

function usage() {
    cat <<'EOF'
#==========================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Created: 20210209
    History: 20230722
Usage:
    thish.sh rename_rule.tsv input_dir output_dir [ext]
Exp:
    ./rename.sh ../sample_data/rename_rule.tsv ../data/01_seq_data_flat/ ../data/03seq_data_metagenome/ fa
#==========================================
EOF
    return 0
}   

rename_file=${1}
inputdir=${2}
outputdir=${3}
ext=""

if [ $# -le 2 ]; then
    echo "option: $#"
    usage
    exit 1
fi

if [ $# -gt 3 ]; then
    ext=.${4}
fi

if [ ! -e ${outputdir} ]; then
    mkdir ${outputdir}
fi

while read line; do

    first_chara=`echo $line|cut -c1`
    if [ ${first_chara} = "#" ]; then  #comment
        continue
    fi

    pre_filename=` echo ${line}|cut -f1 -d " "`   #KR1111-Lander-05
    post_filename=`echo ${line}|cut -f2 -d " "|tr -d '\r\n'`  #ZB_S005.0-007.5
  
    if [ ${pre_filename} = " " ]; then
        continue
    fi

    if [ ${pre_filename} = "0" ]; then
        echo "${post_filename}: no original file"
        continue
    fi


    if [ ! -e  ${outputdir}/${post_filename}${ext} ]; then
        echo ${inputdir}/${pre_filename}${ext} "->" ${outputdir}/${post_filename}${ext}

        cp ${inputdir}/${pre_filename}${ext} ${outputdir}/${post_filename}${ext}

    else
        echo "@@@skip: already exist: ${post_filename}"
    fi

done < ${rename_file}

echo "@@@Done."