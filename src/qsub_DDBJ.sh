#!/bin/bash

function usage() {
    cat <<EOF
========================================================================================================================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Created: 20220204
    History: 20220324
    - Que list:  epyc  intel  medium  gpu 
Usage:
    $(basename ${0})  Que  Threads  Memory(per_thread)  Hour  Commands...
Exp:
    ./qsub.sh medium 6 8 48 XX.sh XX XX ..
========================================================================================================================================
EOF
}

usage_exit() {
    usage
    exit 1
}

if [ $# -le 3 ];   then usage_exit; fi
if [ ! -e QSUB/ ]; then mkdir QSUB; fi

set -e
que=${1}
thread=${2}
memory=${3}
hour=${4}
command=${@:5:($#-1)}

timestamp=`date +%Y%m%d-%H-%M-%S`
a=`echo ${0##*/} | cut -f1 -d "."`;
b=`echo ${5##*/} | cut -f1 -d "."`;
RequestID=${a}.${1##*/}.${b}.${timestamp}.XXXXXX
tmpfile=$(mktemp "QSUB/${RequestID}")
new_command=${command}
current_dir=`pwd`

option=""
if [ ${que} == "gpu" ]; then
    option="-l cuda=1"
fi

cat <<EOF > ${tmpfile}
#----------------<script start>----------------
#!/bin/bash
#$ -pe def_slot ${thread}
#$ ${option}
#$ -l s_vmem=${memory}G
#$ -l mem_req=${memory}G
#$ -l d_rt=${hour}:00:00 
#$ -l s_rt=${hour}:00:00
#$ -l ${que}
#$ -M hiraokas@jamstec.go.jp
#$ -V
#$ -cwd
#$ -S /bin/bash

cd ${current_dir} 
${new_command}
#----------------<script end>----------------
EOF

cat ${tmpfile}
chmod +x ${tmpfile}
b=`basename ${tmpfile}`

qsub_beta -o QSUB/${b}.out -e QSUB/${b}.err -V ${tmpfile} -N ${RequestID}

exit 0
