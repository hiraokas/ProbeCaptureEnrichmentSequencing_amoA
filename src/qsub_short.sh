#!/bin/sh
#==========================================
#created:   20211111
# History: 20220324
#for DDBJ system 
##Satoshi hiraoka
#==========================================

if [ ! -e QSUB/ ]; then mkdir QSUB; fi

set -e
#command=$@
thread=${1}
command=${@:2:($#-1)}
timestamp=`date +%Y%m%d_%H-%M-%S`
RequestID=${0##*/}.${2##*/}.${1##*/}.${timestamp}.XXXXXX
tmpfile=$(mktemp "QSUB/${RequestID}")

echo ${tmpfile}

new_command=${command}

current_dir=`pwd`

cat <<EOF > ${tmpfile}
#--------<script start>--------
#!/bin/bash
#$ -pe def_slot ${thread}
#$ -l s_vmem=6G
#$ -l mem_req=6G
#$ -l d_rt=1:00:00 
#$ -l s_rt=1:00:00
#$ -l short
#$ -M hiraokas@jamstec.go.jp
#$ -V
#$ -cwd
#$ -S /bin/bash
##Satoshi hiraoka

start_time=`date +%s`
cd ${current_dir} #move to currend directory
${new_command}
./runtime.sh "${start_time}"
#--------<script end>--------
EOF

#echo ${tmpfile}
cat ${tmpfile}
chmod +x ${tmpfile}

b=`basename ${tmpfile}`

qsub_beta  -o QSUB/${b}.out -e QSUB/${b}.err -V ${tmpfile} -N ${RequestID}

exit 0
