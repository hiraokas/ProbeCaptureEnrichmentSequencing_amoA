#!/bin/bash

function usage() {
    cat <<EOF
========================================================================================================================================
Description:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    Fork:    20230118
    History: 20230924
    - For jamstec ES system 
    - 64 threads and 128 GB per 1 CPU  (= 1 ResourceSet), 2 CPU           per 1 node, in CPU node (up to 512 slot)
    - 16 threads and 512 GB per 1 unit (= 1 ResourceSet), 8 units (2 CPU) per 1 node, in GPU node (up to   8 slot)
    - Que list:  
        -  cpu
            - [automatically queue into either of the queue]
            - cpu_S (< 1 RS) : <6H  
            - cpu_L (>=2 RS) : <48H
        -  cpu_DEBUG   :       <6H
        -  gpu 
            - (< 2 RS) :       <48H
            - (>=3 RS) :       <48H
Usage:
    $(basename ${0})  Que  ResourceSet_host(1,2,3...) ResourceSet_CPU(1/2 for cpu, 1-8 for gpu)  Hour  Commands...
Exp:
    ./qsub.sh cpu 1 2 48 XX.sh XX XX ..
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
ResourceSet_host=${2}
ResourceSet_CPU=${3}
hour=${4}
command=${@:5:($#-1)}

timestamp=`date +%Y%m%d-%H-%M-%S`
a=`echo ${que##*/} | cut -f1 -d "."`;
b=`echo "${command}" | cut -f1  -d " "| sed -e "s/.\///g" |  cut -f1 -d "."  `;
RequestID=${a}.${b}.${timestamp}.XXXXXX
fileScript=$(mktemp "QSUB/${RequestID}")
fileOutput="${fileScript}.stdout"
fileError="${fileScript}.stderr"

new_command=${command}
current_dir=`pwd`

echo "RequestID: ${RequestID}"

if [ ${que} == "cpu" ] || [ ${que} == "cpu_DEBAG" ]; then 
OMP_NUM_THREADS=` echo "${ResourceSet_host} * ${ResourceSet_CPU} * 64 " | bc `
cat <<EOF > ${fileScript}
#--------<script start>--------
#!/bin/bash
#PBS -q ${que}
#PBS -T openmpi
#PBS -b ${ResourceSet_host}
#PBS --custom cpusetnum-lhost=${ResourceSet_CPU}
#PBS -l elapstim_req=${hour}:00:00 
#PBS -M hiraokas@jamstec.go.jp
#PBS -v OMP_NUM_THREADS=${OMP_NUM_THREADS}

##module load IntelMPI/2021.1.1
##module load Intel_oneAPI/2021.1.1/compiler
#module load OpenMPI/4.0.5

#source ~/.bash_profile
export LANG=en_US

cd ${current_dir}   #cd $PBS_O_WORKDIR
${new_command}   
# 2> ${fileError} 1> ${fileOutput} 
#rm ${fileScript}.err
#rm ${fileScript}.out
#--------<script end>--------
EOF

elif [ ${que} == "gpu" ] || [ ${que} == "gpu_DEBAG" ]; then  
OMP_NUM_THREADS=` echo "${ResourceSet_host} * ${ResourceSet_CPU} * 64 " | bc `

cat <<EOF > ${fileScript}
#--------<script start>--------
#!/bin/bash
#PBS -q ${que}
#PBS -T openmpi
#PBS -b ${ResourceSet_host}
#PBS --custom gpusetnum-lhost=${ResourceSet_CPU}
#PBS -l elapstim_req=${hour}:00:00 
#PBS -M hiraokas@jamstec.go.jp
#PBS -v OMP_NUM_THREADS=${OMP_NUM_THREADS}

##module load IntelMPI/2021.1.1
##module load Intel_oneAPI/2021.1.1/compiler
##module load OpenMPI/4.0.5

#source ~/.bash_profile
export LANG=en_US

cd ${current_dir}   #cd $PBS_O_WORKDIR
${new_command}  
# 2> ${fileError} 1> ${fileOutput} 
#rm ${fileScript}.err
#rm ${fileScript}.out
#--------<script end>--------
EOF

fi



#echo ${fileScript}
cat ${fileScript}
chmod +x ${fileScript}
b=`basename ${fileScript}`

qsub -o QSUB/${b}.out -e QSUB/${b}.err -V ${fileScript} -N ${RequestID}
 
exit 0

