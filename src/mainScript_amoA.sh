

#===============================================================================================================================
# Satoshi Hiraoka
# hiraokas@jamstec.go.jp
#
# Description:
# This is a main scripts used in the "amoA probe capture enrichment sequencing" study.
# 
# Citation:
# **Probe capture enrichment sequencing of amoA genes discloses diverse ammonia-oxidizing archaeal and bacterial populations**
# Satoshi Hiraoka1†*, Minoru Ijichi2†, Hirohiko Takeshima2, Yohei Kumagai2, Ching-Chia Yang2, Yoko Makabe-Kobayashi2, Hideki Fukuda2, Susumu Yoshizawa2, Wataru Iwasaki2,3, Kazuhiro Kogure2, Takuhei Shiozaki2*
# 1. Research Center for Bioscience and Nanoscience (CeBN), Japan Agency for Marine-Earth Science and Technology (JAMSTEC), 2–15 Natsushima-cho, Yokosuka, Kanagawa 237–0061, Japan
# 2. Atmosphere and Ocean Research Institute, the University of Tokyo, 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-8564, Japan
# 3. Department of Integrated Biosciences, Graduate School of Frontier Sciences, the University of Tokyo, 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-0882, Japan.
# † Contributed equally
# * Corresponding author
# Email: hiraokas@jamstec.go.jp, shiozaki@g.ecc.u-tokyo.ac.jp
#===============================================================================================================================



#======================================================================================================
# QC for Illumina short reads
# - standard QC
# - remove chimera
#======================================================================================================
#remove adapter
for f in ../data/01_flat/*_R1.fastq.gz; do ./fastq_adaptor_trim.sh ${f}; done

#remove phiX reads
###for f in ../data/12_QC/amplicon*_R1.fq.gz; do ./qsub_epyc.sh 2 ./fastq_remove_phiX.sh ${f} paired 4; done

#remove reads with low complexity 
for f in ../data/12_QC/*_R1.fq.gz;  do ./qsub_DDBJ.sh epyc 4 4 4 ./fastq_remove_lowcomplexity.sh ${f} paired 4 ; done

#merge paired-end reads
./fastq_pairend_marge.sh  dirNoComp ../data/15_nocomplex_upper100

#convert fastq -> fasta
for f in ../data/16_merged/*.fastq; do
    echo ${f}
    perl fastq2fasta.pl -a ${f}
done
mv *.fa ../data/22_fasta/

# remove chimera (Illumina amplicon) ---------------------------------------
mkdir ../data/41_ampliconIllumina_concatenate_resource
for f in ../data/15_nocomplex_upper100/ampliconIllumina_*_1.fastq.gz; do
    echo ${f}
    P1=${f}
    P2=`echo ${f} | sed -e "s/_1.fastq.gz/_2.fastq.gz/g"`
    Prefix=`basename ${f} | rev | cut -f3- -d "." | cut -f2- -d "_" |rev`
    cat  ${f} |  seqkit fq2fa > ../data/41_ampliconIllumina_concatenate_resource/${Prefix}_1.fasta
    seqkit seq -pr ${P2} |  seqkit fq2fa | seqkit mutate -i 0:NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN --quiet  > ../data/41_ampliconIllumina_concatenate_resource/${Prefix}_2RevComp.fasta
done

for f in ../data/41_ampliconIllumina_concatenate_resource/ampliconIllumina_*_1.fasta; do
    echo ${f}
    P1=${f}
    P2=`echo ${f} | sed -e "s/_1.fasta/_2RevComp.fasta/g"`
    Prefix=`basename ${f} | rev |  cut -f2- -d "_" |rev`
    seqkit concate  ${P1}  ${P2}  > ../data/42_ampliconIllumina_concatenate/${Prefix}.concatenate.fasta
done

for f in ../data/42_ampliconIllumina_concatenate/*; do
    nohup ./fastq_remove_chimera.sh vsearch ${f} &
done

mkdir ../data/43_ampliconIllumina_ChimeraDetection/
mv ../removeChimera/* ../data/43_ampliconIllumina_ChimeraDetection/

mkdir ../data/44_ampliconIllumina_removeChimera_list/
for f in ../data/43_ampliconIllumina_ChimeraDetection/*.fa; do
    echo ${f}
    prefix=`basename ${f} | rev | cut -f3- -d "." | rev`
    seqkit fx2tab -n ${f} | cut -f1 -d ";" > ../data/44_ampliconIllumina_removeChimera_list/${prefix}.tsv
done

#get non-chimeric reads
mkdir ../data/45_ampliconIllumina_removeChimera/
for f in ../data/44_ampliconIllumina_removeChimera_list/*; do
    echo ${f}
    prefix=`basename ${f} | rev | cut -f2- -d "." | rev`
    cat ../data/15_nocomplex_upper100/${prefix}_1.fastq.gz | seqkit grep -f ${f} > ../data/45_ampliconIllumina_removeChimera/${prefix}_1.fastq.gz
done

#nonpareil curve analysis-------------------------------------
for f in ../data/16_merged/*.fastq; do
    FILENAME=${f##*/}
    BASE_FILENAME=${FILENAME%.*}
    
    if [ -e ../nonpareil/${FILENAME}_kmer.npo ]; then
        echo "Already exist. Skip: ${SAMPLE_NAME}"
        continue
    fi

    echo ${SAMPLE_NAME}
    ./qsub_short.sh 6 ./nonpareil.sh ${f} 6
done


#======================================================================================================
# QC for PacBio HiFi reads 
# - remove reads with <400bp and >1000 bp length 
# - remove primer region (20bp)
# - remove chimera
#======================================================================================================
cat ../data/30_PacBio_amplicon/amplicon_DNA.fastq      | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_DNA.fastq
cat ../data/30_PacBio_amplicon/amplicon_OT3-0m.fastq   | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_OT3-0m.fastq
cat ../data/30_PacBio_amplicon/amplicon_OT3-50m.fastq  | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_OT3-50m.fastq
cat ../data/30_PacBio_amplicon/amplicon_OT3-B-5m.fastq | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_OT3-B-5m.fastq

# remove chimera (pacbio amplicon) 
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_DNA.fastq      &
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_OT3-0m.fastq   &
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_OT3-50m.fastq  &
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_OT3-B-5m.fastq &

cp ../removeChimera/amplicon_*.fa ../data/33_removeChimera/

#size check
conda activate py38
for f in ../data/33_removeChimera/*; do
    ./fasta_seqlen_averageWithSD.sh ${f}
done


#======================================================================================================
# assemby using Illumina short reads
#======================================================================================================
#For Metagenomic data
for f in ../data/15_nocomplex_upper100/*DNA_1.fastq.gz; do 
#for f in ../data/15_nocomplex_upper100/2-7-20*_1.fastq.gz; do 
    filename=${f##*/}             
    filenamewe=`echo ${filename} | rev | cut -f2- -d "_"`  
 
    fasta2=`echo ${f} | sed -e "s/_1.fastq/_2.fastq/g"`  
    ./qsub_DDBJ.sh epyc 8 20 40 ./assembly_metagenomic.sh -1 ${f} -2 ${fasta2} -t SPAdes -T 8 -M 150;   # may take longer time
done

#For transcriptomic data
#for f in ../data/15_nocomplex_upper100/[0-2]*_1.fastq.gz; do 
for f in ../data/15_nocomplex_upper100/[b-z]*_1.fastq.gz; do 
    filename=${f##*/}             
    filenamewe=`echo ${filename} | rev | cut -f2- -d "_"`  


    fasta2=`echo ${f} | sed -e "s/_1.fastq/_2.fastq/g"`
    ./qsub_DDBJ.sh epyc 8 12 40 ./assembly_metagenomic.sh -1 ${f} -2 ${fasta2}  -t SPAdes -r -T 8 -M 110; 
done

#file move
cd ../assembly/
for f in metaSPAdes_*/; do 
    ln -s ${f}/contigs.fasta ${f%/}.fasta 
done
for f in rnaSPAdes_*/; do 
    ln -s ${f}/transcripts.fasta ${f%/}.fasta 
done
cd ../src/


#======================================================================================================
# Get CDSs from assembled contigs
#======================================================================================================
for f in ../assembly/*SPAdes*.fasta; do 
    ./qsub_short.sh 1 ./genecall_prodigal.sh ${f} meta; 
done
mkdir ../OTU/0_AllGenes/ -p
for f in ../gene/genecall/*SPAdes_*.fna; do
    filename=${f##*/}             # SPAdes_0-DNA.faa
    filenamewe=${filename%.*}     # SPAdes_0-DNA
    prefix=`echo ${filenamewe} | cut -f2- -d "_"`  #0-DNA

    mv ${f} ../OTU/0_AllGenes/${prefix}.fasta
    rm ../gene/genecall/${filenamewe}.faa
    rm ../gene/genecall/${filenamewe}.gff
done


#======================================================================================================
# similarity search of CDSs against CuMMO gene sequence database
#======================================================================================================
db=../DB/amoA_ijichiReference.dmnd
dbname=${db##*/}
dbnamewe=${dbname%.*}
thread=4

for f in ../OTU/0_AllGenes/*.fasta; do
    filename=${f##*/}
    filenamewe=${filename%.*}
    outputFilename=../blast/${dbnamewe}_${filenamewe}.blast  

    command="diamond blastx --db ${db} --query ${f} --outfmt 6 --max-target-seqs 1 --evalue 1e-30 --threads ${thread} --subject-cover 50 --query-cover 50 --id 60 > ${outputFilename}"
    echo  ${command}
    eval "${command}"
done


#======================================================================================================
# Get potential amoA gene sequences
#======================================================================================================
for f in ../blast/amoA_ijichiReference_*.blast; do 
    echo ${f}
    filename=${f##*/}             # amoA_ijichiReference_0-St2.cds.blast
    prefix=`echo ${filename} | cut -f3- -d "_" | rev | cut -f3- -d "." | rev`  #noncapture_ON1-0m
    ./getseq_blast_output.sh ../OTU/0_AllGenes/${prefix}.fasta  ${f}
done

#merge potential amoA gene sequences in each dataset 
#length filtering
cat ../CandidateSeq/amoA_ijichiReference_[0-2]*DNA.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_MockDNA.fasta
cat ../CandidateSeq/amoA_ijichiReference_rare*DNA.fasta      | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_RareMockDNA.fasta
cat ../CandidateSeq/amoA_ijichiReference_*St-2.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_testSt2.fasta
cat ../CandidateSeq/amoA_ijichiReference_*capture*.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_Real.fasta

#rename seqID for clustering
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_MockDNA.fasta     | awk '/^>/{print ">DNA_"  ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_MockDNA.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_RareMockDNA.fasta | awk '/^>/{print ">Rare_" ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_RareMockDNA.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_testSt2.fasta     | awk '/^>/{print ">St2_"  ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testSt2.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_Real.fasta        | awk '/^>/{print ">Real_" ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta


#======================================================================================================
# remove chimeric sequences
#======================================================================================================
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_MockDNA.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_RareMockDNA.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testSt2.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta
mv ../removeChimera/* ../OTU/3_removeChimera/


#======================================================================================================
# reads mapping to pre-clustering sequences
#======================================================================================================
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa     -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_RareMockDNA.fa -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_testSt2.fa     -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa        -I

for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa     -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_RareMockDNA.fa -i ${f}; done 
for f in ../data/16_merged/*St2.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_testSt2.fa     -i ${f}; done 
for f in ../data/16_merged/*capture*.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa        -i ${f}; done 

for f in ../data/33_removeChimera/ampliconPacBio_DNA.fa;   do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa -i ${f}; done
for f in ../data/33_removeChimera/ampliconPacBio_OT3-*.fa; do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa    -i ${f} ; done
for f in ../data/45_ampliconIllumina_removeChimera/ampliconIllumina_*_1.fastq.gz; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa  -1 ${f} -2 `echo ${f}| sed -e "s/_1.fastq/_2.fastq/g"`  ; done 

mv ../mapping/AmoA_ORF* ../mapping/pre-clustering/

#summarize
for f in ../mapping/pre-chimera_removal/AmoA_ORF_Real/*_summary.tsv; do
    sample=`basename ${f} | cut -f5-6 -d "_"`
    ratio=`cat ${f} | grep "primary mapped (" | sed -e "s/(//g"`
    echo -e "$sample\t$ratio"
done

#blast check 
diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../removeChimera/AmoA_ORF_Real.chimeras --outfmt 6 --max-target-seqs 1 --threads 6 -o ../blast/chimera_ORF_Real.blast


#======================================================================================================
# clustering for OTU definition
#======================================================================================================
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa     -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_RareMockDNA.fa -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_testSt2.fa     -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_Real.fa        -s MMseq2 -t 0.97 -c 0.5
mv ../clustering/AmoA_ORF_*_N_0.97-0.5.MMseq2.*  ../OTU/4_OTUraw/

# Remove singlton clusters
cat ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.tsv | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.tsv        | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.singleton

cat ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.fasta     | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_MockDNA_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.fasta | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.singleton > ../OTU/5_OTUwoSingleton/OTU_RareMockDNA_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.fasta     | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_testSt2_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.fasta        | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.singleton        > ../OTU/5_OTUwoSingleton/OTU_Real_rawSeqID.fasta

#rename OTU seq names
cat ../OTU/5_OTUwoSingleton/OTU_MockDNA_rawSeqID.fasta     | awk '/^>/{print ">OTU_mock" ++i; next}{print}' > ../OTU/6_OTU/OTU_MockDNA.fasta
cat ../OTU/5_OTUwoSingleton/OTU_RareMockDNA_rawSeqID.fasta | awk '/^>/{print ">OTU_rare" ++i; next}{print}' > ../OTU/6_OTU/OTU_RareMockDNA.fasta
cat ../OTU/5_OTUwoSingleton/OTU_testOT6_rawSeqID.fasta     | awk '/^>/{print ">OTU_ot"   ++i; next}{print}' > ../OTU/6_OTU/OTU_testOT6.fasta
cat ../OTU/5_OTUwoSingleton/OTU_testSt2_rawSeqID.fasta     | awk '/^>/{print ">OTU_st"   ++i; next}{print}' > ../OTU/6_OTU/OTU_testSt2.fasta
cat ../OTU/5_OTUwoSingleton/OTU_Real_rawSeqID.fasta        | awk '/^>/{print ">OTU_real" ++i; next}{print}' > ../OTU/6_OTU/OTU_Real.fasta


#======================================================================================================
# taxonomic annotation of OTUs
#======================================================================================================
#annotation using Ijichi reference
seqkit translate ../DB/amoA_ijichiReference.fasta | sed -e "s/\s/_/g"> ../DB/amoA_ijichiReference.faa
diamond makedb --in ../DB/amoA_ijichiReference.faa -d ../DB/amoA_ijichiReference --threads 4
diamond blastx --db ../DB/amoA_ijichiReference.dmnd    --query ../OTU/6_OTU/OTU_MockDNA.fasta     --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_MockDNA.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd    --query ../OTU/6_OTU/OTU_RareMockDNA.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_RareMockDNA.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd    --query ../OTU/6_OTU/OTU_testSt2.fasta     --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_testSt2.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd    --query ../OTU/6_OTU/OTU_Real.fasta        --outfmt 6 --max-target-seqs 1 --threads 6                                             -o ../blast/AnnotationRough_OTU_Real.blast

#annotation using Alves 
diamond makedb --in amoA_Alves_NatCommun2018.faa -d amoA_Alves_NatCommun2018 --threads 4
diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../OTU/6_OTU/OTU_Real.fasta        --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 70 --query-cover 70 --id 97 -o ../blast/AnnotationAlvesAOA_OTU_Real.blast
diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../OTU/6_OTU/OTU_Real.fasta        --outfmt 6 --max-target-seqs 1 --threads 6                                             -o ../blast/AnnotationAlvesRoughAOA_OTU_Real.blast

#make phylogenetic tree
rm ../phylogeny/OTU_* -r
./phylogenetic_tree_construction.sh -i ../OTU/6_OTU/OTU_Real.fasta        -S fasttree_amoA -T 4  #GTR + G which estimated by MEGA
seqkit stat ../OTU/6_OTU/*.fasta


#======================================================================================================
#read mapping to OTU, etc.
#======================================================================================================
#Mock read mapping to the defined plasmids
./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2constracts6.fasta -I
./mapping.sh -t bowtie2 -b ../plasmid/plasmidALL.fasta         -I
./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2.fasta            -I
./mapping.sh -t bowtie2 -b ../plasmid/pET.fasta                -I

for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2constracts6.fasta  -Y   -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2constracts6.fasta  -Y   -i ${f}; done 
for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/plasmidALL.fasta          -Y   -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/plasmidALL.fasta          -Y   -i ${f}; done 
for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pET.fasta                 -Y   -i ${f}; done 

#result summarize
for f in *_summary.tsv; do
    sample=`basename ${f} | cut -f3 -d "_"`
    ratio=`cat ${f} | grep "primary mapped (" | sed -e "s/(//g"`
    echo -e "$sample\t$ratio"
done

#Mock read search using blast
makeblastdb -in ../OTU/6_OTU/OTU_MockDNA.fasta -out ../OTU/6_OTU/OTU_MockDNA -dbtype nucl 
for f in ../data/16_merged/[0-2]*DNA.fa; do 
    prefix=`basename ${f}`
    ./qsub_short.sh 6 blastn -db ../OTU/6_OTU/OTU_MockDNA -query ${f} -outfmt 6 -max_target_seqs 1 -num_threads 6 -out ../blast/blastn_${prefix}.blast
done

#Read mapping to OTUs
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_MockDNA.fasta     -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_RareMockDNA.fasta -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testOT6.fasta     -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testSt2.fasta     -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta        -I

for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_MockDNA.fasta     -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_RareMockDNA.fasta -i ${f}; done 
for f in ../data/16_merged/*St2.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testSt2.fasta     -i ${f}; done 
for f in ../data/16_merged/*capture*.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta        -i ${f}; done 

for f in ../data/33_removeChimera/ampliconPacBio_DNA.fa;   do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_MockDNA.fasta -i ${f}; done
for f in ../data/33_removeChimera/ampliconPacBio_OT3-*.fa; do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta    -i ${f} ; done
for f in ../data/45_ampliconIllumina_removeChimera/ampliconIllumina_*_1.fastq.gz; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta  -1 ${f} -2 `echo ${f}| sed -e "s/_1.fastq/_2.fastq/g"`  ; done 

# mapped read count (will take a bit long minutes) 
# for paired-end reads, only reads that both #1 and #2 were mapped were used
for f in ../mapping/OTU*/*_sorted.bam;  do 
    echo ${f}
    samtools view ${f} -F 0xC -F 0x80 | cut -f3,7 | grep -E "=|\*" | cut -f1  | sort | uniq -c  | awk -v OFS="\t" '{print $2,$1}' | grep -v "^\*" > ${f}.counts  
done

#summarize
for f in *_summary.tsv; do
    sample=`basename ${f} | cut -f4-5 -d "_"`
    ratio=`cat ${f} | grep "primary mapped (" | sed -e "s/(//g"`
    echo -e "$sample\t$ratio"
done


#======================================================================================================
# Make RPKMS table
#======================================================================================================
#calculate_FPKMS
for f in ../mapping/*/bowtie2_OTU_*[0-2]*DNA_sorted.bam.counts;                  do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g"   `; a=`seqkit stat ../data/16_merged/${prefix}.fastq                              -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_MockDNA_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_Rare*_sorted.bam.counts;                       do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g"   `; a=`seqkit stat ../data/16_merged/${prefix}.fastq                              -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_RareMockDNA_seqlen.tsv ${f}; done
for f in ../mapping/*/bowtie2_OTU_*St2_sorted.bam.counts;                        do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g"   `; a=`seqkit stat ../data/16_merged/${prefix}.fastq                              -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_testSt2_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_*capture*_sorted.bam.counts;                   do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g"   `; a=`seqkit stat ../data/16_merged/${prefix}.fastq                              -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv        ${f}; done

for f in ../mapping/*/bowtie2_OTU_MockDNA_ampliconPacBio_DNA_sorted.bam.counts;  do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g"   `; a=`seqkit stat ../data/33_removeChimera/${prefix}.fa                          -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_MockDNA_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_Real_ampliconPacBio_OT3-*_sorted.bam.counts;   do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g"   `; a=`seqkit stat ../data/33_removeChimera/${prefix}.fa                          -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv        ${f}; done
for f in ../mapping/*/bowtie2_OTU_Real_ampliconIllumina_*_sorted.bam.counts;     do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_1_sorted.bam.counts//g" `; a=`seqkit stat ../data/45_ampliconIllumina_removeChimera/${prefix}_1.fastq.gz -T | tail -1 | cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv        ${f}; done

rename "bowtie2_"    ""   ../abundance/*
rename "_sorted.bam" ""   ../abundance/*
rename "_1."         "."  ../abundance/OTU_Real_ampliconIllumina_*

#make tables
mkdir ../RPKMS
cat ../abundance/OTU_MockDNA_*.AlignedReadRatio     | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_MockDNA_//g"               > ../RPKMS/AlignedReadRatio_MockDNA.tsv
cat ../abundance/OTU_RareMockDNA_*.AlignedReadRatio | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_RareMockDNA_//g"           > ../RPKMS/AlignedReadRatio_RareMockDNA.tsv
cat ../abundance/OTU_testSt2_*.AlignedReadRatio     | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_St2_//g"                   > ../RPKMS/AlignedReadRatio_testSt2.tsv
cat ../abundance/OTU_Real_*.AlignedReadRatio        | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_Real_//g" -e "s/_1\t/\t/g" > ../RPKMS/AlignedReadRatio_Real.tsv

cat ../abundance/OTU_MockDNA_*.FPKMseq     | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_MockDNA_//g"                > ../RPKMS/RPKMS_MockDNA.tsv
cat ../abundance/OTU_RareMockDNA_*.FPKMseq | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_RareMockDNA_//g"            > ../RPKMS/RPKMS_RareMockDNA.tsv
cat ../abundance/OTU_testSt2_*.FPKMseq     | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_St2_//g"                    > ../RPKMS/RPKMS_testSt2.tsv
cat ../abundance/OTU_Real_*.FPKMseq        | sed -e "s/OTUcount_//g" | sort -r | uniq | sed -e "s/_sorted.bam//g" -e "s/bowtie2_OTU_Real_//g"  -e "s/_1\t/\t/g" > ../RPKMS/RPKMS_Real.tsv

#summarize
rm     ../RPKMS/countSummary_OTU_Real.tsv
touch  ../RPKMS/countSummary_OTU_Real.tsv
echo -e "sample\tOTU\tcount" > ../RPKMS/countSummary_OTU_Real.tsv
for f in ../mapping/OTU_Real/bowtie2_OTU_Real_*_sorted.bam.counts; do 
    echo ${f}    #OTU_rare1       9127
    prefix=`basename ${f} | sed -e "s/bowtie2_OTU_Real_//g" -e "s/_sorted.bam.counts//g" `
    cat ${f} | sed -e "s/^/${prefix}\t/g" >> ../RPKMS/countSummary_OTU_Real.tsv
done

#comvert to table format
conda activate py38
python3 tidy2table.py ../RPKMS/countSummary_OTU_Real.tsv
