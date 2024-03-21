

#==================================
# Miscellaneous
#==================================

#move work directory
cd /home/hiraoka-s/workspace/_shiozaki_AmoA/src

#convert amoA fasta to faa
seqkit translate amoA_ijichiSelection.fas > amoA_ijichiSelection.faa

#make database for search
diamond makedb --in amoA_ijichiSelection.faa -d amoA_ijichiSelection --threads 4

#rename
for i in `seq 0 9`; do rename "m${i}_" "m_" *; done
for i in `seq 0 9`; do rename "m${i}." "m." *; done

#rename for amplicon Illumina data
./rename_fastafile.sh ../rename_amplicon_Illumina.tsv ../data/40_Illumina_amplicon/ ../data/01_flat/

#change specific sample names
rename OT-6 OT6 *
rename St-2 St2 *

rename amplicon_ ampliconIllumina_ *



#==================================
# Ijichi CuMMO database abalysis
#==================================


./qsub_DDBJ.sh medium 10 6 20 ./phylogenetic_tree_construction.sh -i ../DB/amoA_ijichiReference.fasta        -S fasttree_amoA -T 20  #GTR + G which estimated by MEGA


diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../DB/amoA_ijichiReference.fasta  --outfmt 6 --max-target-seqs 1 --evalue 1e-30 --threads 16 --subject-cover 80 --query-cover 80 --id 80 -o ../blast/ijicjiReference_AlvesNatCommun2018.blast

makeblastdb -in ../DB/amoA_ijichiReference.fasta -out ../DB/amoA_ijichiReference -dbtype nucl 
blastn -db ../DB/amoA_ijichiReference -query ../DB/amoA_plasmids.fas  -outfmt 6 -max_target_seqs 1000 -evalue 1e-30 -num_threads 16  -perc_identity 95 -out ../blast/plasmids_ijicjiReference_95.blast
blastn -db ../DB/amoA_ijichiReference -query ../DB/amoA_plasmids.fas  -outfmt 6 -max_target_seqs 1000 -evalue 1e-30 -num_threads 16  -perc_identity 97 -out ../blast/plasmids_ijicjiReference_97.blast
blastn -db ../DB/amoA_ijichiReference -query ../DB/amoA_plasmids.fas  -outfmt 6 -max_target_seqs 1000 -evalue 1e-30 -num_threads 16  -perc_identity 99 -out ../blast/plasmids_ijicjiReference_99.blast


blastn -db ../DB/amoA_ijichiReference -query ../OTU/6_OTU/OTU_testOT6.fasta   -outfmt 6 -max_target_seqs 1000 -evalue 1e-30 -num_threads 16  -perc_identity 95 -out ../blast/testOT6_ijicjiReference_95.blast
blastn -db ../DB/amoA_ijichiReference -query ../OTU/6_OTU/OTU_testSt2.fasta   -outfmt 6 -max_target_seqs 1000 -evalue 1e-30 -num_threads 16  -perc_identity 95 -out ../blast/testSt2_ijicjiReference_95.blast



#==================================
# read QC for Illumina
#==================================

#remove adapter
for f in ../data/01_flat/*_R1.fastq.gz; do ./fastq_adaptor_trim.sh ${f}; done

#remove phiX reads
###for f in ../data/12_QC/amplicon*_R1.fq.gz; do ./qsub_epyc.sh 2 ./fastq_remove_phiX.sh ${f} paired 4; done

#remove reads with low complexity 
#for f in ../data/14_nophiX/_R1*.fastq.gz;  do ./fastq_remove_lowcomplexity.sh ${f} paired ; done
#for f in ../data/12_QC/*_R1.fq.gz;  do ./fastq_remove_lowcomplexity.sh ${f} paired 16 ; done
for f in ../data/12_QC/*_R1.fq.gz;  do ./qsub_DDBJ.sh epyc 4 4 4 ./fastq_remove_lowcomplexity.sh ${f} paired 4 ; done

#merge paired-end reads
./fastq_pairend_marge.sh  dirNoComp ../data/15_nocomplex_upper100

#convert fastq -> fasta
for f in ../data/16_merged/*.fastq; do
    echo ${f}
    perl fastq2fasta.pl -a ${f}
done
mv *.fa ../data/22_fasta/

#nonpareil-------------------------------------
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


#==================================
# assemby for Illumina
#==================================

#For Metagenomic
for f in ../data/15_nocomplex_upper100/*DNA_1.fastq.gz; do 
#for f in ../data/15_nocomplex_upper100/2-7-20*_1.fastq.gz; do 
    filename=${f##*/}             
    filenamewe=`echo ${filename} | rev | cut -f2- -d "_"`  

    #if [ -e ../assembly/SPAdes_${filenamewe}/transcripts.fasta ]; then 
        #echo "Skipped: ${filenamewe}"; 
        #continue; 
    #fi   
    fasta2=`echo ${f} | sed -e "s/_1.fastq/_2.fastq/g"`  
    #./qsub_DDBJ.sh medium 8 20 40 ./assembly_metagenomic.sh -1 ${f} -2 ${fasta2} -t SPAdes -T 8 -M 150;   # may take longer time
    ./qsub_DDBJ.sh epyc 8 20 40 ./assembly_metagenomic.sh -1 ${f} -2 ${fasta2} -t SPAdes -T 8 -M 150;   # may take longer time
done

#For transcriptome
#for f in ../data/15_nocomplex_upper100/[0-2]*_1.fastq.gz; do 
for f in ../data/15_nocomplex_upper100/[b-z]*_1.fastq.gz; do 
    filename=${f##*/}             
    filenamewe=`echo ${filename} | rev | cut -f2- -d "_"`  

    #if [ -e ../assembly/SPAdes_${filenamewe}/transcripts.fasta ]; then 
        #echo "Skipped: ${filenamewe}"; 
        #continue; 
    #fi   
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


#==================================
# Get CDSs (using Prodigal for all samples)
#==================================

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


#==================================
# blast ORFs against CuMMO (including amoA) sequence DBs
#==================================

#db=../DB/amoA_ijichiSelection.dmnd
#db=../DB/amoA_AlvesNatCommun2018.dmnd
db=../DB/amoA_ijichiReference.dmnd
dbname=${db##*/}
dbnamewe=${dbname%.*}
thread=4

#for f in ../OTU/0_AllGenes/rnaSPAdes/*.fna; do
for f in ../OTU/0_AllGenes/*.fasta; do
#for f in ../gene/genecall_reduce/*.fna; do
    filename=${f##*/}
    filenamewe=${filename%.*}
    outputFilename=../blast/${dbnamewe}_${filenamewe}.blast  

    #command="blastn -db ${db} -query ${f} -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads ${thread}  > ${outputFilename}"  #-num_alignments 1  -qcov_hsp_perc 70
    command="diamond blastx --db ${db} --query ${f} --outfmt 6 --max-target-seqs 1 --evalue 1e-30 --threads ${thread} --subject-cover 50 --query-cover 50 --id 60 > ${outputFilename}"

    echo  ${command}
    eval "${command}"
done


#==================================
# Get seqs using blast output file
#==================================
#for f in ../blast/amoA_fangeneSeed_meta*DNA.blast; do 
#for f in ../blast/amoA_ijichiSelection_*.blast; do 
#for f in ../blast/amoA_AlvesNatCommun2018_*.blast; do 
for f in ../blast/amoA_ijichiReference_*.blast; do 
    echo ${f}
    filename=${f##*/}             # amoA_ijichiReference_0-St2.cds.blast
    prefix=`echo ${filename} | cut -f3- -d "_" | rev | cut -f3- -d "." | rev`  #noncapture_ON1-0m
    ./getseq_blast_output.sh ../OTU/0_AllGenes/${prefix}.fasta  ${f}
done


#==================================
# merge sequences in each set after length filtering
#==================================
cat ../CandidateSeq/amoA_ijichiReference_[0-2]*DNA.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_MockDNA.fasta
cat ../CandidateSeq/amoA_ijichiReference_rare*DNA.fasta      | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_RareMockDNA.fasta
cat ../CandidateSeq/amoA_ijichiReference_*OT-6.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_testOT6.fasta
cat ../CandidateSeq/amoA_ijichiReference_*St-2.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_testSt2.fasta
cat ../CandidateSeq/amoA_ijichiReference_*capture*.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_Real.fasta

#cat ../CandidateSeq/amoA_ijichiSelection_[0-2]*DNA.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_DNA.fasta
#cat ../CandidateSeq/amoA_ijichiSelection_rare*DNA.fasta      | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_Rare.fasta
#cat ../CandidateSeq/amoA_ijichiSelection_*OT-6.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_OT-6.fasta
#cat ../CandidateSeq/amoA_ijichiSelection_*St-2.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_St-2.fasta
#cat ../CandidateSeq/amoA_ijichiSelection_*capture*.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_ORF_Real.fasta

# cat ../CandidateSeq/amoA_AlvesNatCommun2018_[0-2]*DNA.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_Alves_Mock.fasta
# cat ../CandidateSeq/amoA_AlvesNatCommun2018_rare*DNA.fasta      | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_Alves_Dilution.fasta
# cat ../CandidateSeq/amoA_AlvesNatCommun2018_*OT-6.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_Alves_OT6PCR.fasta
# cat ../CandidateSeq/amoA_AlvesNatCommun2018_*St-2.fasta         | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_Alves_St2PCR.fasta
# cat ../CandidateSeq/amoA_AlvesNatCommun2018_*capture*.fasta     | seqkit seq -m 400 -M 1000 > ../OTU/1_ORF_lenTrim/AmoA_Alves_Real.fasta


#rename seqID for clustering
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_MockDNA.fasta     | awk '/^>/{print ">DNA_"  ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_MockDNA.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_RareMockDNA.fasta | awk '/^>/{print ">Rare_" ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_RareMockDNA.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_testOT6.fasta     | awk '/^>/{print ">OT6_"  ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testOT6.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_testSt2.fasta     | awk '/^>/{print ">St2_"  ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testSt2.fasta
cat ../OTU/1_ORF_lenTrim/AmoA_ORF_Real.fasta        | awk '/^>/{print ">Real_" ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta

# cat ../OTU/1_ORF_lenTrim/AmoA_Alves_Mock.fasta     | awk '/^>/{print ">Mock_"     ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_Alves_Mock.fasta
# cat ../OTU/1_ORF_lenTrim/AmoA_Alves_Dilution.fasta | awk '/^>/{print ">Dilution_" ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_Alves_Dilution.fasta
# cat ../OTU/1_ORF_lenTrim/AmoA_Alves_OT6PCR.fasta   | awk '/^>/{print ">OT6PCR_"   ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_Alves_OT6PCR.fasta
# cat ../OTU/1_ORF_lenTrim/AmoA_Alves_St2PCR.fasta   | awk '/^>/{print ">St2PCR_"   ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_Alves_St2PCR.fasta
# cat ../OTU/1_ORF_lenTrim/AmoA_Alves_Real.fasta     | awk '/^>/{print ">Real_"     ++i; next}{print}' > ../OTU/2_ORF_RenameSeqID/AmoA_Alves_Real.fasta



#==================================
# remove chimera
#==================================
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_MockDNA.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_RareMockDNA.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testOT6.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testSt2.fasta
./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta

# ./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_Alves_Mock.fasta
# ./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_Alves_Dilution.fasta
# ./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_Alves_OT6PCR.fasta
# ./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_Alves_St2PCR.fasta
# ./fastq_remove_chimera.sh vsearch ../OTU/2_ORF_RenameSeqID/AmoA_Alves_Real.fasta

mv ../removeChimera/* ../OTU/3_removeChimera/

#mapping reads onth pre-chimera-removal sequences
./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_MockDNA.fasta -I
./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_RareMockDNA.fasta -I
./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testOT6.fasta -I
./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testSt2.fasta -I
./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta -I
for f in ../data/16_merged/*OT6.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testOT6.fasta     -i ${f}; done 
for f in ../data/16_merged/*St2.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_testSt2.fasta     -i ${f}; done 
for f in ../data/16_merged/*capture*.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta        -i ${f}; done 

for f in ../data/33_removeChimera/ampliconPacBio_OT3-*.fa; do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta    -i ${f} ; done
for f in ../data/45_ampliconIllumina_removeChimera/ampliconIllumina_*_1.fastq.gz; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/2_ORF_RenameSeqID/AmoA_ORF_Real.fasta  -1 ${f} -2 `echo ${f}| sed -e "s/_1.fastq/_2.fastq/g"`  ; done 

mv ../mapping/AmoA_ORF* ../mapping/pre-chimera_removal/

#get ummappted reads
for f in ../mapping/pre-chimera_removal/AmoA_ORF_Real/bowtie2_AmoA_ORF_Real_ampliconIllumina_*_sorted.bam; do
    filename=`basename ${f} | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_AmoA_ORF_Real_//g"`
    bedtools bamtobed -i ${f} | cut -f4 > ${f}_ReadSeqID.tsv
    seqkit grep -f ${f}_ReadSeqID.tsv -v ../data/45_ampliconIllumina_removeChimera/${filename}.fastq.gz -o ${f}.ummapped.fastq
done

#mapping reads to pre-clustering sequences
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_RareMockDNA.fa -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_testOT6.fa -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_testSt2.fa -I
./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa -I

for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa     -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_RareMockDNA.fa -i ${f}; done 
for f in ../data/16_merged/*OT6.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_testOT6.fa     -i ${f}; done 
for f in ../data/16_merged/*St2.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_testSt2.fa     -i ${f}; done 
for f in ../data/16_merged/*capture*.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa        -i ${f}; done 

for f in ../data/33_removeChimera/ampliconPacBio_DNA.fa;   do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa -i ${f}; done
for f in ../data/33_removeChimera/ampliconPacBio_OT3-*.fa; do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa    -i ${f} ; done
for f in ../data/45_ampliconIllumina_removeChimera/ampliconIllumina_*_1.fastq.gz; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/3_removeChimera/AmoA_ORF_Real.fa  -1 ${f} -2 `echo ${f}| sed -e "s/_1.fastq/_2.fastq/g"`  ; done 

mv ../mapping/AmoA_ORF* ../mapping/pre-clustering/


#for f in ../mapping/pre-chimera_removal/AmoA_ORF_Real/*_summary.tsv; do
for f in *_summary.tsv; do
    sample=`basename ${f} | cut -f5-6 -d "_"`
    #ratio=`cat ${f} | grep "primary mapped (" | cut -f6 -d " " | cut -c2-`
    ratio=`cat ${f} | grep "primary mapped (" | sed -e "s/(//g"`
    echo -e "$sample\t$ratio"
done


#blast check 
diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../removeChimera/AmoA_ORF_Real.chimeras --outfmt 6 --max-target-seqs 1 --threads 6                                             -o ../blast/chimera_ORF_Real.blast




#==================================
#clustering
#==================================
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_MockDNA.fa     -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_RareMockDNA.fa -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_testOT6.fa     -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_testSt2.fa     -s MMseq2 -t 0.97 -c 0.5
./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_ORF_Real.fa        -s MMseq2 -t 0.97 -c 0.5
mv ../clustering/AmoA_ORF_*_N_0.97-0.5.MMseq2.*  ../OTU/4_OTUraw/

# ./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_Alves_Mock.fa     -s MMseq2 -t 0.97 -c 0.5
# ./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_Alves_Dilution.fa -s MMseq2 -t 0.97 -c 0.5
# ./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_Alves_OT6PCR.fa   -s MMseq2 -t 0.97 -c 0.5
# ./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_Alves_St2PCR.fa   -s MMseq2 -t 0.97 -c 0.5
# ./clustering_seq.sh -i ../OTU/3_removeChimera/AmoA_Alves_Real.fa     -s MMseq2 -t 0.97 -c 0.5
#mv ../clustering/AmoA_Alves_*_N_0.97-0.5.MMseq2.*  ../OTU/4_OTUraw/



# Remove singlton clusters
cat ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.tsv | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_testOT6_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_testOT6_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.singleton
cat ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.tsv        | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.singleton

cat ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.fasta     | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_MockDNA_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_MockDNA_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.fasta | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_RareMockDNA_N_0.97-0.5.MMseq2.singleton > ../OTU/5_OTUwoSingleton/OTU_RareMockDNA_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_testOT6_N_0.97-0.5.MMseq2.fasta     | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_testOT6_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_testOT6_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.fasta     | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_testSt2_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_testSt2_rawSeqID.fasta
cat ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.fasta        | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.singleton        > ../OTU/5_OTUwoSingleton/OTU_Real_rawSeqID.fasta

# cat ../OTU/4_OTUraw/AmoA_Alves_Mock_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_Alves_Mock_N_0.97-0.5.MMseq2.singleton
# cat ../OTU/4_OTUraw/AmoA_Alves_Dilution_N_0.97-0.5.MMseq2.tsv | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_Alves_Dilution_N_0.97-0.5.MMseq2.singleton
# cat ../OTU/4_OTUraw/AmoA_Alves_OT6PCR_N_0.97-0.5.MMseq2.tsv   | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_Alves_OT6PCR_N_0.97-0.5.MMseq2.singleton
# cat ../OTU/4_OTUraw/AmoA_Alves_St2PCR_N_0.97-0.5.MMseq2.tsv   | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_Alves_St2PCR_N_0.97-0.5.MMseq2.singleton
# cat ../OTU/4_OTUraw/AmoA_Alves_Real_N_0.97-0.5.MMseq2.tsv     | cut -f1 | sort | uniq -c | sort -nr | sed -e "s/^ *//g" | grep "^1 " |cut -f2 -d " " > ../OTU/4_OTUraw/AmoA_Alves_Real_N_0.97-0.5.MMseq2.singleton

# cat ../OTU/4_OTUraw/AmoA_Alves_Mock_N_0.97-0.5.MMseq2.fasta    | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_Alves_Mock_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_Alves_Mock_rawSeqID.fasta
# cat ../OTU/4_OTUraw/AmoA_Alves_Dilution_N_0.97-0.5.MMseq2.fasta| seqkit grep -v -f ../OTU/4_OTUraw/AmoA_Alves_Dilution_N_0.97-0.5.MMseq2.singleton > ../OTU/5_OTUwoSingleton/OTU_Alves_Dilution_rawSeqID.fasta
# cat ../OTU/4_OTUraw/AmoA_Alves_OT6PCR_N_0.97-0.5.MMseq2.fasta  | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_Alves_OT6PCR_N_0.97-0.5.MMseq2.singleton   > ../OTU/5_OTUwoSingleton/OTU_Alves_OT6PCR_rawSeqID.fasta
# cat ../OTU/4_OTUraw/AmoA_Alves_St2PCR_N_0.97-0.5.MMseq2.fasta  | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_Alves_St2PCR_N_0.97-0.5.MMseq2.singleton   > ../OTU/5_OTUwoSingleton/OTU_Alves_St2PCR_rawSeqID.fasta
# cat ../OTU/4_OTUraw/AmoA_Alves_Real_N_0.97-0.5.MMseq2.fasta    | seqkit grep -v -f ../OTU/4_OTUraw/AmoA_Alves_Real_N_0.97-0.5.MMseq2.singleton     > ../OTU/5_OTUwoSingleton/OTU_Alves_Real_rawSeqID.fasta


#ckeck blast similarity
#blastn -db ../DB/amoA_plasmids.fas -query ../clustering/AmoA_ORF_DNA_N_0.97-0.3.MMseq2.fasta -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 4 
# blastn -db ../DB/amoA_plasmids.fas        -query ../OTU/4_OTUraw/AmoA_ORF_DNA_N_0.97-0.5.MMseq2.fasta -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 4 
# seqkit stat ../OTU/4_OTUraw/AmoA_ORF_DNA_N_0.97-0.5.MMseq2.fasta
# blastn -db ../DB/amoA_plasmids.fas        -query ../OTU/4_OTUraw/AmoA_ORF_Rare_N_0.97-0.5.MMseq2.fasta -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 4 
# blastn -db ../DB/amoA_ijichiSelection.fas -query ../OTU/4_OTUraw/AmoA_ORF_Real_N_0.97-0.5.MMseq2.fasta -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 4 


# #check phylogenetic tree
# cat ../OTU/4_OTUraw/AmoA_ORF_DNA_N_0.97-0.5.MMseq2.fasta ../DB/amoA_plasmids.fas | cut -f1 -d ";" > ../OTU/6_OTUwithReference/DNAwithPlasmid_rawOTU.fasta
# rm -r ../phylogeny/DNAwithPlasmid_rawOTU/
# ./phylogenetic_tree_construction.sh -i ../OTU/6_OTUwithReference/DNAwithPlasmid_rawOTU.fasta -S fasttree


#rename OTU seq names
#AmoA_ORF_DNA_N_0.97-0.5.MMseq2.fasta (Mock): manuary changed, file: OTU_Mock.fasta

cat ../OTU/5_OTUwoSingleton/OTU_MockDNA_rawSeqID.fasta     | awk '/^>/{print ">OTU_mock" ++i; next}{print}' > ../OTU/6_OTU/OTU_MockDNA.fasta
cat ../OTU/5_OTUwoSingleton/OTU_RareMockDNA_rawSeqID.fasta | awk '/^>/{print ">OTU_rare" ++i; next}{print}' > ../OTU/6_OTU/OTU_RareMockDNA.fasta
cat ../OTU/5_OTUwoSingleton/OTU_testOT6_rawSeqID.fasta     | awk '/^>/{print ">OTU_ot"   ++i; next}{print}' > ../OTU/6_OTU/OTU_testOT6.fasta
cat ../OTU/5_OTUwoSingleton/OTU_testSt2_rawSeqID.fasta     | awk '/^>/{print ">OTU_st"   ++i; next}{print}' > ../OTU/6_OTU/OTU_testSt2.fasta
cat ../OTU/5_OTUwoSingleton/OTU_Real_rawSeqID.fasta        | awk '/^>/{print ">OTU_real" ++i; next}{print}' > ../OTU/6_OTU/OTU_Real.fasta

#check phylogenetic tree
# cat ../OTU/6_OTU/OTU_Mock.fasta ../DB/amoA_plasmids.fas  > ../OTU/6_OTUwithReference/DNAwithPlasmid.fasta
# rm -r ../phylogeny/DNAwithPlasmid/
# ./phylogenetic_tree_construction.sh -i ../OTU/6_OTUwithReference/DNAwithPlasmid.fasta -S fasttree



#==================================
# PacBio QC 
# - remove reads with <400bp and >1000 bp length 
# - remove primer region (20bp)
#==================================

cat ../data/30_PacBio_amplicon/amplicon_DNA.fastq      | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_DNA.fastq
cat ../data/30_PacBio_amplicon/amplicon_OT3-0m.fastq   | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_OT3-0m.fastq
cat ../data/30_PacBio_amplicon/amplicon_OT3-50m.fastq  | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_OT3-50m.fastq
cat ../data/30_PacBio_amplicon/amplicon_OT3-B-5m.fastq | seqkit seq -m 400 -M 1000 | seqkit subseq -r 21:-21 > ../data/31_QC/amplicon_OT3-B-5m.fastq


#==================================
# remove chimera (pacbio amplicon) 
#==================================
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_DNA.fastq &
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_OT3-0m.fastq &
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_OT3-50m.fastq &
nohup ./fastq_remove_chimera.sh vsearch ../data/31_QC/amplicon_OT3-B-5m.fastq &

#fastq -> fasta
#for f in ../removeChimera/*.fastq; do perl fastq2fasta.pl -a ${f}; done

cp ../removeChimera/amplicon_*.fa ../data/33_removeChimera/

#size check
conda activate py38
for f in ../data/33_removeChimera/*; do
    ./fasta_seqlen_averageWithSD.sh ${f}
done

#==================================
# remove chimera (Illumina amplicon) 
#==================================
#~/workspace/software/FLASH-1.2.11-Linux-x86_64/flash -f 650 -r 300 -m 0 -o test -d ../data amplicon_ON1-0m_1.fastq.gz amplicon_ON1-0m_2.fastq.gz  -t 10

mkdir ../data/41_ampliconIllumina_concatenate_resource
for f in ../data/15_nocomplex_upper100/ampliconIllumina_*_1.fastq.gz; do
    echo ${f}
    P1=${f}
    P2=`echo ${f} | sed -e "s/_1.fastq.gz/_2.fastq.gz/g"`
    Prefix=`basename ${f} | rev | cut -f3- -d "." | cut -f2- -d "_" |rev`
    #cp ${f} ../data/41_ampliconIllumina_concatenate_resource
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
    cat ../data/15_nocomplex_upper100/${prefix}_2.fastq.gz | seqkit grep -f ${f} > ../data/45_ampliconIllumina_removeChimera/${prefix}_2.fastq.gz
    #break
done



#==================================
# OTU annotation
#==================================
# lastal ../DB/amoA_ijichiSelection ../OTU/6_OTU/OTU_Mock.fasta  -P 4 -f TAB > ../last/OTUs_search/Annotation_OTU_Mock.last
# lastal ../DB/amoA_ijichiSelection ../OTU/6_OTU/OTU_Rare.fasta  -P 4 -f TAB > ../last/OTUs_search/Annotation_OTU_Rare.last
# lastal ../DB/amoA_ijichiSelection ../OTU/6_OTU/OTU_OT-6.fasta  -P 4 -f TAB > ../last/OTUs_search/Annotation_OTU_OT-6.last
# lastal ../DB/amoA_ijichiSelection ../OTU/6_OTU/OTU_St-2.fasta  -P 4 -f TAB > ../last/OTUs_search/Annotation_OTU_St-2.last
# lastal ../DB/amoA_ijichiSelection ../OTU/6_OTU/OTU_Real.fasta  -P 4 -f TAB > ../last/OTUs_search/Annotation_OTU_Real.last

# python3 last_getBestScoreAlignment.py ../last/OTUs_search/Annotation_OTU_Mock.last
# python3 last_getBestScoreAlignment.py ../last/OTUs_search/Annotation_OTU_Rare.last
# python3 last_getBestScoreAlignment.py ../last/OTUs_search/Annotation_OTU_OT-6.last
# python3 last_getBestScoreAlignment.py ../last/OTUs_search/Annotation_OTU_St-2.last
# python3 last_getBestScoreAlignment.py ../last/OTUs_search/Annotation_OTU_Real.last

# mv ../last/OTUannotation_Annotation_OTU* ../last/OTUs_annotation/

#phylogenetic annotation using Ijichi reference
seqkit translate ../DB/amoA_ijichiReference.fasta | sed -e "s/\s/_/g"> ../DB/amoA_ijichiReference.faa
diamond makedb --in ../DB/amoA_ijichiReference.faa -d ../DB/amoA_ijichiReference --threads 4
diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ../OTU/6_OTU/OTU_MockDNA.fasta     --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_MockDNA.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ../OTU/6_OTU/OTU_RareMockDNA.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_RareMockDNA.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ../OTU/6_OTU/OTU_testOT6.fasta     --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_testOT6.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ../OTU/6_OTU/OTU_testSt2.fasta     --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_testSt2.blast

#diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_Real.blast
diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ../OTU/6_OTU/OTU_Real.fasta        --outfmt 6 --max-target-seqs 1 --threads 6                                             -o ../blast/AnnotationRough_OTU_Real.blast


#phylogenetic annotation using Ijichi selection
#blastn -db ../DB/amoA_ijichiSelection.fas -query ../OTU/6_OTU/OTU_Mock.fasta -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 4 

# diamond blastx --db ../DB/amoA_ijichiSelection.dmnd --query ../OTU/6_OTU/OTU_Mock.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_Mock.blast
# diamond blastx --db ../DB/amoA_ijichiSelection.dmnd --query ../OTU/6_OTU/OTU_Rare.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_Rare.blast
# diamond blastx --db ../DB/amoA_ijichiSelection.dmnd --query ../OTU/6_OTU/OTU_OT-6.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_OT-6.blast
# diamond blastx --db ../DB/amoA_ijichiSelection.dmnd --query ../OTU/6_OTU/OTU_St-2.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_St-2.blast
# diamond blastx --db ../DB/amoA_ijichiSelection.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/Annotation_OTU_Real.blast
# diamond blastx --db ../DB/amoA_ijichiSelection.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 50 --threads 6 --subject-cover 80 --query-cover 80 --id 80 -o ../blast/Annotation_OTU_WeakReal.blast
# mv ../blast/Annotation_OTU_*.blast ../blast/OTUs/


#phylogenetic annotation using Alves 
#diamond makedb --in amoA_Alves_NatCommun2018.faa -d amoA_Alves_NatCommun2018 --threads 4
#diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 80 --query-cover 80 --id 97 -o ../blast/amoAphylogeny_OTU_Real.blast
diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 1 --threads 6 --subject-cover 70 --query-cover 70 --id 97 -o ../blast/AnnotationAlvesAOA_OTU_Real.blast
diamond blastx --db ../DB/amoA_AlvesNatCommun2018.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 1 --threads 6                                             -o ../blast/AnnotationAlvesRoughAOA_OTU_Real.blast

#phylogenetic annotation using nr
# diamond blastx --db ${HOME}/database/nr/nr_20220929/nr.dmnd --query ../OTU/6_OTU/OTU_Mock.fasta --outfmt 6 --max-target-seqs 40 --threads 6 --subject-cover 50 --query-cover 50 --id 60 -o ../blast/nr_OTU_Mock.blast
# diamond blastx --db ${HOME}/database/nr/nr_20220929/nr.dmnd --query ../OTU/6_OTU/OTU_Rare.fasta --outfmt 6 --max-target-seqs 40 --threads 6 --subject-cover 50 --query-cover 50 --id 60 -o ../blast/nr_OTU_Rare.blast
# diamond blastx --db ${HOME}/database/nr/nr_20220929/nr.dmnd --query ../OTU/6_OTU/OTU_OT-6.fasta --outfmt 6 --max-target-seqs 40 --threads 6 --subject-cover 50 --query-cover 50 --id 60 -o ../blast/nr_OTU_OT-6.blast
# diamond blastx --db ${HOME}/database/nr/nr_20220929/nr.dmnd --query ../OTU/6_OTU/OTU_St-2.fasta --outfmt 6 --max-target-seqs 40 --threads 6 --subject-cover 50 --query-cover 50 --id 60 -o ../blast/nr_OTU_St-2.blast
# diamond blastx --db ${HOME}/database/nr/nr_20220929/nr.dmnd --query ../OTU/6_OTU/OTU_Real.fasta --outfmt 6 --max-target-seqs 40 --threads 6 --subject-cover 50 --query-cover 50 --id 60 -o ../blast/nr_OTU_Real.blast


# estimate original taxonomy of each sequence using BASTA
# seqkit translate ../OTU/6_OTU/OTU_Real.fasta > ../OTU/6_OTU/OTU_Real.faa
# diamond blastp --query ../OTU/6_OTU/OTU_Real.faa --db ${HOME}/database/UniRef/uniref90.dmnd  --outfmt 6  --threads 12 --very-sensitive  --max-target-seqs 100  -o ../blast/OTU_Real_uniref90.blast
# sed 's/UniRef90_//' ../blast/OTU_Real_uniref90.blast > ../blast/OTU_Real_uniref90_edit.blast

# source ${HOME}/anaconda3/etc/profile.d/conda.sh
# conda activate py39
# basta sequence -d ${HOME}/database/BASTA/ ../blast/OTU_Real_uniref90_edit.blast ../BASTA/OTU_Real_uniref90_edit.basta prot --minimum 1 -p 60 --identity 90  #p:Percentage of hits


#make phylogenetic tree within the OTUs
rm ../phylogeny/OTU_* -r
# ./phylogenetic_tree_construction.sh -i ../OTU/6_OTU/OTU_MockDNA.fasta     -S fasttree -T 4
# ./phylogenetic_tree_construction.sh -i ../OTU/6_OTU/OTU_RareMockDNA.fasta -S fasttree -T 4
# ./phylogenetic_tree_construction.sh -i ../OTU/6_OTU/OTU_testOT6.fasta     -S fasttree -T 4
# ./phylogenetic_tree_construction.sh -i ../OTU/6_OTU/OTU_testSt2.fasta     -S fasttree -T 4
./phylogenetic_tree_construction.sh -i ../OTU/6_OTU/OTU_Real.fasta        -S fasttree_amoA -T 4  #GTR + G which estimated by MEGA

seqkit stat ../OTU/6_OTU/*.fasta


#==================================
# blast ampliconIllumina against amoA DB
#==================================
for f in ../data/43_ampliconIllumina_removeChimera/*.fa ; do
    echo ${f}
    prefix=`basename ${f} | rev | cut -f2- -d "." | rev`
     ./qsub_short.sh 6  diamond blastx --db ../DB/amoA_ijichiReference.dmnd --query ${f}        --outfmt 6 --max-target-seqs 1 --threads 6 -o ../blast/reads_${prefix}.blast 
done


#==================================
# setup for OTU analysis
#==================================
#make last db
lastdb ../OTU/6_OTU/OTU_MockDNA     ../OTU/6_OTU/OTU_MockDNA.fasta  
lastdb ../OTU/6_OTU/OTU_RareMockDNA ../OTU/6_OTU/OTU_RareMockDNA.fasta  
lastdb ../OTU/6_OTU/OTU_testOT6     ../OTU/6_OTU/OTU_testOT6.fasta 
lastdb ../OTU/6_OTU/OTU_testSt2     ../OTU/6_OTU/OTU_testSt2.fasta  
lastdb ../OTU/6_OTU/OTU_Real        ../OTU/6_OTU/OTU_Real.fasta  

mkdir ../last/OTUs_search/ -p 
mkdir ../last/OTUs_annotation/ -p 

#get gene length 
conda activate py38
python3 fasta_seqlen_eachSeqID.py  ../OTU/6_OTU/OTU_MockDNA.fasta
python3 fasta_seqlen_eachSeqID.py  ../OTU/6_OTU/OTU_RareMockDNA.fasta
python3 fasta_seqlen_eachSeqID.py  ../OTU/6_OTU/OTU_testOT6.fasta
python3 fasta_seqlen_eachSeqID.py  ../OTU/6_OTU/OTU_testSt2.fasta
python3 fasta_seqlen_eachSeqID.py  ../OTU/6_OTU/OTU_Real.fasta


#==================================
#read  alignment on OTU
#==================================
#using alignment approach (last)
# for f in ../data/22_fasta/[0-2]*DNA.fa;    do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh medium 4 4 10 bash -c \"lastal ../OTU/6_OTU/OTU_MockDNA     ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done
# for f in ../data/22_fasta/rare*.fa;        do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh medium 4 4 10 bash -c \"lastal ../OTU/6_OTU/OTU_RareMockDNA ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done
# for f in ../data/22_fasta/*OT-6.fa;        do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh medium 4 4 10 bash -c \"lastal ../OTU/6_OTU/OTU_testOT6     ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done
# for f in ../data/22_fasta/*St-2.fa;        do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh medium 4 4 10 bash -c \"lastal ../OTU/6_OTU/OTU_testSt2     ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done
# for f in ../data/22_fasta/*capture*.fa;    do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh medium 4 4 10 bash -c \"lastal ../OTU/6_OTU/OTU_Real        ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done

# #map pacbio amplicon to OTU
# for f in ../data/33_removeChimera/amplicon_DNA.fa;   do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh epyc 4 4 2 bash -c \"lastal ../OTU/6_OTU/OTU_MockDNA ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done
# for f in ../data/33_removeChimera/amplicon_OT3-*.fa; do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_DDBJ.sh epyc 4 4 2 bash -c \"lastal ../OTU/6_OTU/OTU_Real    ${f}  -P 4 -f TAB  \> ../last/${a}.last\" ; done

# #count
# for f in ../last/*.last; do python3 last_getBestScoreAlignment.py ${f}; done
# #for f in ../last/a*.last; do python3 last_getBestScoreAlignment.py ${f}; done

# #calculate FPKMseq
# for f in ../last/OTUcount_[0-2,a]*DNA.tsv; do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/22_fasta/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_MockDNA_seqlen.tsv     ${f}; done
# for f in ../last/OTUcount_rare*.tsv;       do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/22_fasta/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_RareMockDNA_seqlen.tsv ${f}; done
# for f in ../last/OTUcount_*OT-6.tsv;       do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/22_fasta/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_testOT6_seqlen.tsv     ${f}; done
# for f in ../last/OTUcount_*St-2.tsv;       do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/22_fasta/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_testSt2_seqlen.tsv     ${f}; done
# for f in ../last/OTUcount_*capture*.tsv;   do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/22_fasta/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv        ${f}; done
# for f in ../last/OTUcount_amplicon_DNA.tsv;      do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/33_removeChimera/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_MockDNA_seqlen.tsv ${f}; done
# for f in ../last/OTUcount_amplicon_OT3-*.tsv;    do prefix=`basename ${f} | cut -f2- -d "_" | cut -f1 -d "."`; a=`grep ">" ../data/33_removeChimera/${prefix}.fa | wc -l`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv    ${f}; done


#Mock read mapping on plasmid
./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2constracts6.fasta     -I
./mapping.sh -t bowtie2 -b ../plasmid/plasmidALL.fasta -I
./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2.fasta -I
./mapping.sh -t bowtie2 -b ../plasmid/pET.fasta -I

for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2constracts6.fasta  -Y   -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2constracts6.fasta  -Y   -i ${f}; done 
#for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2.fasta  -Y   -i ${f}; done 
#for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pTAKN-2.fasta  -Y   -i ${f}; done 
for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/plasmidALL.fasta  -Y   -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/plasmidALL.fasta  -Y   -i ${f}; done 
for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../plasmid/pET.fasta  -Y   -i ${f}; done 


for f in *_summary.tsv; do
    sample=`basename ${f} | cut -f3 -d "_"`
    #ratio=`cat ${f} | grep "primary mapped (" | cut -f6 -d " " | cut -c2-`
    ratio=`cat ${f} | grep "primary mapped (" | sed -e "s/(//g"`
    echo -e "$sample\t$ratio"
done


#Mock read search using blast
makeblastdb -in ../OTU/6_OTU/OTU_MockDNA.fasta -out ../OTU/6_OTU/OTU_MockDNA -dbtype nucl 
for f in ../data/16_merged/[0-2]*DNA.fa; do 
    prefix=`basename ${f}`
    ./qsub_short.sh 6 blastn -db ../OTU/6_OTU/OTU_MockDNA -query ${f} -outfmt 6 -max_target_seqs 1 -num_threads 6 -out ../blast/blastn_${prefix}.blast
done


#using mapping base (bowtie2)
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_MockDNA.fasta     -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_RareMockDNA.fasta -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testOT6.fasta     -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testSt2.fasta     -I
./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta        -I

for f in ../data/16_merged/[0-2]*DNA.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_MockDNA.fasta     -i ${f}; done 
for f in ../data/16_merged/rare*.fastq;     do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_RareMockDNA.fasta -i ${f}; done 
for f in ../data/16_merged/*OT6.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testOT6.fasta     -i ${f}; done 
for f in ../data/16_merged/*St2.fastq;      do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_testSt2.fasta     -i ${f}; done 
for f in ../data/16_merged/*capture*.fastq; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta        -i ${f}; done 

for f in ../data/33_removeChimera/ampliconPacBio_DNA.fa;   do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_MockDNA.fasta -i ${f}; done
for f in ../data/33_removeChimera/ampliconPacBio_OT3-*.fa; do a=`basename ${f} | cut -f1 -d "."`; echo ${a}; ./qsub_short.sh 6  ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta    -i ${f} ; done
#for f in ../data/15_nocomplex_upper100/ampliconIllumina_*_1.fastq.gz; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta  -1 ${f} -2 `echo ${f}| sed -e "s/_1.fastq/_2.fastq/g"`  ; done 
for f in ../data/45_ampliconIllumina_removeChimera/ampliconIllumina_*_1.fastq.gz; do ./qsub_short.sh 6 ./mapping.sh -t bowtie2 -b ../OTU/6_OTU/OTU_Real.fasta  -1 ${f} -2 `echo ${f}| sed -e "s/_1.fastq/_2.fastq/g"`  ; done 


#extract mapped reads that only one time mapped on unique reference reads. (filter multiple-hits)
# for f in ../mapping/*/*_sorted.bam; do
#     samtools view -F 0x100  
# done

#count (take bit minutes) (for paired-end reads, only reads #1 that both paired-ends were mapped were used)
#samtools view ../mapping/OTU_Real/bowtie2_OTU_Real_capture_OT3-0m_sorted.bam | cut -f3 |sort | uniq -c  | awk -v OFS="\t" '{print $2,$1}'
for f in ../mapping/OTU*/*_sorted.bam;  do 
    echo ${f}
    #if [ -e ${f}.counts  ]; then echo "Already exist. skip."; continue; fi
    #samtools view ${f} -F 0xC | cut -f3 |sort | uniq -c  | awk -v OFS="\t" '{print $2,$1}' | grep -v "^\*" > ${f}.counts  
    samtools view ${f} -F 0xC -F 0x80 | cut -f3,7 |grep -E "=|\*" | cut -f1  | sort | uniq -c  | awk -v OFS="\t" '{print $2,$1}' | grep -v "^\*" > ${f}.counts  
done


for f in *_summary.tsv; do
    sample=`basename ${f} | cut -f4-5 -d "_"`
    #ratio=`cat ${f} | grep "primary mapped (" | cut -f6 -d " " | cut -c2-`
    ratio=`cat ${f} | grep "primary mapped (" | sed -e "s/(//g"`
    echo -e "$sample\t$ratio"
done


#calculate_FPKMS
for f in ../mapping/*/bowtie2_OTU_*[0-2]*DNA_sorted.bam.counts;          do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/16_merged/${prefix}.fastq     -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_MockDNA_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_Rare*_sorted.bam.counts;               do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/16_merged/${prefix}.fastq     -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_RareMockDNA_seqlen.tsv ${f}; done
for f in ../mapping/*/bowtie2_OTU_*OT6_sorted.bam.counts;                do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/16_merged/${prefix}.fastq     -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_testOT6_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_*St2_sorted.bam.counts;                do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/16_merged/${prefix}.fastq     -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_testSt2_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_*capture*_sorted.bam.counts;           do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/16_merged/${prefix}.fastq     -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv        ${f}; done

for f in ../mapping/*/bowtie2_OTU_MockDNA_ampliconPacBio_DNA_sorted.bam.counts;  do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/33_removeChimera/${prefix}.fa -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_MockDNA_seqlen.tsv     ${f}; done
for f in ../mapping/*/bowtie2_OTU_Real_ampliconPacBio_OT3-*_sorted.bam.counts;   do prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_sorted.bam.counts//g" `; a=`seqkit stat ../data/33_removeChimera/${prefix}.fa -T | tail -1 |cut -f4`;  python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv    ${f}; done
for f in ../mapping/*/bowtie2_OTU_Real_ampliconIllumina_*_sorted.bam.counts; do 
    prefix=`basename ${f} | cut -f4- -d "_" | sed -e "s/_1_sorted.bam.counts//g" `; 
    echo ${prefix}
    #a=`seqkit stat ../data/15_nocomplex_upper100/${prefix}_1.fastq.gz -T | tail -1 |cut -f4`;  
    a=`seqkit stat ../data/45_ampliconIllumina_removeChimera/${prefix}_1.fastq.gz -T | tail -1 |cut -f4`;  
    python3 calculate_FPKMS.py ${a} ../OTU/6_OTU/OTU_Real_seqlen.tsv    ${f}; 
done

rename "bowtie2_"    ""   ../abundance/*
rename "_sorted.bam" ""   ../abundance/*
rename "_1."         "."  ../abundance/OTU_Real_ampliconIllumina_*


#==================================
# make table
#==================================
mkdir ../RPKMS
cat ../abundance/OTU_MockDNA_*.AlignedReadRatio     |sed -e "s/OTUcount_//g"| sort -r|uniq | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_MockDNA_//g"     > ../RPKMS/AlignedReadRatio_MockDNA.tsv
cat ../abundance/OTU_RareMockDNA_*.AlignedReadRatio |sed -e "s/OTUcount_//g"| sort -r|uniq | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_RareMockDNA_//g" > ../RPKMS/AlignedReadRatio_RareMockDNA.tsv
cat ../abundance/OTU_testOT6_*.AlignedReadRatio     |sed -e "s/OTUcount_//g"| sort -r|uniq | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_OT6_//g"         > ../RPKMS/AlignedReadRatio_testOT6.tsv
cat ../abundance/OTU_testSt2_*.AlignedReadRatio     |sed -e "s/OTUcount_//g"| sort -r|uniq | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_St2_//g"         > ../RPKMS/AlignedReadRatio_testSt2.tsv
cat ../abundance/OTU_Real_*.AlignedReadRatio        |sed -e "s/OTUcount_//g"| sort -r|uniq | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_Real_//g" | sed -e "s/_1\t/\t/g"  > ../RPKMS/AlignedReadRatio_Real.tsv


cat ../abundance/OTU_MockDNA_*.FPKMseq      |sed -e "s/OTUcount_//g"| sort -r|uniq  | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_MockDNA_//g"                        > ../RPKMS/RPKMS_MockDNA.tsv
cat ../abundance/OTU_RareMockDNA_*.FPKMseq  |sed -e "s/OTUcount_//g"| sort -r|uniq  | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_RareMockDNA_//g"                    > ../RPKMS/RPKMS_RareMockDNA.tsv
cat ../abundance/OTU_testOT6_*.FPKMseq      |sed -e "s/OTUcount_//g"| sort -r|uniq  | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_OT6_//g"                            > ../RPKMS/RPKMS_testOT6.tsv
cat ../abundance/OTU_testSt2_*.FPKMseq      |sed -e "s/OTUcount_//g"| sort -r|uniq  | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_St2_//g"                            > ../RPKMS/RPKMS_testSt2.tsv
cat ../abundance/OTU_Real_*.FPKMseq         |sed -e "s/OTUcount_//g"| sort -r|uniq  | sed -e "s/_sorted.bam//g" | sed -e "s/bowtie2_OTU_Real_//g" | sed -e "s/_1\t/\t/g"    > ../RPKMS/RPKMS_Real.tsv


#rarefaction from Real samples
rm     ../RPKMS/countSummary_OTU_Real.tsv
touch  ../RPKMS/countSummary_OTU_Real.tsv
echo -e "sample\tOTU\tcount" > ../RPKMS/countSummary_OTU_Real.tsv
for f in ../mapping/OTU_Real/bowtie2_OTU_Real_*_sorted.bam.counts; do 
    echo ${f}
    #OTU_rare1       9127
    prefix=`basename ${f} | sed -e "s/bowtie2_OTU_Real_//g" | sed -e "s/_sorted.bam.counts//g" `
    cat ${f} | sed -e "s/^/${prefix}\t/g" >> ../RPKMS/countSummary_OTU_Real.tsv
done

conda activate py38
python3 tidy2table.py ../RPKMS/countSummary_OTU_Real.tsv
