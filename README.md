# Probe capture enrichment sequencing of amoA genes

The ammonia monooxygenase subunit A (amoA) gene has been used to investigate the phylogenetic diversity, spatial distribution, and activity of ammonia-oxidizing archaeal (AOA) and bacterial (AOB), which contribute significantly to the nitrogen cycle in various ecosystems. Amplicon sequencing of amoA is a widely used method; however, it produces inaccurate results owing to the lack of a ‘universal’ primer set. Moreover, currently available primer sets suffer from amplification biases, which can lead to severe misinterpretation. Although shotgun metagenomic and metatranscriptomic analyses are alternative approaches without amplification bias, the low abundance of target genes in heterogeneous environmental DNA restricts a comprehensive analysis to a realizable sequencing depth. In this study, we developed a probe set and bioinformatics workflow for amoA enrichment sequencing using a hybridization capture technique.

Please note that the scripts were prepared for the analysis of the study. This means that the scripts were not implemented as a stand-alone research software. The analysis was performed on [the supercomputing system at National Institute of Genetics (NIG), Research Organization of Information and Systems (ROIS)](https://sc.ddbj.nig.ac.jp/en/), and [the Earth Simulator systems at JAMSTEC](https://www.jamstec.go.jp/es/en/).

Main script mainScript_amoA.sh and related modules are stored under "src" directory. Also, amoA/CuMMO gene sequenecs used in the study are stored in "data" directory.
The official source code repository is at https://github.com/hiraokas/ProbeCaptureEnrichmentSequencing_amoA.

## Code
The codes are written in shell script and python.

| File                    | Description |
----|---- 
| clustering_seq.sh       | Sequence clustering |
| fasta_seqlen_averageWithSD.sh | Analysis of sequence length distribution |
| fastq_pairend_marge.sh  | Merge paired-end short-read sequences to single reads |
| fastq_remove_chimera.sh | Remove chimeric sequences |
| getseq_blast_output.sh  | Get sequences using blast output file |
| mainScript_amoA.sh      | Main script in this study |
| mapping.sh              | Read mapping |
| phylogenetic_tree_construction.sh | Phylogenetic tree estimation |
| qsub_DDBJ.sh            | Autogenerate script for grid engine installed in supercomputing system at DDBJ, Japan |
| qsub_ES.sh              | Autogenerate script for grid engine installed in Earth simulator system at JAMSTEC, Japan |
| qsub_short.sh           | Wrapper of qsub_DDBJ.sh |
| rename_fastafile.sh     | Convert filenames according to given correspondence table |
| tidy2table.py           | Convert tidy data to table format |

## Dataset
| File                    | Description |
----|---- 
|amoA_MockPlasmid.fasta| amoA gene sequences used for the Mock sample |
|CuMMO_DB.fasta| CuMMO gene sequence database used for probe design|
|CuMMO_Selected.fasta| Selected 20 CuMMO gene sequences used as queries to retrieve CuMMO gene sequences from public gene databases (NCBI nt and env_nt) |

## Dependencies
- fastq2fasta.pl [(direct link)](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwj2jeSvru6GAxW6slYBHah0DWUQFnoECBAQAQ&url=http%3A%2F%2Fbrianknaus.com%2Fsoftware%2Fsrtoolbox%2Ffastq2fasta.pl&usg=AOvVaw0jCezWC5YM0DBNzazLZyxs&opi=89978449)
- [TrimGalore](https://github.com/FelixKrueger/TrimGalore) - adapter trim
- [PRINSEQ++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) - remove low complexity sequences
- [FLASH](https://github.com/dstreett/FLASH2) - merge paired-end reads
- [metaSPAdes](https://github.com/ablab/spades) - Metagenomic assembly
- [rnaSPAdes](https://github.com/ablab/spades) - Transcriptomic assembly
- [Prodigal](https://github.com/hyattpd/Prodigal) - CDS prediction
- [DIAMOND](https://github.com/bbuchfink/diamond) - Similarity search
- [VSEARCH](https://github.com/torognes/vsearch) - Chimeric read prediction and removal
- [SeqKit](https://bioinf.shenwei.me/seqkit/) - Sequence manipulation including length filtering
- [MMseq2](https://github.com/soedinglab/MMseqs2) - Sequence clustering
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) - Sequence alignment
- [FastTree2](https://www.microbesonline.org/fasttree/) - Phylogenetic tree prediction
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - Read mapping
- [Nonpareil3](https://github.com/lmrodriguezr/nonpareil) - Metagenomic coverage estimation

Also we used some tools and databases for detailed data analysis in this study.
- [Ocean Data View](https://odv.awi.de/) - Oceanic diversity analysis
- [MEGA X](https://www.megasoftware.net/) - Phylogenetic tree analysis
- AOA amoA sequence database defined by Alves et al. [(Ref.)](https://www.nature.com/articles/s41467-018-03861-1) - Taxonomic assignments of AOA OTUs 



## Citation 

Hiraoka S. (2024) **Probe capture enrichment sequencing of amoA genes discloses diverse ammonia-oxidizing archaeal and bacterial populations**. *bioRxiv*. doi:[10.1101/2023.04.10.536224](https://www.biorxiv.org/content/10.1101/2023.04.10.536224v2)

```
**Probe capture enrichment sequencing of amoA genes discloses diverse ammonia-oxidizing archaeal and bacterial populations**

Satoshi Hiraoka1†*, Minoru Ijichi2†, Hirohiko Takeshima2, Yohei Kumagai2, Ching-Chia Yang2, Yoko Makabe-Kobayashi2, Hideki Fukuda2, Susumu Yoshizawa2, Wataru Iwasaki2,3, Kazuhiro Kogure2, Takuhei Shiozaki2*

1. Research Center for Bioscience and Nanoscience (CeBN), Japan Agency for Marine-Earth Science and Technology (JAMSTEC), 2–15 Natsushima-cho, Yokosuka, Kanagawa 237–0061, Japan
2. Atmosphere and Ocean Research Institute, the University of Tokyo, 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-8564, Japan
3. Department of Integrated Biosciences, Graduate School of Frontier Sciences, the University of Tokyo, 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-0882, Japan.

† Contributed equally
* Corresponding author

Email: hiraokas@jamstec.go.jp, shiozaki@g.ecc.u-tokyo.ac.jp
```
