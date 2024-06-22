# Probe capture enrichment sequencing of amoA genes

The scripts were prepared for probe capture enrichment sequencing of amoA genes.

The official source code repository is at https://github.com/hiraokas/ProbeCaptureEnrichmentSequencing_amoA.

The scripts are prepared for the analysis of the study. This means that the scripts are not implemented as a stand-alone research software. 

Main script mainScript_amoA.sh and related modules are stored under "src" directory. Also, sequenecs data used in the study are atored in "data" directory

## Code
The codes are written mostly in shell script.

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
|CuMMO_Selected.fasta| Selected 20 CuMMO gene sequences used to retrieve CuMMO gene sequences from public gene database |

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
