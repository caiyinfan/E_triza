# Chromosome-level genome assembly of the Tristram's Bunting (*Emberiza tristrami*)

## 📖 About The Project

This repository contains the software configurations, analysis pipelines, and supplementary data files used for the de novo genome assembly and annotation of the Tristram's bunting (*Emberiza tristrami*). 

The genome was assembled using a combination of **PacBio HiFi** long reads and **Hi-C** sequencing data, resulting in a high-quality, chromosome-level reference genome spanning 1.33 Gb, with a scaffold N50 of 64.94 Mb. This work has been accepted in principle at *Scientific Data*.

## 📂 Repository Contents

* `pipeline.sh`: The core command-line workflows for genome survey, assembly, quality assessment, repeat annotation, and protein-coding gene prediction.
* `software.txt` / `Software.xls`: Comprehensive lists of all bioinformatics software tools, versions, and key parameters used in the analyses.
* `pipline.png`: A visual flowchart illustrating the complete bioinformatics workflow.
* `quality_statistics.dat.txt`: Detailed statistics evaluating the quality, contiguity, and completeness of the genome assembly.
* `whfs-xs-242144_Emberiza_tristrami.order`: The Hi-C scaffolding order file used for anchoring contigs.

## ✉️ Contact
For any questions regarding the genome assembly or data, please contact:
* Baowei Zhang (zhangbw@ahu.edu.cn)
* Tingli Hu (hutingli@scsio.ac.cn)
