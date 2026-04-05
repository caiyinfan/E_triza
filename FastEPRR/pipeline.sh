####Genome Survey
# 数据质控
SOAPnuke filter -l 150 -q 0.5 -n 0.1 -f read1.fq.gz -r read2.fq.gz -1 clean_1.fq.gz -2 clean_2.fq.gz

# K-mer 频率统计及基因组特征评估 (K=17)
kmerfreq -k 17 -t 10 -p prefix clean_1.fq.gz clean_2.fq.gz > kmer_freq.stat
gce -g genome_size_estimate -f kmer_freq.stat -M 2000 -m 1 -b 0
 -H 1 > gce_output.txt

####De novo Assembly & Scaffolding
# HiFi 序列组装 (使用 -l3 参数消除冗余)
hifiasm -o E_triza.asm -t 32 -l3 HiFi_reads.fastq.gz

# 将组装图转换为 fasta 格式
awk '/^S/{print ">"$2;print $3}' E_triza.asm.p_ctg.gfa > E_triza.asm.p_ctg.fa

# Hi-C 数据比对与互作图谱构建
juicer.sh -d juicer_dir -D juicer_scripts -s MboI -z E_triza.asm.p_ctg.fa -p chrom.sizes -t 32

# 染色体挂载
3d-dna/run-asm-pipeline.sh -r 2 E_triza.asm.p_ctg.fa merged_nodups.txt
# 之后在 Juicebox v1.11.08 中进行人工纠错

####Assembly Quality Assessment
# BGI 短序列比对率评估
bwa index E_triza_genome.fa
bwa mem -t 32 E_triza_genome.fa clean_1.fq.gz clean_2.fq.gz | samtools sort -@ 8 -o bgi_mapped.bam
samtools depth bgi_mapped.bam > bgi_depth.txt

# PacBio HiFi 序列比对率评估
minimap2 -ax map-hifi -t 32 E_triza_genome.fa HiFi_reads.fastq.gz | samtools sort -@ 8 -o hifi_mapped.bam

# Merqury 评估组装准确性 (QV)
meryl count k=21 memory=64 threads=32 clean_*.fq.gz output meryl_db
merqury.sh meryl_db E_triza_genome.fa merqury_out

# BUSCO 核心基因完整度评估 (针对雀形目 passeriformes_odb10)
busco -i E_triza_genome.fa -l passeriformes_odb10 -o busco_genome -m geno -c 32

####Repeat Annotation
# 从头预测 LTR
LTR_FINDER -w 2 -C -t 32 E_triza_genome.fa > ltr_finder.out

# 构建物种特异性重复序列库
BuildDatabase -name Etriza_db E_triza_genome.fa
RepeatModeler -database Etriza_db -threads 32 -LTRStruct > run.out

# 同源比对 (结合自定义库和 Repbase)
RepeatMasker -pa 32 -nolow -no_is -norna -lib combined_repeat_library.fasta E_triza_genome.fa

# 串联重复序列注释 (使用指定参数)
trf E_triza_genome.fa 2 7 7 80 10 50 2000 -d -h

####Protein-coding Gene Annotation
# 1. 基于转录组数据的预测
hisat2-build E_triza_genome.fa E_triza_idx
hisat2 -p 32 -x E_triza_idx -1 rna_1.fq.gz -2 rna_2.fq.gz -S rna_mapped.sam
samtools sort -@ 8 rna_mapped.sam -o rna_mapped.bam
stringtie -p 32 -o stringtie_assembled.gtf rna_mapped.bam
# 使用 PASA 整合转录本
Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g E_triza_genome.fa -t stringtie_assembled.fasta --ALIGNERS blat,gmap

# 2. 基于同源蛋白的预测
tblastn -query related_species_proteins.fa -db E_triza_genome.db -evalue 1e-5 -num_threads 32 -out tblastn.out
exonerate --model protein2genome --query related_species_proteins.fa --target E_triza_genome.fa --showtargetgff yes > exonerate.gff

# 3. Ab initio 从头预测
augustus --species=chicken --noInFrameStop=true --strand=both E_triza_genome.fa > augustus.gff
genscan HumanIso.smat E_triza_genome.fa > genscan.out

# 4. 使用 MAKER 整合基因集
maker -cpus 32 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

####ncRNA & Functional Annotation
# tRNA 和 rRNA 注释
tRNAscan-SE -o tRNA.out -f tRNA.ss -m tRNA.stats E_triza_genome.fa
rnammer -S euk -m lsu,ssu,tsu -xml rnammer.xml -gff rnammer.gff -f rnammer.fasta E_triza_genome.fa

# miRNA 和 snRNA 注释
cmsearch --cpu 32 -Z 1000 --outfmt dfam Rfam.cm E_triza_genome.fa > rfam_results.txt

# 蛋白质功能数据库比对
diamond blastp --db swissprot.dmnd --query E_triza_proteins.fa --evalue 1e-5 --threads 32 --outfmt 6 --out swissprot_anno.txt
interproscan.sh -i E_triza_proteins.fa -f TSV -goterms -pa -cpu 32 -d interpro_output/
