WORK = `pwd`
#11.software、准备工作
mkdir 11.software
cd 11.software
axel -n 100 -o ./ http://ftp.genek.cn:8888/xzwjhh_20211129/software/lnctools.sif
axel -n 100 -o ./ https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
axel -n 100 -o ./ http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip *.gz

#12.ref
cd $WORK
mkdir 12.ref -p
cd 12.ref
axel -n 100 -o ./ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/653/495/GCA_023653495.1_Mcan_v2/GCA_023653495.1_Mcan_v2_genomic.fna.gz
awk '{print $1}' GCA_023653495.1_Mcan_v2_genomic.fna |singularity exec ../11.software/lnctools.sif seqkit seq -u >genome.fa
axel -n 10 -o ./ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/653/495/GCA_023653495.1_Mcan_v2/GCA_023653495.1_Mcan_v2_genomic.gff.gz
mv GCA_023653495.1_Mcan_v2_genomic.gff genome.gff3
singularity exec ../11.software/lnctools.sif gffread genome.gff3 -g genome.fa -x cds.fa -y protein.fa
singularity exec ../11.software/lnctools.sif gffread genome.gff3 -T -o genome.gtf
awk '{print $0" transcript_biotype \"protein_coding\";"}' genome.gtf > genomeplas.gtf
#13.funtion
cd $WORK
mkdir 13.function
cd 13.function
git clone http://git.genek.cn:3333/zhxd2/emcp
#protein.fa提交到http://eggnog-mapper.embl.de/
axel -n 100 -o ./ http://eggnog-mapper.embl.de/MM_owju2zl0/out.emapper.annotations
singularity exec ../11.software/lnctools.sif Rscript ./emcp/emapperx.R out.emapper.annotations ../12.ref/protein.fa ./

#14.data
cd $WORK
mkdir 14.data/raw -p
mkdir 14.data/clean -p
cd 14.data/raw
#找到ncbi的数据,创建数据信息表
awk '{print "singularity exec ../../11.software/lnctools.sif kingfisher get -r "$2 " -m ena-ascp aws-http prefetch 1>"$2".log 2>&1"}' ../samples.txt >step1.download.sh
cd $WORK/14.data/clean
awk '{print "singularity exec ../../11.software/lnctools.sif fastp -w 16 -z 4 -i ../raw/"$2 "_1.fastq.gz -I ../raw/"$2 "_2.fastq.gz -o "$2 "_1.fq.gz -O "$2 "_2.fq.gz -h "$2 ".html -j "$2 ".json"}' ../samples.txt >step2.fastp.sh

#比对
cd $WORK
mkdir 21.Mapping -p
cd 21.Mapping
#构建index
singularity exec ../11.software/lnctools.sif hisat2-build -p 40 ../12.ref/genome.fa ../12.ref/genome 1>hisat2-build.log 2>&1
#比对
awk '{print "singularity exec ../11.software/lnctools.sif hisat2 --new-summary --rna-strandness RF -p 8 -x ../12.ref/genome -1 "$3" -2 "$4" -S "$2".sam 1>"$2".log 2>&1 "}' ../14.data/samples1.txt >step2.run_hisat2.sh
#sort排序
awk '{print "singularity exec ../11.software/lnctools.sif samtools sort -@ 10 -o "$2".bam "$2".sam "}' ../14.data/samples1.txt >step3.sam2bam.sh
#建立bam的index
awk '{print "singularity exec ../11.software/lnctools.sif samtools index "$2".bam "}' ../14.data/samples1.txt >step4.bamindex.sh
#统计比对率
singularity exec ../11.software/lnctools.sif Rscript step5.statistics.R

#转录本重构
#第一步，单样本转录本重构
cd $WORK
mkdir 22.Assembly -p
cd 22.Assembly
awk '{print "singularity exec ../11.software/lnctools.sif stringtie -p 8 --rf -G ../12.ref/genome.gff3 -o "$2".gtf ../21.Mapping/"$2".bam 1>"$2".log 2>&1"}' ../14.data/samples1.txt > step1.run_stringtie.sh
#第二部，合并转录本
singularity exec ../11.software/lnctools.sif stringtie --merge -o merged.gtf -G ../12.ref/genome.gff3 *.gtf



#第二步：预测编码能力
cd $WORK
mkdir 23.Identification -p
cd 23.Identification
#lncRNA鉴定
#第一步：过滤
singularity exec ../11.software/lnctools.sif FEELnc_filter.pl -i ../22.Assembly/merged.gtf -a ../12.ref/genome.gtf --monoex=-1 -s 200 -f 0 -p 40 > candidate_lncRNA.gtf 2>FEELnc_filter.log
#提取lncRNA的cds序列
singularity exec ../11.software/lnctools.sif gffread -w candidate_lncRNA.fa -g ../12.ref/genome.fa candidate_lncRNA.gtf
#提取候选 lncRNA ID
singularity exec ../11.software/lnctools.sif gtftk get_attr_value_list -i candidate_lncRNA.gtf -k transcript_id -o candidate_lncRNA.txt
## FEELnc 耗时（30min+）需要已知lncRNA序列
singularity exec ../11.software/lnctools.sif FEELnc_codpot.pl -i candidate_lncRNA.gtf -a ../12.ref/known_mRNA.gtf -l ../12.ref/known_lncRNA.gtf -g ../12.ref/genome.fa

## CPC2 根据算法
singularity exec ../11.software/lnctools.sif CPC2.py -i candidate_lncRNA.fa -o cpc2output

##PLEK 根据算法
singularity exec ../11.software/lnctools.sif python /opt/PLEK.1.2/PLEK.py -fasta candidate_lncRNA.fa -out plekoutput.txt -thread 50

## CNCI 根据算法
singularity exec ../11.software/lnctools.sif CNCI.py -f candidate_lncRNA.fa -o cnci -m pl -p 50

#### lncADeep 仅支持人和小鼠
#singularity exec ../11.software/lnctools.sif LncADeep.py -MODE lncRNA -m partial -s human -f candidate_lncRNA.fa -o LncADeep -th 50 --HMMthread 50

## Blast 对Swissprot数据库进行比对找出编码蛋白
singularity exec ../11.software/lnctools.sif diamond makedb --in ../11.software/swissprot -d ../11.software/swissprot
singularity exec ../11.software/lnctools.sif diamond blastx -p 50 -d ../11.software/swissprot -q candidate_lncRNA.fa -o diamondoutput.txt
## HMMER 60min,基于结构域进行编码蛋白的比对
singularity exec ../11.software/lnctools.sif hmmpress ../11.software/Pfam-A.hmm
singularity exec ../11.software/lnctools.sif hmmscan --tblout hmmscanoutput.tbl -o hmmscanoutput.txt --noali -E 1e-5 --cpu 50 ../11.software/Pfam-A.hmm candidate_lncRNA.fa

#取交集
singularity exec ../11.software/lnctools.sif Rscript step3.intersect.R

#提取 CDS 和 GTF
singularity exec ../11.software/lnctools.sif gtftk select_by_key -i candidate_lncRNA.gtf -k transcript_id -f lncRNA.txt -o lncRNA.gtf
singularity exec ../11.software/lnctools.sif gffread -w lncRNA.fa -g ../12.ref/genome.fa lncRNA.gtf

#合并结果
awk '{print $0" transcript_biotype \"lncRNA\";"}' lncRNA.gtf >lncRNA1.gtf
awk '{print $1}' lncRNA.fa >lncRNA1.fa
cat ../12.ref/genomeplas.gtf lncRNA1.gtf > genes.gtf

#定量
cd $WORK
mkdir 24.Quantification -p
cd 24.Quantification
#第一步：使用 featureCounts 进行表达定量
awk '{print "singularity exec ../11.software/lnctools.sif run-featurecounts.R -a transcript_id --strandSpecific 2 -b ../21.Mapping/"$2".bam -g ../23.Identification/genes.gtf -o "$2" 1>"$2".log 2>&1 "}' ../14.data/samples1.txt >step1.run_Quantification.sh
#第二步：合并成表达矩阵并进行 TMM 标准化
ls *.count >quant_files.txt
singularity exec ../11.software/lnctools.sif abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files quant_files.txt --out_prefix trans

#第三步：定量统计
singularity exec ../11.software/lnctools.sif Rscript step3.statistics.R

#第四步，按照转录本和lncRNA进行分割表达定量
singularity exec ../11.software/lnctools.sif Rscript step4.split.R


#差异表达分析
cd $WORK
mkdir 25.DE_analysis/mRNA -p
mkdir 25.DE_analysis/lncRNA -p
cd 25.DE_analysis/mRNA
#建议lncRNA 和 mRNA 分开做，分开做更准确。
#写好差异的分组contrasts.txt文件
#mRNA
singularity exec ../../11.software/lnctools.sif run_DE_analysis.pl --matrix ../../24.Quantification/mRNA.counts.matrix --method DESeq2 --samples_file ../../14.data/samples1.txt --contrasts ../contrasts.txt

#
#lncRNA
cd ../lncRNA
singularity exec ../../11.software/lnctools.sif run_DE_analysis.pl --matrix ../../24.Quantification/lncRNA.counts.matrix --method DESeq2 --samples_file ../../14.data/samples1.txt --contrasts ../contrasts.txt







































