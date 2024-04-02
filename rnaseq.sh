#rnaseq,安装conda，使用conda装hisat2\samtools\gffread....
#sam
#bam
WORK=`pwd`
#11.ref
cd $WORK
mkdir 11.ref -p
cd 11.ref
#下载好基因组、注释信息（包括gff3、protein、fa），并传输到服务器上
mkdir -p rawdata
cd rawdata
#下载测序数据,https://sra-explorer.info/
cd ..
mkdir -p cleandata;cd cleandata
ls ../rawdata/*.fq.gz |cut -d"/" -f 3|cut -d"_" -f 1> ../sample.txt
awk '{print "fastp -w 16 -z 4 -i ../rawdata/"$1 "_1.fq.gz -I ../rawdata/"$1 "_2.fq.gz -o "$1 "_1.fq.gz -O "$1 "_2.fq.gz -h "$1 ".html -j "$1 ".json &>" $1".log"}' ../sample.txt >step1.fastp.sh
#12.index
cd $WORK
mkdir -p 12.index;cd 12.index
hisat2-build ../11.ref/genome.fa ./genome 1>hisat2-build.log 2>&1
#13.bam
cd $WORK
mkdir -p 13.bam;cd 13.bam
#hisat2
awk '{print "hisat2 --new-summary -p 8 -x ../12.index/genome -1 ../11.ref/cleandata/"$1 "_1.fq.gz -2 ../11.ref/cleandata/"$1 "_2.fq.gz -S "$1 ".sam 1>"$1 ".log 2>&1"}' ../11.ref/sample.txt > step2.hisat2.sh
#samtools
awk '{print "samtools sort -@ 20 -o "$1 ".bam "$1 ".sam"}' ../11.ref/sample.txt > step3.sort.sh
#samtools index 
#samtools index *.bam
awk '{print "samtools index -@ 20 "$1 ".bam "}' ../11.ref/sample.txt > step3.index.sh
#回传bam文件和bam.bai文件
#计算couts数
#run-featurecounts.R
gffread maker.gff3 -T -o maker.gtf
awk '{print "Rscript /home/data/ssy28/workplace/genek/转录组/workplace/rnaseq-apple-training-share/data/rnaseq-apple/1.RNASeq-ref/2.Quantification/script/run-featurecounts.R -b "$1 ".bam -g ./maker.gtf -o "$1 " &"}' ./sample.txt > step4.counts.sh
#合并
ls *.count > genes.quant_files.txt
perl /home/data/ssy28/workplace/genek/转录组/rnaseq-apple-training-share/data/rnaseq-apple/1.RNASeq-ref/3.Merge_result/script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes



