#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to assemble genome 
        Usage: sh assemble.sh long.fastq fastq1 fastq2 genome-size thread
    ";  
fi

echo "flye assemble";
mkdir -p 04.assemble/flye;
flye --pacbio-raw $1 --out-dir 04.assemble/flye --genome-size $4 --threads 20 -i 3 --plasmids;
echo "polish1 genome";
mkdir -p 04.assemble/polish/index;
bwa index -p 04.assemble/polish/index/rawgenome 04.assemble/flye/assembly.fasta;
cd 04.assemble/polish/;
axel -n 100 -o ./ https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar;
chmod 700 pilon-1.24.jar;
bwa mem -t $5 index/rawgenome $2 $3 | samtools view -@ $5 -Sb - > align1.bam;
samtools sort -@ $5 -o align1.sorted.bam align1.bam;
samtools index -@ $5 align1.sorted.bam;
sambamba markdup -t $5 align1.sorted.bam align1.sorted.markup.bam;
rm align1.bam align1.sorted.bam -f;
samtools view -@ $5 -q $5 -b align1.sorted.markup.bam > align1.sorted.markup.filter.bam;
samtools index -@ $5 align1.sorted.markup.filter.bam;
java -Xms100g -Xmx100g -XX:ParallelGCThreads=$5 -jar pilon-1.24.jar --genome ../flye/assembly.fasta --threads $5 --frags align1.sorted.markup.filter.bam --fix all --changes --output pilon_polished1 --vcf &> pilon1.log;
echo "polish2 genome";
#plion2
bwa index pilon_polished1.fasta;
bwa mem -t $5 pilon_polished1.fasta $2 $3 | samtools view -@ $5 -Sb - > align2.bam;
samtools sort -@ $5 -o align2.sorted.bam align2.bam;
samtools index -@ $5 align2.sorted.bam;
sambamba markdup -t $5 align2.sorted.bam align2.sorted.markup.bam;
rm align2.bam align2.sorted.bam -f;
samtools view -@ $5 -q $5 -b align2.sorted.markup.bam > align2.sorted.markup.filter.bam;
samtools index -@ $5 align2.sorted.markup.filter.bam;
java -Xms100g -Xmx100g -XX:ParallelGCThreads=$5 -jar pilon-1.24.jar --genome ./pilon_polished1.fasta --threads $5 --frags align2.sorted.markup.filter.bam --fix all --changes --output pilon_polished2 --vcf &> pilon2.log;
echo "polish3 genome";
#plion3
bwa index pilon_polished2.fasta;
bwa mem -t $5 pilon_polished2.fasta $2 $3 | samtools view -@ $5 -Sb - > align3.bam;
samtools sort -@ $5 -o align3.sorted.bam align3.bam;
samtools index -@ $5 align3.sorted.bam;
sambamba markdup -t $5 align3.sorted.bam align3.sorted.markup.bam;
rm align3.bam align3.sorted.bam -f;
samtools view -@ $5 -q $5 -b align3.sorted.markup.bam > align3.sorted.markup.filter.bam;
samtools index -@ $5 align3.sorted.markup.filter.bam;
java -Xms100g -Xmx100g -XX:ParallelGCThreads=$5 -jar pilon-1.24.jar --genome ./pilon_polished2.fasta --threads $5 --frags align3.sorted.markup.filter.bam --fix all --changes --output pilon_polished3 --vcf &> pilon3.log;
echo "mv genome";
#由大到小按照序列排序
seqkit sort -r --by-length pilon_polished3.fasta -o genome_1.fasta;
#seqkit替换id为指定字符串(直接按照序列顺序，把序列id替换为Chr01,Chr02)
seqkit replace --ignore-case --pattern .+ -r "Scaffold{nr}" --nr-width 2 genome_1.fasta -o ../genome.fa;
echo "done";
#转大写seqkit seq -u