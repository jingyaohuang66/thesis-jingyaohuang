#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to assess genome size
        Usage: sh genomesize.sh fastq1 fastq2 kmer p thread
    ";  
fi

echo "fastp rawdata";
	mkdir -p 01.cleandata;
	fastp -w 16 -i $1 -I $2 -o ./01.cleandata/illumina_1.clean.fq.gz -O ./01.cleandata/illumina_2.clean.fq.gz -D -z 4 -q 20 -u 30 -h ./01.cleandata/illumina_clean.html -j ./01.cleandata/illumina_clean.json;	
    mkdir -p 02.fastqc;
	echo "fastqc cleandata";
	fastqc -o 02.fastqc ./01.cleandata/*.fq.gz $1 $2 ;
	multiqc 02.fastqc -o 02.fastqc;
	mkdir -p 03.assessment;
	echo "assessment genome size";
	unpigz -p $5 ./01.cleandata/illumina_1.clean.fq.gz ./01.cleandata/illumina_2.clean.fq.gz;
	jellyfish count -C -s 50G -m $3 -t $5 -o ./03.assessment/kmer_counts.jf ./01.cleandata/illumina_1.clean.fq ./01.cleandata/illumina_2.fastq;
	jellyfish histo -h 65535 -f -t $5 -o ./03.assessment/kmer.freq.stat.2colum ./03.assessment/kmer_counts.jf;
	genomescope2 -i ./03.assessment/kmer.freq.stat.2colum -o ./03.assessment/genomescope -k $3 -p $4 -m 100000;
