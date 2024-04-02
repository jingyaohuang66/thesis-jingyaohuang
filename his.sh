#/usr/bin/bash
	WORK=`pwd`
	#mkdir -p star/rawdata;
	#cd star/rawdata;
	#下载转录组数据
	#cd $WORK
	echo "fastp";
	mkdir -p 01.cleandata;
	/home/data/ssy28/miniconda3/envs/assemble/bin/fastp -w 16 -i $1 -I $2 -o ./01.cleandata/illumina_1.clean.fq.gz -O ./01.cleandata/illumina_2.clean.fq.gz -z 4 -h ./01.cleandata/illumina_clean.html -j ./01.cleandata/illumina_clean.json;	
    mkdir -p 02.fastqc;
	echo "fastqc cleandata";
	/home/data/ssy28/miniconda3/envs/assemble/bin/fastqc -t 40 -o 02.fastqc ./01.cleandata/*.fq.gz $1 $2 ;
	/home/data/ssy28/miniconda3/envs/assemble/bin/multiqc 02.fastqc -o 02.fastqc;
	mkdir -p 03.star/index;
	cp ../genome.fa 03.star/index/genome.fa;
	ulimit -n 999999;
	/home/data/ssy28/miniconda3/envs/star/bin/STAR --runThreadN 40 --runMode genomeGenerate --genomeDir 03.star/index --genomeFastaFiles 03.star/index/genome.fa
	cd 03.star;
	/home/data/ssy28/miniconda3/envs/star/bin/STAR --runThreadN 40 --runMode alignReads --genomeDir index --readFilesIn ../01.cleandata/illumina_1.clean.fq.gz ../01.cleandata/illumina_2.clean.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 40 --limitBAMsortRAM 2842479818
	/home/data/ssy28/miniconda3/envs/his-string-trans/bin/stringtie Aligned.sortedByCoord.out.bam -p 40 -o sample.gtf
	/home/data/ssy28/miniconda3/envs/his-string-trans/bin/gffread -E sample.gtf -o - | sed -e "s#transcript#match#g" -e "s#exon#match_part#g" > sample.gff
	awk 'BEGIN{OFS="\t"} $7 >= 2{if($4==1){st="+"}else{st="-"} print $1,"STAR","intron",$2,$3,$7,st,".","."}' SJ.out.tab > STAR.gff
	/home/data/ssy28/software/gmes_linux_64_4/gmes_petap.pl --sequence index/genome.fa --ET STAR.gff --fungus --cores 60