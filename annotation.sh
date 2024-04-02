#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to predict gene
        Usage: sh genepredict.sh genome.fa fastq1 fastq2 thread protein.fasta
    ";  
fi	
	echo "rna-seq annotation";
	WORK=`pwd`;
	mkdir -p 06.genepredict/rna-seq/his-str-trans/index;
	cd 06.genepredict/rna-seq/his-str-trans/;
	mkdir -p bam;
	#第一个方法基于HISAT2 + StringTie+TransDecoder,推荐屏蔽重复序列的基因组
	hisat2-build $1 index/genome -p $4;
	hisat2 --dta -p $4 -x index/genome -1 $2 -2 $3 | samtools sort -@ $4 > bam/rna-seq.bam ;
	#多个转录组数据则：
	samtools merge -@ $4 bam/merged.bam bam/*.bam;
	stringtie -p $4 -o bam/merged.gtf bam/rna-seq.bam;
	gtf_genome_to_cdna_fasta.pl bam/merged.gtf $1 > transcripts.fasta;
	gtf_to_alignment_gff3.pl bam/merged.gtf > transcripts.gff3;
	TransDecoder.LongOrfs -t transcripts.fasta;
	diamond blastp --threads $4 -d /home/data/ssy28/workplace/db/uniprot_sprot/uniprot_sprot.fasta -q transcripts.fasta.transdecoder_dir/longest_orfs.pep --evalue 1e-5 --max-target-seqs 20 --out blast.xml --outfmt 5 --sensitive --id 20 --tmpdir /dev/shm --index-chunks 1;
	parsing_blast_result.pl --no-header --max-hit-num 20 --query-coverage 0.2 --subject-coverage 0.1 blast.xml > blastp.outfmt6;
	hmmscan --cpu $4 --domtblout pfam.domtbl -E 1e-5 --domE 1e-5 /home/data/ssy28/workplace/db/Pfam-A/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep;
	TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtbl --retain_blastp_hits blastp.outfmt6;
	cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3;
	echo "genewise annotation";
	cd $WORK;
	mkdir -p 06.genepredict/genewise;
	cd 06.genepredict/genewise;
	#到NCBI准备近缘物种注释好的蛋白质序列，合并到protein-all.fasta
	#perl -p -i -e 'if (m/^>/) {s/\s+.*//; s/\./_/g; }' protein-all.fasta
	homolog_genewise --cpu $4 --coverage_ratio 0.4 --evalue 1e-9 --max_gene_length 20000 $5 $1;
	homolog_genewiseGFF2GFF3 --genome $1 --min_score 15 --gene_prefix genewise genewise.gff > genewise.gff3;
