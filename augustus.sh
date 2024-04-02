#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to augustus 
        Usage: sh augustus.sh fastq1 fastq2 genome thread
    ";  
fi
#注：fastq不可压缩，否则Trinity报错；
echo "Trinity De novo";
WORK=`pwd`
mkdir -p 06.genepredict/denovo/augustus/pasa;
Trinity --seqType fq --max_memory $5 --left $1  --right $2 --CPU $4 --output 06.genepredict/denovo/augustus/Trinity_Denovo &> 06.genepredict/denovo/augustus/trinity_denovo.logs;
#统计Trinity 组装结果:
TrinityStats.pl 06.genepredict/denovo/augustus/Trinity_Denovo/Trinity.fasta > 06.genepredict/denovo/augustus/Trinity_Denovo/Trinity.fasta.stats;
#提取最长的 Unigene:
extract_longest_isoforms_from_TrinityFasta.pl 06.genepredict/denovo/augustus/Trinity_Denovo/Trinity.fasta > 06.genepredict/denovo/augustus/Trinity_Denovo/unigene.longest.fasta;
echo "Trinity genome ";


echo "PASA piline";
cd 06.genepredict/denovo/augustus/pasa;
cp /home/data/ssy28/miniconda3/envs/tri-pasa-aug/opt/pasa-2.5.2/pasa_conf/pasa.alignAssembly.Template.txt ./alignAssembly.config;
echo "perl -p -i -e 's/DATABASE=<__DATABASE__>/DATABASE=./denovo.sqlite/' alignAssembly.config" | sh;
echo "perl -p -i -e 's/MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>/MIN_PERCENT_ALIGNED=80/' alignAssembly.config" | sh;
echo "perl -p -i -e 's/MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>/MIN_AVG_PER_ID=80/' alignAssembly.config" | sh;
/home/data/ssy28/miniconda3/envs/tri-pasa-aug/opt/pasa-2.5.2/bin/seqclean $WORK/06.genepredict/denovo/augustus/Trinity_Denovo/Trinity.fasta -v /home/data/ssy28/workplace/db/UniVec/UniVec;
/home/data/ssy28/miniconda3/envs/tri-pasa-aug/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $3 -t transcript.fasta.clean -T -u transcript.fasta --ALIGNERS blat,gmap --CPU $4 &> &> pasa_$(date +%Y-%m-%d-%H-%M).log;