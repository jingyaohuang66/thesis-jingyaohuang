#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to ncRNA genome 
        Usage: sh ncrna.sh fasta thread genomesize
    ";  
fi

echo "RepeatModeler";
WORK=`pwd`;
mkdir -p 05.genomefeature/ncrna/tRNA;
cd 05.genomefeature/ncrna;
tRNAscan-SE -o ./tRNA/tRNA.out -f ./tRNA/tRNA.ss -m ./tRNA/tRNA.stats -j ./tRNA/tRNA.gff --thread $2 $1;
convert_tRNAScanSE_to_gff3.pl -i ./tRNA/tRNA.out > ./tRNA/tRNAscan.gff3;
mkdir -p rRNA;
barrnap --kingdom euk --threads $2 --quiet --outseq rRNA/rRNA.fa $1 > rRNA/rrna.gff3;
#--kingdom参数指定物种类型，bac代表细菌，arc代表古菌，euk代表真核生物，mito代表后生动物线粒体；--threads指定并行的线程数。
#--quiet 安静模式 --outseq 输出rrna的序列
#预测结果以GFF3格式保存
#其他ncRNA注释http://t.zoukankan.com/jessepeng-p-15392809.html
genome_total=$3  #基因组大小esl-seqstat $Genome
CMnumber=`grep "ACC" /home/data/ssy28/workplace/db/ncrna/Rfam.cm | wc -l`;
Z=`echo $genome_total*2*$CMnumber/1000000 | bc`;
#echo $Z
cmscan -Z $Z --cut_ga --rfam --nohmmonly --tblout genome.tblout --fmt 2 --cpu $2 --clanin /home/data/ssy28/workplace/db/ncrna/Rfam.clanin /home/data/ssy28/workplace/db/ncrna/Rfam.cm $1 > genome.cmscan
infernal-tblout2gff.pl --cmscan --fmt2 genome.tblout >genome.infernal.ncRNA.gff3
