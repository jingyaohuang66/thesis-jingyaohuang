#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to repeat genome 
        Usage: sh repeat.sh fasta genus thread
    ";  
fi

echo "RepeatModeler";
WORK=`pwd`;
mkdir -p 05.genomefeature/repeat/RepeatModeler;
cd 05.genomefeature/repeat/RepeatModeler;
BuildDatabase -name mydb $1 > genome.fasta.log;
RepeatModeler -database mydb -pa $3 -LTRStruct > run.out;
#查看物种在RepeatMaskerLib.h5中的情况
#进入RepeatMasker目录中
cd /home/data/ssy28/miniconda3/envs/repeat/share/RepeatMasker/;
./famdb.py -i Libraries/RepeatMaskerLib.h5 lineage -ad $2 > $2-families.txt;
#查找并导出红酵母及其上祖先节点、其下所有类群的repeat families
./famdb.py -i Libraries/RepeatMaskerLib.h5 families -f embl -ad $2 > $2_ad.embl;
#只查找并导出红酵母及其上祖先节点（可以避免衍生节点物种分布不均匀问题。当有多个物种需要注释，并且需要对比注释的比例及总量时，这样可以控制变量）：
./famdb.py -i Libraries/RepeatMaskerLib.h5 families -f embl -a $2 > $2_a.embl;
#将embl格式转化为fasta格式，方便后期与RepeatModeler2结果合并
cd util
./buildRMLibFromEMBL.pl ../$2_ad.embl > $2_ad.fasta;
./buildRMLibFromEMBL.pl ../$2_a.embl > $2_a.fasta;
#合并RepeatModeler2结果
cd $WORK;
cat /home/data/ssy28/miniconda3/envs/repeat/share/RepeatMasker/util/$2_ad.fasta ./05.genomefeature/repeat/RepeatModeler/mydb-families.fa > ./05.genomefeature/repeat/RepeatModeler/all_$2_final.fasta;
cd 05.genomefeature/repeat;
mkdir RepeatMasker -p;
echo "RepeatMasker";
RepeatMasker -nolow -no_is -norna -parallel $3 -s -gff -lib ./RepeatModeler/all_$2_final.fasta -dir ./RepeatMasker $1 ;
echo "done";
