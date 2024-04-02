mkdir SNAP1
cd SNAP1

gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log

grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3
gff3_clear.pl --prefix gene genome.maker.gff3 > maker.gff3
grep -P "\t" maker.gff3 >maker1.gff3
/pub/software/SNAP-master/gff3_to_zff.pl ../genome.fasta ./maker1.gff3 > rho.zff
/pub/software/SNAP-master/fathom -validate rho.zff ../genome.fasta > validate.log
/pub/software/SNAP-master/fathom -gene-stats rho.zff ../genome.fasta > gene-stats.log

#maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 genome.all.gff
#fathom genome.ann genome.dna -gene-stats &> gene-stats.log

#fathom genome.ann genome.dna -validate &> validate.log

#perl -ne 'print "$1\n" if /.*:\s+(\S+)\s+OK/' validate.log >zff2keep.txt

#filterGenes.pl zff2keep.txt rho.zff > tmp 
#mv tmp good.zff

/pub/software/SNAP-master/fathom -categorize 1000 rho.zff ../genome.fasta

/pub/software/SNAP-master/fathom uni.ann uni.dna -export 1000 -plus
fathom -validate export.ann export.dna

forge export.ann export.dna

hmm-assembler.pl rho . > ../rho1.hmm

cp genome.all.gff ../genome1.all.gff
cd ..
#导入maker中。重复两三次

cp maker_opts.ctl maker_opts.ctl_backup_rho1

vi maker_opts.ctl 

# maker_gff= pyu_rnd1.all.gff
# est_pass=1 # 使用第一轮的EST比对结果
# protein_pass=1 #使用第一轮的protein比对结果
# rm_pass=1 # 使用gff文件中的repeats
# snaphmm=pyu1.hmm
# est= # 删除est文件，这一步不需要再跑EST比对了
# protein= # 同上
# model_org= #同上
# rmlib= # 同上
# repeat_protein= #同上
# est2genome=0 # 不需要再构建基于EST证据的基因模型
# protein2genome=0 #同上
# pred_stats=1 #report AED stats
# alt_splice=0 # 0: keep one isoform per gene; 1: identify splicing variants of the same gene
# keep_preds=1 # keep genes even without evidence support, set to 0 if no

nohup mpiexec -n 34 maker &> maker.log &
pstree -ap [mpiexec进程ID]更细致的了解MAKER的进程树
##如果只是用snap预测基因
#snap species.hmm  species.genome > species.zff

#zff2gff3.pl species.zff > species.gff3

mkdir SNAP2
cd SNAP2


gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log

maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 genome.all.gff

fathom genome.ann genome.dna -gene-stats &> gene-stats.log

fathom genome.ann genome.dna -validate &> validate.log

fathom -categorize 1000 genome.ann genome.dna

fathom uni.ann uni.dna -export 1000 -plus

forge export.ann export.dna

hmm-assembler.pl rho . > ../rho2.hmm

mkdir SNAP3
cd SNAP3


gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log

maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 genome.all.gff

fathom genome.ann genome.dna -gene-stats &> gene-stats.log

fathom genome.ann genome.dna -validate &> validate.log

fathom -categorize 1000 genome.ann genome.dna

fathom uni.ann uni.dna -export 1000 -plus

forge export.ann export.dna

hmm-assembler.pl rho . > ../rho3.hmm

