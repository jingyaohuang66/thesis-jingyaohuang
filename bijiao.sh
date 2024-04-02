cp ../maker_* ./
cp ../protein.fasta ./
cp star/03.star/sample.gff ./
cp repeat/RepeatModeler/mydb-families.fa ./repeat.fa
WORK=`pwd`
perl -p -i -e 's/^genome=.*/genome= genome.fa/; 
			   s/^protein=.*/protein=protein.fasta /; 
			   s/^est2genome=.*/est2genome=1/;
			   s/^protein2genome=.*/protein2genome=1/;
			   s/^est_gff=.*/est_gff=sample.gff/;
			   s/^model_org=.*/model_org=/;
			   s/^rmlib=.*/rmlib=repeat.fa/;
			   ' maker_opts.ctl;
mpiexec -n 40 maker &> maker.log;
mkdir SNAP1
cd SNAP1

gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log

grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3
gff3_clear.pl --prefix gene genome.maker.gff3 > maker.gff3
grep -P "\t" maker.gff3 >maker1.gff3
/home/data/ssy28/software/SNAP/gff3_to_zff.pl ../genome.fa ./maker1.gff3 > rho.zff
/home/data/ssy28/software/SNAP/fathom -validate rho.zff ../genome.fa > validate.log
/home/data/ssy28/software/SNAP/fathom -gene-stats rho.zff ../genome.fa > gene-stats.log
/home/data/ssy28/software/SNAP/fathom -categorize 1000 rho.zff ../genome.fa

/home/data/ssy28/software/SNAP/fathom uni.ann uni.dna -export 1000 -plus
/home/data/ssy28/software/SNAP/fathom -validate export.ann export.dna
/home/data/ssy28/software/SNAP/forge export.ann export.dna
hmm-assembler.pl rho . > ../snap1.hmm

cp genome.all.gff ../genome1.all.gff;
cd ..;
cp maker_opts.ctl maker_opts.ctl_backup;
perl -p -i -e 's/^protein=.*/protein= /; 
			   s/^rmlib=.*/rmlib= /;
			   s/^est_gff=.*/est_gff= /;
			   s/^est2genome=.*/est2genome=0/;
			   s/^protein2genome=.*/protein2genome=0/;
			   s/^maker_gff=.*/maker_gff=genome1.all.gff/;
			   s/^est_pass=.*/est_pass=1/;
			   s/^protein_pass=.*/protein_pass=1/;
			   s/^rm_pass=.*/rm_pass=1/;
			   s/^snaphmm=.*/snaphmm=snap1.hmm/;
			   s/^pred_stats=.*/pred_stats=1/;
			   s/^alt_splice=.*/alt_splice=0/;
			   s/^keep_preds=.*/keep_preds=1/;
			   ' maker_opts.ctl;
mpiexec -n 40 maker &> maker1.log;

mkdir SNAP2;
cd SNAP2;
gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log;
grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3;
gff3_clear.pl --prefix gene genome.maker.gff3 > maker.gff3;
grep -P "\t" maker.gff3 >maker2.gff3;
/home/data/ssy28/software/SNAP/gff3_to_zff.pl ../genome.fa ./maker2.gff3 > rho.zff
/home/data/ssy28/software/SNAP/fathom -validate rho.zff ../genome.fa > validate.log
/home/data/ssy28/software/SNAP/fathom -gene-stats rho.zff ../genome.fa > gene-stats.log
/home/data/ssy28/software/SNAP/fathom -categorize 1000 rho.zff ../genome.fa

/home/data/ssy28/software/SNAP/fathom uni.ann uni.dna -export 1000 -plus
/home/data/ssy28/software/SNAP/fathom -validate export.ann export.dna
/home/data/ssy28/software/SNAP/forge export.ann export.dna
hmm-assembler.pl rho . > ../snap2.hmm
cp genome.all.gff ../genome2.all.gff;
cd ..;
cp maker_opts.ctl maker_opts.ctl_backup_snap1;
perl -p -i -e 's/^maker_gff=.*/maker_gff=genome2.all.gff/;
			   s/^snaphmm=.*/snaphmm=snap2.hmm/;
			   ' maker_opts.ctl;
		
mpiexec -n 40 maker &> maker2.log;

mkdir SNAP3;
cd SNAP3;
gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log;
grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3;
gff3_clear.pl --prefix gene genome.maker.gff3 > maker.gff3;
grep -P "\t" maker.gff3 >maker3.gff3;
/home/data/ssy28/software/SNAP/gff3_to_zff.pl ../genome.fa ./maker3.gff3 > rho.zff
/home/data/ssy28/software/SNAP/fathom -validate rho.zff ../genome.fa > validate.log
/home/data/ssy28/software/SNAP/fathom -gene-stats rho.zff ../genome.fa > gene-stats.log
/home/data/ssy28/software/SNAP/fathom -categorize 1000 rho.zff ../genome.fa

/home/data/ssy28/software/SNAP/fathom uni.ann uni.dna -export 1000 -plus
/home/data/ssy28/software/SNAP/fathom -validate export.ann export.dna
/home/data/ssy28/software/SNAP/forge export.ann export.dna
hmm-assembler.pl rho . > ../snap3.hmm
cp genome.all.gff ../genome3.all.gff;
cd ..;
cp maker_opts.ctl maker_opts.ctl_backup_snap2;
#genemark
cp star/03.star/gmhmm.mod ./
perl -p -i -e 's/^maker_gff=.*/maker_gff=genome3.all.gff/;
			   s/^snaphmm=.*/snaphmm=snap3.hmm/;
			   s/^gmhmm=.*/gmhmm=gmhmm.mod/;
			   ' maker_opts.ctl;
mpiexec -n 40 maker &> maker3.log;
cp maker_opts.ctl maker_opts.ctl_backup_genemark;

#augustus
mkdir augustus
cd augustus
gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log
grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3
gff3_clear.pl --prefix gene genome.maker.gff3 > maker.gff3
grep -P "\t" maker.gff3 >maker4.gff3
cp ../genome.fa ./
gff2gbSmallDNA.pl maker4.gff3 genome.fa 1000 trans.gene.raw.gb
new_species.pl  --species=for_bad_genes_removing
#训练，捕捉错误
etraining --species=for_bad_genes_removing --stopCodonExcludedFromCDS=false trans.gene.raw.gb 2> train.err
#过滤错误
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes_list.txt
filterGenes.pl badgenes_list.txt trans.gene.raw.gb > trans.genes.gb
grep 'BASE COUNT' trans.genes.gb |wc -l
randomSplit.pl trans.genes.gb 200
new_species.pl --species=Phaffia_rhodozyma
etraining --species=Phaffia_rhodozyma trans.genes.gb.train > train.out

perl -e 'open IN, "train.out"; while (<IN>) { $tag = $1  if m/tag:.*\((.*)\)/; $taa = $1  if m/taa:.*\((.*)\)/;$tga = $1  if m/tga:.*\((.*)\)/; } while (<>) { s#/Constant/amberprob.x#/Constant/amberprob   $tag#; s#/Constant/ochreprob.*#/Constant/ochreprob   $taa#; s#/Constant/opalprob.*#/Constant/opalprob   $tga#;print }' ~/miniconda3/envs/maker/config/species/Phaffia_rhodozyma/Phaffia_rhodozyma_parameters.cfg > 11
mv 11 ~/miniconda3/envs/maker/config/species/Phaffia_rhodozyma/Phaffia_rhodozyma_parameters.cfg

augustus --species=Phaffia_rhodozyma trans.genes.gb.test | tee firsttest.out
grep -A 22 Evaluation firsttest.out #查看报告
grep 'LOCUS' trans.genes.gb.train |wc
#由于trans.genes.gb.train有7116个，所以
randomSplit.pl trans.genes.gb.train 316
ln -s trans.genes.gb.train.train training.gb.onlytrain
optimize_augustus.pl --species=Phaffia_rhodozyma --rounds=5 --cpus=48 --kfold=55 --onlytrain=training.gb.onlytrain  trans.genes.gb.train.test > optimize.out
etraining --species=Phaffia_rhodozyma trans.genes.gb.train
augustus --species=Phaffia_rhodozyma trans.genes.gb.test | tee secondtest.out

#一般CRF久而且精度不太好
#比较CRF与非CRF的精确性，注意augustus在3.3.3版本才能很好的运行CRF
cd ~/miniconda3/envs/maker/config/species/Phaffia_rhodozyma
cp Phaffia_rhodozyma_exon_probs.pbl Phaffia_rhodozyma_exon_probs.pbl.withoutCRF
cp Phaffia_rhodozyma_igenic_probs.pbl Phaffia_rhodozyma_igenic_probs.pbl.withoutCRF
cp Phaffia_rhodozyma_intron_probs.pbl Phaffia_rhodozyma_intron_probs.pbl.withoutCRF
cd -
etraining --species=Phaffia_rhodozyma trans.genes.gb.train --CRF=1
augustus --species=Phaffia_rhodozyma trans.genes.gb.test | tee secondtest.out.withoutCRF
#如果CRF的精确性比非CRF的精确性低，那么就备份回去
cd ~/miniconda3/envs/maker/config/species/Phaffia_rhodozyma
cp Phaffia_rhodozyma_exon_probs.pbl Phaffia_rhodozyma_exon_probs.pbl.withCRF
cp Phaffia_rhodozyma_igenic_probs.pbl Phaffia_rhodozyma_igenic_probs.pbl.withCRF
cp Phaffia_rhodozyma_intron_probs.pbl Phaffia_rhodozyma_intron_probs.pbl.withCRF
cd -

#last
cp ../maker_opts.ctl ./
perl -p -i -e 's/^genome=.*/genome= genome.fa/; 
			   s/^protein=.*/protein=protein.fasta /; 
			   s/^est2genome=.*/est2genome=1/;
			   s/^protein2genome=.*/protein2genome=1/;
			   s/^est_gff=.*/est_gff=sample.gff/;
			   s/^model_org=.*/model_org=/;
			   s/^rmlib=.*/rmlib=repeat.fa/;
			   s/^snaphmm=.*/snaphmm=snap3.hmm/;
			   s/^augustus_species=.*/augustus_species=Rhodotorula_kratochvilovae1/;
			   s/^gmhmm=.*/gmhmm=gmhmm.mod/;
			   ' maker_opts.ctl;
			   	   
mpiexec -n 34 maker &> maker4.log;			   
mkdir final_gff;
cd final_gff;			   
gff3_merge -d ../genome.maker.output/genome_master_datastore_index.log;
grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3;
gff3_clear.pl --prefix gene genome.maker.gff3 > maker1.gff3;
grep -P "\t" maker1.gff3 >maker.gff3;			   
			   
			   
			   
			   
			   
			   
			   
			   
			   
			   
			   
#转蛋白出来
/home/data/ssy28/miniconda3/envs/his-string-trans/bin/gffread maker.gff3 -g ../genome.fa -y protein.fa
