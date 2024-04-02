#/usr/bin/bash 
if [ $# == 0 ];then 
    echo "This shell is used to repeat genome 
        Usage: sh maker.sh genome species genus thread
    ";  
fi

echo "cp";
WORK=`pwd`;
cp ~/workplace/bijiao/maker/maker_* ./;
cp repeat/$2/all_{$2}_final.fasta;
cp ../protein.fasta ./;
perl -p -i -e 's/^genome=.*/genome={$1}_genome.fa/;
			   s/^protein=.*/protein=protein.fasta/; 
			   s/^rmlib=.*/rmlib=all_{$2}_final.fasta/; 
			   s/^est2genome=0/est2genome=1/;
			   s/^protein2genome=0/protein2genome=1/
			   ' maker_opts.ctl;

mpiexec -n 34 maker &> maker1.log;
mkdir SNAP1;
cd SNAP1;
gff3_merge -d ../{$2}_genome.maker.output/{$2}_genome_master_datastore_index.log;

grep -P "\tmaker\t" {$2}_genome.all.gff > genome.maker.gff3;
gff3_clear.pl --prefix gene genome.maker.gff3 > maker1.gff3;
grep -P "\t" maker1.gff3 >maker.gff3;
/home/data/ssy28/software/SNAP/gff3_to_zff.pl ../{$1}_genome.fa ./maker.gff3 > rho.zff;
/home/data/ssy28/software/SNAP/fathom -validate rho.zff ../{$1}_genome.fa > validate.log;
/home/data/ssy28/software/SNAP/fathom -gene-stats rho.zff ../{$1}_genome.fa > gene-stats.log;
/home/data/ssy28/software/SNAP/fathom -categorize 1000 rho.zff ../{$1}_genome.fa;
/home/data/ssy28/software/SNAP/fathom uni.ann uni.dna -export 1000 -plus;
/home/data/ssy28/software/SNAP/fathom -validate export.ann export.dna;
/home/data/ssy28/software/SNAP/forge export.ann export.dna;
/home/data/ssy28/software/SNAP/hmm-assembler.pl rho . > ../rho1.hmm;
cp maker.gff3 ../genome1.all.gff;
cd $WORK;
cp maker_opts.ctl maker_opts.ctl_backup_rho1;
perl -p -i -e 's/^protein=.*/protein= /; 
			   s/^rmlib=.*/rmlib= /; 
			   s/^est2genome=.*/est2genome=0/;
			   s/^protein2genome=.*/protein2genome=0/;
			   s/^maker_gff=.*/maker_gff=genome1.all.gff/;
			   s/^est_pass=.*/est_pass=1/;
			   s/^protein_pass=.*/protein_pass=1/;
			   s/^rm_pass=.*/rm_pass=1/;
			   s/^snaphmm=.*/snaphmm=rho1.hmm/;
			   s/^pred_stats=.*/pred_stats=1/;
			   s/^alt_splice=.*/alt_splice=0/;
			   s/^keep_preds=.*/keep_preds=1/;
			   ' maker_opts.ctl;

mpiexec -n 34 maker &> maker2.log;
mkdir SNAP2;
cd SNAP2;
gff3_merge -d ../{$2}_genome.maker.output/{$2}_genome_master_datastore_index.log;

grep -P "\tmaker\t" {$2}_genome.all.gff > genome.maker.gff3;
gff3_clear.pl --prefix gene genome.maker.gff3 > maker1.gff3;
grep -P "\t" maker1.gff3 >maker.gff3;
/home/data/ssy28/software/SNAP/gff3_to_zff.pl ../{$1}_genome.fa ./maker.gff3 > rho.zff;
/home/data/ssy28/software/SNAP/fathom -validate rho.zff ../{$1}_genome.fa > validate.log;
/home/data/ssy28/software/SNAP/fathom -gene-stats rho.zff ../{$1}_genome.fa > gene-stats.log;
/home/data/ssy28/software/SNAP/fathom -categorize 1000 rho.zff ../{$1}_genome.fa;
/home/data/ssy28/software/SNAP/fathom uni.ann uni.dna -export 1000 -plus;
/home/data/ssy28/software/SNAP/fathom -validate export.ann export.dna;
/home/data/ssy28/software/SNAP/forge export.ann export.dna;
/home/data/ssy28/software/SNAP/hmm-assembler.pl rho . > ../rho2.hmm;
cp maker.gff3 ../genome2.all.gff;
cd $WORK;
cp maker_opts.ctl maker_opts.ctl_backup_rho2;
perl -p -i -e 's/^protein=.*/protein= /; 
			   s/^rmlib=.*/rmlib= /; 
			   s/^est2genome=.*/est2genome=0/;
			   s/^protein2genome=.*/protein2genome=0/;
			   s/^maker_gff=.*/maker_gff=genome2.all.gff/;
			   s/^est_pass=.*/est_pass=1/;
			   s/^protein_pass=.*/protein_pass=1/;
			   s/^rm_pass=.*/rm_pass=1/;
			   s/^snaphmm=.*/snaphmm=rho2.hmm/;
			   s/^pred_stats=.*/pred_stats=1/;
			   s/^alt_splice=.*/alt_splice=0/;
			   s/^keep_preds=.*/keep_preds=1/;
			   ' maker_opts.ctl;
			   
mpiexec -n 34 maker &> maker3.log;

mkdir SNAP3;
cd SNAP3;
gff3_merge -d ../{$2}_genome.maker.output/{$2}_genome_master_datastore_index.log;

grep -P "\tmaker\t" {$2}_genome.all.gff > genome.maker.gff3;
gff3_clear.pl --prefix gene genome.maker.gff3 > maker1.gff3;
grep -P "\t" maker1.gff3 >maker.gff3;
/home/data/ssy28/software/SNAP/gff3_to_zff.pl ../{$1}_genome.fa ./maker.gff3 > rho.zff;
/home/data/ssy28/software/SNAP/fathom -validate rho.zff ../{$1}_genome.fa > validate.log;
/home/data/ssy28/software/SNAP/fathom -gene-stats rho.zff ../{$1}_genome.fa > gene-stats.log;
/home/data/ssy28/software/SNAP/fathom -categorize 1000 rho.zff ../{$1}_genome.fa;
/home/data/ssy28/software/SNAP/fathom uni.ann uni.dna -export 1000 -plus;
/home/data/ssy28/software/SNAP/fathom -validate export.ann export.dna;
/home/data/ssy28/software/SNAP/forge export.ann export.dna;
/home/data/ssy28/software/SNAP/hmm-assembler.pl rho . > ../rho3.hmm;
cp maker.gff3 ../genome3.all.gff;
cd $WORK;
cp maker_opts.ctl maker_opts.ctl_backup_rho3;
perl -p -i -e 's/^protein=.*/protein= /; 
			   s/^rmlib=.*/rmlib= /; 
			   s/^est2genome=.*/est2genome=0/;
			   s/^protein2genome=.*/protein2genome=0/;
			   s/^maker_gff=.*/maker_gff=genome3.all.gff/;
			   s/^est_pass=.*/est_pass=1/;
			   s/^protein_pass=.*/protein_pass=1/;
			   s/^rm_pass=.*/rm_pass=1/;
			   s/^snaphmm=.*/snaphmm=rho3.hmm/;
			   s/^pred_stats=.*/pred_stats=1/;
			   s/^alt_splice=.*/alt_splice=0/;
			   s/^keep_preds=.*/keep_preds=1/;
			   ' maker_opts.ctl;
			   
mpiexec -n 34 maker &> maker4.log;