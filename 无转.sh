cp ../maker_* ./;
cp ../all_protain.fa ./protein.fasta
cp repeat/RepeatModeler/all_Rhodotorula_mucilaginosa_final.fasta ./repeat.fa
WORK=`pwd`
perl -p -i -e 's/^genome=.*/genome= genome.fa/; 
			   s/^protein=.*/protein=protein.fasta /; 
			   s/^est2genome=.*/est2genome=0/;
			   s/^protein2genome=.*/protein2genome=1/;
			   s/^model_org=.*/model_org=/;
			   s/^rmlib=.*/rmlib=repeat.fa/;
			   ' maker_opts.ctl;

nohup mpiexec -n 40 maker &> maker.log &
cp ../maker_opts.ctl ./
perl -p -i -e 's/^genome=.*/genome= genome.fa/; 
			   s/^protein=.*/protein=protein.fasta /; 
			   s/^est2genome=.*/est2genome=0/;
			   s/^protein2genome=.*/protein2genome=1/;
			   s/^model_org=.*/model_org=/;
			   s/^rmlib=.*/rmlib=repeat.fa/;
			   s/^snaphmm=.*/snaphmm=snap3.hmm/;
			   s/^augustus_species=.*/augustus_species=Rhodotorula_babjevae/;
			   ' maker_opts.ctl;