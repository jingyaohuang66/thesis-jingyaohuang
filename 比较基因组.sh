#比较基因组分析
#转换cds和pep
/home/data/ssy28/miniconda3/envs/his-string-trans/bin/gffread maker.gff3 -g ../genome.fa -x Rhkra.cds.fa -y Rhkra.pep.fa

mkdir data;
cd data;
ln -s */final_gff/*.pep.fa ./
cp ~/workplace/bijiao/maker/*/final_gff/*.pep.fa ./
cp ~/workplace/genome/05.genomefeature/Rhmuc1.pep.fa ./
cp ~/workplace/bijiao/p/*.pep.fa ./
cp ~/workplace/bijiao/p/dio/Rhdio.pep.fa ./
for i in `ls *.fa |cut -d"." -f 1`;do mv ${i}.pep.fa ${i}.fa;done
ln -s ../maker/Rhodotorula_kratochvilovae/final_gff/Rhkra.pep.fa ./Rhkra.fa
orthofinder -f data -S diamond -M msa -T fasttree -t 40
#基因家族统计
在结果目录中的Comparative_Genomics_Statistics文件夹中的Statistics_PerSpecies.tsv
统计出每个copy_number的orthogroups中的基因，以及未匹配的基因
然后统一放到进化树中进行可视化



#串联树
mkdir tree;
cp ../data/OrthoFinder/Results_May03/MultipleSequenceAlignments/SpeciesTreeAlignment.fa ./
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif trimal -in SpeciesTreeAlignment.fa -out SpeciesTreeAlignment_trim.fa -fasta -gt 0.6 -cons 60
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif  raxmlHPC-PTHREADS -T 50 -m PROTGAMMAJTT -f a -p 123 -x 123 -# 1000 -n out -s SpeciesTreeAlignment.fa 1>tree.log 2>tree.err

singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif iqtree -s SpeciesTreeAlignment.fa -m MFP -bb 1000 -nt 50 -pre iqtree.all_FAD.muscle.fasta
#合并树
mkdir tree1;
orthDir=../data/OrthoFinder/Results_May03
cat $orthDir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read aa;do cat $orthDir/Gene_Trees/$aa\_tree.txt |awk '{print $0}' ;done > SingleCopy.trees
#base中的sed
sed -r 's/([(,]{1}[A-z0-9]+)_[^:]*/\1/g' SingleCopy.trees | sed 's/_rho//g' > Astral_input.trees

singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif java -Xms200g -Xmx200g -XX:ParallelGCThreads=48 -jar /opt/ASTRAL/astral.jar -i Astral_input.trees -o Astral_output.tree 2>out.log

#分化时间_基于4d位点做分析，准备cds和pep
mkdir -p timetree/data;
cd timetree/data;
## 拷贝cds和pep文件
cp ../../data/*.fa ./
cp 
## 对cds和pep文件ID加物种简写前缀
for i in `ls *.cds.fa |cut -d "." -f 1` ;do echo "sed 's/^>/>${i}_/' $i.cds.fa > $i.cds.fasta.change" >> add_prefix.sh; echo "sed 's/^>/>${i}_/' $i.pep.fa > $i.pep.fasta.change" >> add_prefix.sh; done
sh add_prefix.sh
cat *.cds.fasta.change > ../all.cds.fasta
cat *.pep.fasta.change > ../all.pep.fasta
#
orthoDir=~/workplace/bijiao/timetree/data/OrthoFinder/Results_May05
cp $orthoDir/Orthogroups/Orthogroups.tsv ./
cp $orthoDir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ./
dos2unix Orthogroups.tsv
# 生成单拷贝基因家族列表文件
## 提取单拷贝基因家族成员方法1
awk '{if(NR==FNR){A[$1]=1}else{ if($1 in A || FNR==1){print $0}} }' \
Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv \
> single_copy.txt
## 对基因ID添加物种前缀
awk '{if(NR==1){ n=split($0,A, "\t")} else { for(i=2; i<= n; i++){ \
printf A[i]"_"$i"\t" }; printf "\n" } }' single_copy.txt | \
sed 's/\s\+$//' > single_copy.txt.change
# 生成cds水平多序列比对结果
## 生成ParaAT线程数文件
echo "60" > proc.txt
#orthofinder的结果会把:变成_,所以要替换下all文件
sed -i 's/:/_/g' all.*
## 生成cds水平比对结果
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif ParaAT_mdf.pl -h single_copy.txt.change -a all.pep.fasta -n all.cds.fasta -o Para_out -p proc.txt
singularity exec ../../software/PhyloTools.sif ParaAT_mdf.pl \
-h single_copy.txt.change \ #基因家族信息，每一行Wie一个家族
-a all.pep.fasta \ #全部蛋白序列
-n all.cds.fasta \ #全部cds序列
-o Para_out \ #输出结果目录
-p proc.txt #指定线程数文件
# 合并单个家族成supergene phylip格式文件
## 去掉基因ID，仅保留物种简写
ls Para_out/*.cds_aln.fasta | while read aa; do \
awk -F "_" '{print $1}' $aa > $aa.change; done
## 对单个家族进行连接
seqkit concat Para_out/*.cds_aln.fasta.change \
> single_copy.cds_msa.fasta
## fasta格式转phylip格式
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif trimal -in single_copy.cds_msa.fasta -out single_copy.cds_msa.phy -phylip_paml
singularity exec ../../software/PhyloTools.sif trimal \
-in single_copy.cds_msa.fasta \
-out single_copy.cds_msa.phy \
-phylip_paml
# 提取4d位点
perl /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/extract_4d_phy.pl single_copy.cds_msa.phy single_copy.cds_msa.4d.phy
# 进行分化时间计算
cp /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/03.Divtime/test4_4d/mcmctree.ctl ./
sed 's/:[^,)(]\+//g' Astral_output.tree | sed 's/)1/)/g' > input.tree


singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif mcmctree mcmctree.ctl
## 以下为手动生成单拷贝基因家族supergene蛋白序列方法
ls Para_out/*.pep_aln |while read aa;do \
awk -F "_" '{print $1}' $aa > $aa.change;done
seqkit concat Para_out/*.pep_aln.change > single_copy.pep_msa.fasta
singularity exec ../../software/PhyloTools.sif trimal \
-in single_copy.pep_msa.fasta \
-out single_copy.pep_msa.phy \
-phylip_paml




基于蛋白做，不推荐
cp ../tree/SpeciesTreeAlignment.fa ./
cp ../tree1/Astral_output.tree ./
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif trimal -in SpeciesTreeAlignment.fa -out supergene.phy -phylip_paml
sed 's/:[^,)(]\+//g' Astral_output.tree | sed 's/)1/)/g' > input.tree


分歧时间树至少要做两次；并且使用脚本评估两次结果的差异度<0.1%
#看趋势趋于平稳
perl -n -e 'my $out; while (s/(.*?(\d\. \d+))//){ my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/； $out .= $info; } $out .= $_; print $out' t2/FigTree.tre > tree01.nex
perl -n -e 'my $out; while (s/(.*?(\d\. \d+))//){ my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/； $out .= $info; } $out .= $_; print $out' t3/FigTree.tre > tree02.nex
perl -e 'while(<>) { if (s/\S*UTREE.*?=\s*//){ s/\s*\[.*?\]//g; print; }}' tree01.nex >tree01.txt
perl -e 'while(<>) { if (s/\S*UTREE.*?=\s*//){ s/\s*\[.*?\]//g; print; }}' tree02.nex >tree02.txt

calculating_branchLength_bias_percentage_of_two_trees.pl --no_normalization_of_total_branch_length tree01.txt tree02.txt > bias_of_2runs.txt
rm tree*
perl -n -e 'my $out; while (s/(.*?(\d\.\d+))//) ( my $info = $1; my $value = $2; my $new = $value * 100;$info =~ s/$value/$new/； $out .= $info; } $out .= $_; print $out' t3/FigTree.tre > tree_abbr.mcmctree

###基因家族收缩扩张

# 拷贝mcmctree分化时间tree文件
cp ../03.Divtime/test4_SingleCopy/FigTree.tre ./
# 准备orthofinder聚类结果统计表
cp ../01.GeneFamilyCluster/data/OrthoFinder/Results_May14/Orthogroups/Orthogroups.GeneCount.tsv ./
# 替换分化时间，以百万年为单位, 即所有枝长乘以100倍
cat FigTree.tre | sed 's/\[[^]]\+\]//g'|awk -F "=" '/UTREE/{print $2}' > input.tree.nwk
sed -e 's/:/\n:/g' -e 's/\([),]\)/\n\1/g' input.tree.nwk | awk '{if($1~/:$/){printf ":"100*$2} else {printf $0}}' | sed 's/\s\+//g' > input.tree
# orthofinder genecount表格转成cafe输入格式
dos2unix Orthogroups.GeneCount.tsv
sed 's/[a-zA-Z0-9]\+$//' Orthogroups.GeneCount.tsv | \awk '{print $1"\t"$0}' > input.tab
# 运行cafe
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif cafe5 --infile input.tab --tree input.tree --output_prefix cafe_ortho --cores 80
singularity exec ../software/PhyloTools.sif cafe5 \
--infile input.tab \
--tree input.tree \
--output_prefix cafe_ortho \
--cores 5
# 生成itol绘图用tree文件和饼图注释信息文件
awk '/TREE OG0000002 = /{print $NF}' cafe_ortho/Base_asr.tre | sed -r 's/_[0-9]+//g ; s/([^)])<[0-9]+>/\1/g' | sed 's/*//g' > out_ortho.tree
Tnum=`grep "TREE " cafe_ortho/Base_asr.tre|wc -l`
echo -e "DATASET_PIECHART
SEPARATOR\tTAB
DATASET_LABEL\ttest
FIELD_COLORS\t#339900\t#FF0033\t#0099FF
FIELD_LABELS\tIncrease\tDecrease\tRemain
DATA" > cafe_itol_ortho.txt
sed -r 's/^([^>]+)<[0-9]+>/\1/g' cafe_ortho/Base_clade_results.txt | \
awk '!/^#/{ p=0.5; s=10; C='$Tnum'-$2-$3; if(!/^</){p=-1}; \
print $1 "\t" p "\t" s "\t" $2 "\t" $3 "\t" C }' >> cafe_itol_ortho.txt

#树加上文字,
echo -e "DATASET_TEXT
SEPARATOR\tTAB
DATASET_LABEL\ttest
DATA" > cafe_itol_ortho_text.txt
sed -r 's/^([^>]+)<[0-9]+>/\1/g' cafe_ortho/Base_clade_results.txt | awk '!/^#/{p=0.5;if(!/^</){p=-1};print $1"\t+"$2"/-"$3"\t" p "\t#ff0000\tnormal\t1\t0"}' >> cafe_itol_ortho_text.txt

#树加上条形图；
利用orthofinder软件生成的结果进行可视化

#基因家族收缩扩张富集分析
le Base_branch_probabilities.tab |awk '$7<0.05||NR==1' |cut -f 1 |sed 1d >genefam.id
awk '{ if (NR==FNR){A[$1] = 1} else{ if(FNR==1 || A[$1]){print $0}}}' genefam.id Base_change.tab >Base_change.tab.genefam
less Base_change.tab.genefam |awk '$7>0' > Base_change.tab.genefam.increase
less Base_change.tab.genefam |awk '$7<0' > Base_change.tab.genefam.decrease
cp ../../time1/data/OrthoFinder/Results_May06/Orthogroups/Orthogroups.tsv ./
dos2unix Orthogroups.tsv
awk '{ if (NR==FNR){A[$1] = 1} else{ if(FNR==1 || A[$1]){print $0}}}' Base_change.tab.genefam.increase Orthogroups.tsv > Base_change.tab.genefam.increase.tab
awk '{ if (NR==FNR){A[$1] = 1} else{ if(FNR==1 || A[$1]){print $0}}}' Base_change.tab.genefam.decrease Orthogroups.tsv > Base_change.tab.genefam.decrease.tab

#富集分析
统计好gene list
然后参考go_kegg_enrich.R


#共线性分析
#1、推断wgd全基因组复制事件，是不是加倍
#2、祖先染色体推断，根据已有的染色体，推断出祖先的染色体
cp ../time1/data1/Rhmucp.pep.fasta.change ./Rhmucp.fasta
cp ../time1/data1/Rhmuc.pep.fasta.change ./Rhmuc.fasta
cp ../time1/data1/Phrho.pep.fasta.change ./Phrho.fasta
cat Rhmuc.fasta Rhmucp.fasta > Rhmuc.Rhmucp.fasta
makeblastdb -in Rhmuc.Rhmucp.fasta  -dbtype prot
blastp -query Rhmuc.Rhmucp.fasta -db Rhmuc.Rhmucp.fasta  -out Rhmuc.Rhmucp_out.tab -outfmt 6 -evalue 1e-10 -num_threads 50
cat Rhmucp.fasta Phrho.fasta > Rhmucp.Phrho.fasta
makeblastdb -in Rhmucp.Phrho.fasta  -dbtype prot
blastp -query Rhmucp.Phrho.fasta -db Rhmucp.Phrho.fasta  -out Rhmucp.Phrho_out.tab -outfmt 6 -evalue 1e-10 -num_threads 50
# 提取mscanx输入的xxx.blast文件
awk '$1~"Rhmucp" && $2~"Rhmucp" && $1 != $2' Rhmuc.Rhmucp_out.tab | sed 's/Rhmucp_//g' >Rhmucp.blast
awk '$1~"Rhmuc_" && $2~"Rhmucp_"' Rhmuc.Rhmucp_out.tab | sed 's/Rhmuc_//g' | sed 's/Rhmucp_//g' >Rhmuc_Rhmucp.blast
awk '$1~"Rhmucp_" && $2~"Rhmuc_"' Rhmuc.Rhmucp_out.tab | sed 's/Rhmucp_//g' | sed 's/Rhmuc_//g' >>Rhmuc_Rhmucp.blast
awk '$1~"Rhmuc_" && $2~"Rhmuc_" && $1 != $2' Rhmuc.Rhmucp_out.tab | sed 's/Rhmuc_//g' >Rhmuc.blast

awk '$1~"Phrho" && $2~"Phrho" && $1 != $2' Rhmucp.Phrho_out.tab | sed 's/Phrho_//g' >Phrho.blast
awk '$1~"Rhmucp" && $2~"Phrho"' Rhmucp.Phrho_out.tab | sed 's/Rhmucp_//g' | sed 's/Phrho_//g' >Rhmucp_Phrho.blast
awk '$1~"Phrho" && $2~"Rhmucp"' Rhmucp.Phrho_out.tab | sed 's/Phrho_//g' | sed 's/Rhmucp_//g' >>Rhmucp_Phrho.blast



# 提取mcscanx输入的xxx.gff文件
cp ~/workplace/bijiao/maker/Rhodotorula_mucilaginosa/final_gff/maker.gff3 ./Rhmuc.gff3
cp ~/workplace/bijiao/maker/Phaffia_rhodozyma/final_gff/maker.gff3 ./Phrho.gff3
cp ~/workplace/genome/05.genomefeature/maker.gff3 ./Rhmucp.gff3
awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' Rhmuc.gff3 | sed 's/ID=//' > Rhmuc.gff
awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' Rhmucp.gff3 | sed 's/ID=//' > Rhmucp.gff
awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' Phrho.gff3 | sed 's/ID=//' > Phrho.gff
cat Rhmuc.gff Rhmucp.gff >Rhmuc_Rhmucp.gff
cat Rhmucp.gff Phrho.gff >Rhmucp_Phrho.gff
#共线性分析
## 种内
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif MCScanX Rhmucp
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif MCScanX Rhmuc
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif MCScanX Phrho
## 种间
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif MCScanX Rhmuc_Rhmucp
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif MCScanX Rhmucp_Phrho

#可视化
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA Rhmuc.gff3 -o Rhmuc.bed
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA Rhmucp.gff3 -o Rhmucp.bed
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA Phrho.gff3 -o Phrho.bed

singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny fromaligns Rhmuc_Rhmucp.collinearity > Rhmuc.Rhmucp.anchors1
awk '{if(NF>=2){print $0"\t1000"} else {print $0}}' Rhmuc.Rhmucp.anchors1 | awk '{print $2"\t"$1"\t"$3 }' |sed  's/\t###/###/g'  > Rhmuc.Rhmucp.anchors
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny screen --minspan=30 --simple Rhmuc.Rhmucp.anchors Rhmuc.Rhmucp.anchors.new

singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny fromaligns Rhmucp_Phrho.collinearity > Rhmucp.Phrho.anchors1
awk '{if(NF>=2){print $0"\t1000"} else {print $0}}' Rhmucp.Phrho.anchors1 > Rhmucp.Phrho.anchors
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny screen --minspan=30 --simple Rhmucp.Phrho.anchors Rhmucp.Phrho.anchors.new

#绘图
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.karyotype seqids2 layout2


#ka/ks的计算
singularity exec /home/data/ssy28/workplace/bijiao/ka-ks/wgd.sif wgd dmd --eval 1e-10  -o dmd_out --nostrictcds Rhmucp.cds.fa Rhmuc.cds.fa
cp dmd_out/*.rbh ./homo_pairs.txt
cat *.cds.fa > input_cds.fasta
cat *.pep.fa > input_pep.fasta
echo "60" > proc.txt
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif ParaAT_mdf.pl -h homo_pairs.txt  -a input_pep.fasta -n input_cds.fasta -p proc.txt  -o  align_out  -m muscle -f axt
cat align_out/*.axt > merge_align.axt
## calculate kaks 
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m YN -i ./merge_align.axt  -o  result.txt
cut -f 1,2,3,4,5 result.txt  > result_KaKs.txt
singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m MA -i ./merge_align.axt  -o  result.txt1
cut -f 1,2,3,4,5 result.txt1  > result_KaKs.txt1

#直系同源基因提取做ka/ks
cut -f 8 single_copy.txt |sed 1d > muc_id.txt
seqkit grep --pattern-file muc_id.txt Rhmuc.cds.fa > Rhmuc_single.cds.fa
seqkit grep --pattern-file muc_id.txt Rhmuc.pep.fa > Rhmuc_single.pep.fa
cut -f 9 single_copy.txt |sed 1d > mucp_id.txt
seqkit grep --pattern-file mucp_id.txt Rhmucp.cds.fa > Rhmucp_single.cds.fa
seqkit grep --pattern-file mucp_id.txt Rhmucp.pep.fa > Rhmucp_single.pep.fa
cd align_out
for i in `ls *.axt |cut -d"." -f 1`;do echo "singularity exec /home/data/ssy28/workplace/genek/基因课比较基因组/PART3_Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m MA -i ${i}.axt  -o  ${i}.txt &">> run.sh;done
split -a3 -l 50 -d run.sh run.






