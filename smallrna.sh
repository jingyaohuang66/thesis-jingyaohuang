#smallRNA分析
#data下载
#软件安装
#安装singularity(root安装)
# Ensure repositories are up-to-date
sudo apt-get update
# Install debian packages for dependencies
sudo apt-get install -y \
    build-essential \
    libseccomp-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup \
    curl wget git
export GOVERSION=1.17.3 OS=linux ARCH=amd64  # change this as you need
wget -O /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${GOVERSION}.${OS}-${ARCH}.tar.gz
sudo tar -C /usr/local -xzf /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz
echo 'export PATH=$PATH:/usr/local/go/bin' >> /etc/profile.d/env.sh
source ~/.bashrc
curl -sSfL https://raw.githubusercontent.com/golangci/golangci-lint/master/install.sh | sh -s -- -b $(go env GOPATH)/bin v1.43.0
echo 'export PATH=$PATH:$(go env GOPATH)/bin' >> /etc/profile.d/env.sh
source ~/.bashrc
git clone https://github.com/hpcng/singularity.git
cd singularity
./mconfig
cd ./builddir
make
sudo make install
singularity --version
cd ..
./mconfig -b ./buildtree -p /usr/local
#安装singularity(非root安装)
mamba create -n singularity singularity -y

#构建samllRNA分析容器
singularity pull --arch amd64 library://library/default/ubuntu:20.04
singularity build --sandbox smallRNA ubuntu_20.04.sif
singularity shell --writable --fakeroot smallRNA
#安装mamba下载安装脚本到本地
apt-get update
apt-get upgrade -y
apt-get dist-upgrade -y
apt install wget vim make gcc g++ axel cutadapt zip git libpdf-api2-perl -y
cd /opt
#安装软件
#r
apt update -qq
apt install --no-install-recommends software-properties-common dirmngr -y
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt install --no-install-recommends r-base
#R包
## tidyverse 用于数据处理
install.packages("tidyverse")
## BiocManager 用于安装其他R包
install.packages("BiocManager")
## 安装绘图相关R包
install.packages("ggplot2")
install.packages("cowplot")
install.packages("ggrepel")
install.packages("ggsci")
## pheatmap用于绘制热图
install.packages("pheatmap")
## VennDiagram用于绘制venn图
install.packages("VennDiagram")
## PCAtools 用于PCA分析
BiocManager::install("PCAtools")
## DESeq2和edgeR均为差异表达分析软件
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
## clusterProfiler用于富集分析
BiocManager::install("clusterProfiler")
## 富集分析结果绘图R包
BiocManager::install("enrichplot")
## AnnotationForge 用于orgdb构建
BiocManager::install("AnnotationForge")
#安装mamba
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh
# 添加执行权限
chmod +x Mambaforge-pypy3-Linux-x86_64.sh
# 运行脚本
./Mambaforge-pypy3-Linux-x86_64.sh
#添加.condarc
cd conda
vim .condarc
channels:
  - main
  - free
  - conda-forge
  - r
  - bioconda
  - defaults
#添加到环境变量中
echo 'export PATH=/opt/conda/bin:$PATH' >> /.singularity.d/env/99-base.sh
#安装程序
mamba install -y fastqc bowtie seqkit RepeatMasker RepeatModeler seqtk bedtools gffread perl-bioperl targetfinder mirdeep2 r-tidyverse r r-ggplot2 r-cowplot r-ggrepel r-ggsci r-pheatmap r-venndiagram bioconductor-pcatools bioconductor-deseq2 bioconductor-edger bioconductor-clusterprofiler bioconductor-enrichplot bioconductor-annotationforge r-argparser r-formattable r-seqinr r-pkgbuild
#echo 'export PATH=$PATH:/opt/conda/envs/small/bin' >> /etc/profile.d/env.sh
#更新perl-bioperl
#安装trim-galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
# Run Trim Galore
./TrimGalore-0.6.10/trim_galore
echo 'export PATH=$PATH:/opt/TrimGalore-0.6.10' >> /.singularity.d/env/99-base.sh
#安装miRDeep-P2
mkdir miRDeep-P2
cd miRDeep-P2
axel -n 10 -o ./ https://sourceforge.net/projects/mirdp2/files/version%201.1.4/miRDP2-v1.1.4.tar.gz
axel -n 10 -o ./ https://sourceforge.net/projects/mirdp2/files/latest_version/ncRNA_rfam.tar.gz
tar -zxvf  *.gz
rm *.gz
cd 1.1.4/scripts/index/
ln -s ../../../ncRNA_rfam.fa
bowtie-build ncRNA_rfam.fa rfam_index
echo 'export PATH=$PATH:/opt/miRDeep-P2/1.1.4/scripts' >> /.singularity.d/env/99-base.sh
echo 'export PATH=$PATH:/opt/miRDeep-P2/1.1.4' >> /.singularity.d/env/99-base.sh

#安装miRDeep
## 手动安装
## 软件下载
wget https://github.com/rajewsky-lab/mirdeep2/archive/v0.1.3.tar.gz
## 解压压缩包生成mirdeep2-0.1.3目录
tar -zxvf v0.1.3.tar.gz
cd mirdeep2-0.1.3
source ~/.bashrc
perl install.pl
## 添加环境
echo 'export PATH=$PATH:/opt/mirdeep2-0.1.3/bin' >> /.singularity.d/env/99-base.sh

#安装miRanda
cd /opt
axel -n 10 -o ./ http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz
tar -zxvf miRanda-aug2010.tar.gz
cd miRanda-3.3a
./configure
make
make install

#安装RNAhybrid
axel -n 10 -o ./ https://bibiserv.cebitec.uni-bielefeld.de/applications/rnahybrid/resources/downloads/RNAhybrid-2.1.2.tar.gz
tar -zxvf RNAhybrid-2.1.2.tar.gz
cd RNAhybrid-2.1.2
./configure
make
make install

#安装targetscan
cd /opt
mkdir targetscan_70
cd targetscan_70
wget https://www.targetscan.org/vert_80/vert_80_data_download/targetscan_70.zip
unzip targetscan_70.zip
rm targetscan_70.zip
echo 'export PATH=$PATH:/opt/targetscan_70' >> /.singularity.d/env/99-base.sh

#安装psRobot
cd /opt
wget http://omicslab.genetics.ac.cn/psRobot/program/WebServer/psRobot_v1.2.tar.gz
tar -zxvf psRobot_v1.2.tar.gz
cd psRobot_v1.2
./configure
make
make install
#orgde
cd /opt
git clone http://git.genek.cn:3333/zhxd2/emcp
mamba install r-argparser r-formattable r-seqinr r-pkgbuild RepeatModeler RepeatMasker
echo 'export PATH=$PATH:/opt/emcp' >> /.singularity.d/env/99-base.sh
cd /opt
rm Mambaforge-pypy3-Linux-x86_64.sh trim_galore.tar.gz miRanda-aug2010.tar.gz RNAhybrid-2.1.2.tar.gz psRobot_v1.2.tar.gz v0.1.3.tar.gz
conda clean --all
exit
#打包沙箱
singularity build smallRNA.sif smallRNA


#
mamba create -n small fastqc trim-galore bowtie seqkit seqtk bedtools gffread perl-bioperl targetfinder mirdeep2 r-tidyverse r r-ggplot2 r-cowplot r-ggrepel r-ggsci r-pheatmap r-venndiagram bioconductor-pcatools bioconductor-deseq2 bioconductor-edger bioconductor-clusterprofiler bioconductor-enrichplot bioconductor-annotationforge


WORK=`pwd`
mkdir data/fastq -p
cd data/fastq
#下载smallrna测序得到的fastq序列
#改名
bash mv.sh 


#pmiren数据库下载
cd $WORK
mkdir data/pmiren -p
cd data/pmiren
#下载自己指定物种
wget -np -nH -r https://www.pmiren.com/ftp-download/Citrus_sinensis_Csi/


#miRNA数据准备
#下载mirBase数据库
cd $WORK
mkdir data/mirBase -p
cd data/mirBase
axel -n 100 -o ./  https://mirbase.org/download//hairpin.fa
axel -n 100 -o ./  https://mirbase.org/download/mature.fa
axel -n 100 -o ./  https://mirbase.org/download/miRNA.dat
#提取自己物种的mirna(如果没有无需提取)
## 基于前缀提取鲤鱼miRNA数据
spc=ccr
mature=../../data/mirBase/mature.fa
hairpin=../../data/mirBase/hairpin.fa
## step1 提取miRNA成熟体序列
### 提取成熟体序列id
grep "^>$spc" $mature | sed 's/^>//' | \
awk '{print $1}' > mature.$spc.id
### 基于id提取成熟体序列
seqkit grep -f mature.$spc.id $mature | \
awk '{print $1}' > mature.$spc.fa.tmp
### 将序列中的U替换成T
sed '/^>/!s/U/T/gi' mature.$spc.fa.tmp > mature.$spc.fa
## step1 提取miRNA前体序列
### 提取前体序列id
grep "^>$spc" $hairpin | sed 's/^>//' | \
awk '{print $1}' > hairpin.$spc.id
### 基于id提取序列
seqkit grep -f hairpin.$spc.id $hairpin | \
awk '{print $1}' > hairpin.$spc.fa.tmp
### 替换U成T
sed '/^>/!s/U/T/gi' hairpin.$spc.fa.tmp > hairpin.$spc.fa
### 构建bowtie index
bowtie-build hairpin.$spc.fa hairpin.$spc

#提取近缘物种的mirna
## 提取miRBase数据库中所有Teleostei分类的物种简写名称，去除鲤鱼
grep Teleostei ../../data/mirBase/organisms.txt | \
awk '{print $1}'|grep -v ccr > Teleostei_pre.txt
## 基于物种简写名称提取成熟体序列id
cat Teleostei_pre.txt | while read line ; do \
grep "^>$line" ../../data/mirBase/mature.fa ; done | \
awk '{print $1}' |sed 's/^>//' > mature.Teleostei.id
## 基于序列id提取成熟体序列
seqkit grep -f mature.Teleostei.id ../../data/mirBase/mature.fa | \
awk '{print $1}' > mature.Teleostei.fa.tmp
## 将U替换成T，方便后续比对分析
sed '/^>/!s/U/T/gi' mature.Teleostei.fa.tmp > mature.Teleostei.fa
## 基于物种简写名称提取前体序列id
cat Teleostei_pre.txt | while read line ; do \
grep "^>$line" ../../data/mirBase/hairpin.fa ; done | \
awk '{print $1}' |sed 's/^>//' > hairpin.Teleostei.id
## 基于序列id提取前体序列
seqkit grep -f hairpin.Teleostei.id ../../data/mirBase/hairpin.fa | \
awk '{print $1}' > hairpin.Teleostei.fa.tmp
## 将U替换成T
sed '/^>/!s/U/T/gi' hairpin.Teleostei.fa.tmp > hairpin.Teleostei.fa


#ncRNA 数据准备
cd $WORK
mkdir data/Rfam -p
cd data/Rfam
#下载
axel -n 100 -o ./ https://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/database_files/family.txt.gz
wget -np -nH -r https://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/fasta_files/
#
gunzip family.txt.gz
awk '{print "axel -n 10 -o ./ https://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/fasta_files/"$1".fa.gz"}' family.txt >download.sh
nohup bash download.sh &
#拆分ncRNA数据库，按照mirna进行拆分
## 将family文件中的type和家族信息进行提取
cut -f 1,2,19 family.txt > family.txt1
## 整理Rfam中非miRNA家族数据
### 删除其中的miRNA家族
awk '$0 !~ /Gene; miRNA/' family.txt1 > family.txt.rmMIR
### 对家族列表中每一个家族fasta序列文件
### 在其序列id前添加家族编号，方便后续分析
cut -f 1 family.txt.rmMIR | while read line ;
do
zcat $line.fa.gz| sed "s/^>/>$line:/"
done > Rfam_rmMIR.fasta
## 构建bowtie数据库
bowtie-build Rfam_rmMIR.fasta Rfam_rmMIR.fasta
## 整理Rfam中miRNA家族数据
### 提取miRNA家族系信息
awk '$0 ~ /Gene; miRNA/' family.txt1 > family.txt.MIR
### 循环处理miRNA家族序列，添加家族编号并合并
cut -f 1 family.txt.MIR | while read line ;
do
zcat $fastadir/$line.fa.gz| sed "s/^>/>$line:/"
done > Rfam_MIR.fasta
### 构建Rfam blast 库
makeblastdb -in Rfam_MIR.fasta -dbtype nucl

#基因组准备
cd $WORK
mkdir data/genome -p
cd data/genome
#浏览器下载http://citrus-hzau-edu-cn.vpn.scau.edu.cn:443/download.php，然后传到服务器中
#处理文件


































