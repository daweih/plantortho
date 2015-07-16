Data pipeline for Plant Homolog Database.
==========

#English version

## Rawdata download

## Method

### OrthoMCL (Homolog group)

### Phylogenetic Analysis

### Ortholog inference

### dN/dS Analysis

### Annotation using CV

## Reference
* [Rouard,M. et al. (2011) GreenPhylDB v2.0: comparative and functional genomics in plants. Nucleic Acids Research, 39, D1095–102.](http://www.ncbi.nlm.nih.gov/pubmed/?term=20864446)

# 中文版

本分析流程是 [IC4R 数据库](http://www.ic4r.org)的子库，[Plant Homolog Database](http://homolog.ic4r.org) 的数据分析流程。这个分析流程是在[北京基因组研究所](http://www.big.ac.cn)的服务器上运行的。其中的一些细节受限于服务器的集群环境。生物学的部分，则是通用的。如果使用本流程，发表文章请引用 IC4R 数据库的文章（pubmed链接会在文章发表后列出）。

## 原始数据

采用 emsembl release-24 的数据，路径：

`my_working_dir_at_leofs/DataBase/ensembl/release-24`

### 下载数据

`wget -nH -m –ftp-user=anonymous –ftp-password=anonymous ftp://ftp.ensemblgenomes.org/pub/release-24/plants/`

### 数据整理和统计

由于做 homolog 分析要用到蛋白序列比对，所以需要对下载好的蛋白数据做整合。对于一个基因有多个 isoform，因此对应多个蛋白的，选择其中一条进入下游分析（*.pep.all.fa.gz）。聚类的依据是 ensembl 的注释信息，例子如下：

`>OS01T0100100-01 pep:known chromosome:IRGSP-1.0:1:2983:10815:1 gene:OS01G0100100 transcript:OS01T0100100-01 description:"Note\x3dRabGAP/TBC domain containing protein., Transcript_evidence\x3dAK242339 (DDBJ, antisense transcript), ORF_evidence\x3dQ655M0 (UniProt), NIAS_FLcDNA\x3dJ075199P03,"`
`MSSAAGQDNGDTAGDYIKWMCGAGGRAGGAMANLQRGVGSLVRDIGDPCLNPSPVKGSKM`

其中，gene:OS01G0100100 标明这个蛋白对应的基因，只有一个转录本。而 gene:OS08G0564300 有 6 条转录本，所有需要写一个 perl 脚本选择一个。其他几个的 ID 也同时保存，用于后面数据的搜索用。

## Method

这个部分是数据的分析流程。可以分成两大部分。OrthoMCL 部分（OrthoMCL），脚本生成/提交部分。由于这个研究使用了集群，所有任务采用 dsub 命令提交。

### OrthoMCL (Homolog group)

### Phylogenetic Analysis

#### mafft

#### trimal

#### PhyML

### Ortholog inference

#### RIO

#### GSDI


### dN/dS Analysis

### Annotation using cv

For more information, visit:
- [Plant Homolog Database](http://homolog.ic4r.org)
- [DaWei Huang's Homepage at Zhang Zhang Lab](http://cbb.big.ac.cn/Dawei_Huang)
- [Li Yang's Home page at Zhang Zhang Lab](http://cbb.big.ac.cn/Li_Yang)
- [Xingjian Xu's Home page at Zhang Zhang Lab](http://cbb.big.ac.cn/Xingjian_Xu)
