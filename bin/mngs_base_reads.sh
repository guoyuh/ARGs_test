#-rw-rw-r-- 1 huanggy huanggy 3.3G 5月  26 18:49 /mnt/home/huanggy/raw_fq/PRJNA626477/SRR14693230_1.fastq
#-rw-rw-r-- 1 huanggy huanggy 3.3G 5月  26 18:49 /mnt/home/huanggy/raw_fq/PRJNA626477/SRR14693230_2.fastq

##1.1 filter
fastp -i reads.1.fq.gz -I reads.2.fq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z 4 -q 20 -u 30 -n 5 --length_required 80 -w 16 
workdir=/mnt/home/huanggy/project/PRJNA626477
cd $workdir && mkdir 01.filter
fastp=/mnt/db/software/fastp/fastp
tail -n+2 metadata.txt |cut -f1 | sort | uniq |rush -j 4 \
"$fastp \
-i /mnt/home/huanggy/raw_fq/PRJNA626477/{}_1.fastq \
-I /mnt/home/huanggy/raw_fq/PRJNA626477/{}_2.fastq \
-o 01.filter/{}.R1.fq.gz -O 01.filter/{}.R2.fq.gz \
-q 20 -u 30 --length_required 80 --thread 4 \
--json 01.filter/{}.json \
--html 01.filter/{}.html"

## 1.2 host remove
host=/mnt/db/kneaddata/human_genome/hg37dec_v0.1
tail -n+2 metadata.txt |cut -f1 | sort | uniq |rush -j 8 \
"/mnt/project/tools/bowtie2-2.4.4/bowtie2 --mm --sensitive --threads 8 \
-x ${host} \
-1 01.filter/{}.R1.fq.gz \
-2 01.filter/{}.R2.fq.gz \
--un-conc-gz 01.filter/{}.clean.fq.gz \
--no-unal --no-head -k 1 -S 01.filter/{}.host.sam \
1> 01.filter/{}.bowtie.log 2>&1 "

#########################qc stat##############################################
tail -n+2 metadata.txt |cut -f1 | sort | uniq |rush -j 4 \
"python ~/idseq_wkf/scripts/qc_paired.py  01.filter/{}.json 01.filter/{}.clean.fq.1.gz 01.filter/{}.clean.fq.2.gz {} 01.qc/{}.qc.xls"

###
## 2.1双端合并为单个文件
mkdir -p temp/concat
#for i in `cat metadata.txt | cut -f1  | sort | uniq `;do zcat  01.filter/${i}.clean.fq.*.gz > temp/concat/${i}.fq ; done
tail -n+2 metadata.txt | cut -f1  | sort | uniq | rush -j 8  "zcat 01.filter/{}.clean.fq.*.gz > temp/concat/{}.fq"

## 2.2 HUMAnN2计算物种和功能组成
多样本并行计算，本次测试数据约耗时2小时
mkdir -p temp/humann
nohup tail -n +2 metadata.txt | cut -f1  | sort | uniq |rush -j 4 'humann --threads 16 --input temp/concat/{1}.fq  --output temp/humann/' 1> humann.o 2> humann.e &

## 2.3 物种组成表
### 2.3.1 样品结果合并
mkdir -p 02.metaphlan3
merge_metaphlan_tables.py temp/humann/*/*_metaphlan_bugs_list.tsv | sed 's/_metaphlan_bugs_list//g' > 02.metaphlan3/taxonomy.tsv
tail -n +2 02.metaphlan3/taxonomy.tsv | head -n 1 | sed 's/clade_name/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/' > 02.metaphlan3/taxonomy.spf.header
#提取出门纲目科属种
#grep -E "s__|clade" taxonomy.tsv | sed 's/^.*s__//g' | cut -f1,3-33  > taxonomy.S.tsv
#python /mnt/home/huanggy/idseq_wkf/scripts/metaphlan3_out_split.py
### 2.3.2 转换为stamp的spf格式

metaphlan_to_stamp.pl 02.metaphlan3/taxonomy.tsv > 02.metaphlan3/taxonomy.spf.tmp
cat 02.metaphlan3/taxonomy.spf.header 02.metaphlan3/taxonomy.spf.tmp | sed '2d' > 02.metaphlan3/taxonomy.spf


## 2.4 功能组成分析
### 2.4.1 功能组成合并、标准化和分层

合并通路丰度(pathabundance)，含功能和对应物种组成。
可选基因家族(genefamilies 太多)，通路覆盖度(pathcoverage)。
注意看屏幕输出`# Gene table created: 03.humann/pathabundance.tsv`

mkdir -p 03.humann
humann_join_tables \
--input temp/humann \
--file_name pathabundance \
--output 03.humann/pathabundance.tsv
# 样本名调整：删除列名多余信息
head 03.humann/pathabundance.tsv
sed -i 's/_Abundance//g' 03.humann/pathabundance.tsv
# 预览和统计
head 03.humann/pathabundance.tsv
/mnt/project/tools/utilis/csvtk/csvtk -t stat 03.humann/pathabundance.tsv

#标准化为相对丰度relative abundance(1)或百万比cpm(1,000,000)
humann_renorm_table \
--input 03.humann/pathabundance.tsv \
--units relab \
--output 03.humann/pathabundance_relab.tsv
head -n5 03.humann/pathabundance_relab.tsv

#分层结果：功能和对应物种表(stratified)和功能组成表(unstratified)
humann_split_stratified_table \
--input 03.humann/pathabundance_relab.tsv \
--output 03.humann/

########################################################4.1  kraken2 物种组成##############################################
mkdir -p temp/kraken2
kraken2=/mnt/project/tools/kraken2/kraken2
db=/mnt/data/kraken2_mask
#并行
tail -n+2 metadata.txt|cut -f1| rush -j 4 \
"/mnt/project/tools/kraken2/kraken2 --db /mnt/data/kraken2_mask --paired 01.filter/{1}.clean.fq.*.gz \
--threads 16 --use-names --report-zero-counts \
--report temp/kraken2/{1}.report \
--output temp/kraken2/{1}.output"


mkdir -p 04.kraken2
#使用krakentools转换report为mpa格式
for i in `tail -n+2 metadata.txt|cut -f1`;do
~/idseq_wkf/KrakenTools/kreport2mpa.py -r temp/kraken2/${i}.report \
--display-header \
-o temp/kraken2/${i}.mpa;done

#合并样本为表格
# 输出结果行数相同，但不一定顺序一致，要重新排序
for i in `tail -n+2 metadata.txt|cut -f1`;do
tail -n+2 temp/kraken2/${i}.mpa | LC_ALL=C sort | cut -f 2 | sed "1 s/^/${i}\n/" > temp/kraken2/${i}_count;done
# 提取第一样本品行名为表行名
header=`tail -n 1 metadata.txt | cut -f 1`
echo $header
tail -n+2 temp/kraken2/${header}.mpa | LC_ALL=C sort | cut -f 1 | \
sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
head -n3 temp/kraken2/0header_count
# paste合并样本为表格
ls temp/kraken2/*count
paste temp/kraken2/*count > 04.kraken2/tax_count.mpa
# 检查表格及统计
/mnt/project/tools/utilis/csvtk/csvtk -t stat 04.kraken2/tax_count.mpa
grep -E "Taxonomy|s_" taxonomy.tsv | sed 's/^.*s__//g' | cut -f1,3-33  > taxonomy.S.tsv
#################################################4.2 Bracken估计丰度########################################
循环重新估计每个样品的丰度 ,并对reads 矫正
# 设置估算的分类级别D,P,C,O,F,G,S，常用 P和S
#####################提取估算后reads  第六列#######################
tax=P
mkdir -p temp/bracken_${tax}
for i in `tail -n+2 metadata.txt|cut -f1`;do
# i=C1
/mnt/project/tools/Bracken/bracken -d /mnt/data/kraken2_mask \
-i temp/kraken2/${i}.report \
-r 125 -l ${tax} -t 0 \
-o temp/bracken_${tax}/${i};done

# 结果描述：共7列，分别为物种名、ID、分类级、读长计数、补充读长计数、**总数、百分比**
# name    taxonomy_id     taxonomy_lvl    kraken_assigned_reads   added_reads     new_est_reads        fraction_total_reads
# Capnocytophaga sputigena        1019    S       4798    996     5794    0.23041
# Capnocytophaga sp. oral taxon 878       1316596 S       239     21      260     0.01034

#bracken结果new_est_reads合并成表
for i in `tail -n+2 metadata.txt|cut -f1`;do tail -n+2 temp/bracken_${tax}/${i} | sort -k 1 | cut -f6 | sed "1 s/^/${i}\n/" > temp/bracken_${tax}/${i}.count;done
# 提取第一样本品行名为表行名
h=`tail -n1 metadata.txt|cut -f1`
tail -n+2 temp/bracken_${tax}/${h}|sort -k 1 |cut -f1 | sed "1 s/^/Taxonomy\n/" > temp/bracken_${tax}/0header.count
# 检查文件数，为n+1
ls temp/bracken_${tax}/*count | wc
# paste合并样本为表格，并删除非零行
paste temp/bracken_${tax}/*count > 04.kraken2/bracken.${tax}.count.txt
# 统计行列，默认去除表头
/mnt/project/tools/utilis/csvtk/csvtk -t stat 04.kraken2/bracken.${tax}.count.txt


# 设置估算的分类级别D,P,C,O,F,G,S，常用 P和S
######################提取相对丰度  第七列#######################
for i in `tail -n+2 metadata.txt|cut -f1`;do tail -n+2 temp/bracken_${tax}/${i} | sort -k 1 | cut -f7 | sed "1 s/^/${i}\n/" > temp/bracken_${tax}/${i}.abu;done
# 提取第一样本品行名为表行名
h=`tail -n1 metadata.txt|cut -f1`
tail -n+2 temp/bracken_${tax}/${h}| sort -k 1 | cut -f1 | sed "1 s/^/Taxonomy\n/" > temp/bracken_${tax}/0header.abu
# 检查文件数，为n+1
ls temp/bracken_${tax}/*abu | wc
# paste合并样本为表格，并删除非零行
paste temp/bracken_${tax}/*abu > 04.kraken2/bracken.${tax}.abu.txt
# 统计行列，默认去除表头
/mnt/project/tools/utilis/csvtk/csvtk -t stat 04.kraken2/bracken.${tax}.abu.txt
#awk -F '\t' '{if ($2 > 0.0 && $3 >0.0 ) {print $0}}' 04.kraken2/bracken.${tax}.abu.txt > 04.kraken2/bracken.${tax}.abu.gt0.txt

######heatmap######
Rscript ~/bracken.hclust_heatmap.R   -i 04.kraken2/bracken.G.count.txt   -t Genus -n 25   -w 183 -e 118   -o 04.kraken2/Heatmap


######################################################################
# humann_regroup_table -h 
# usage: humann_regroup_table [-h] [-i INPUT]
#                             [-g {uniref90_rxn,uniref50_rxn,uniref50_go,uniref90_go,uniref50_ko,uniref90_ko,uniref50_level4ec,uniref90_level4ec,uniref50_pfam,uniref90_pfam,uniref50_eggnog,uniref90_eggnog}]
#                             [-c CUSTOM] [-r] [-f {sum,mean}] [-e PRECISION]
#                             [-u {Y,N}] [-p {Y,N}] [-o OUTPUT]

# HUMAnN utility for regrouping table features
# =============================================
# Given a table of feature values and a mapping
# of groups to component features, produce a 
# new table with group values in place of 
# feature values.

# optional arguments:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Original output table (tsv or biom format); default=[TSV/STDIN]
#   -g {uniref90_rxn,uniref50_rxn,uniref50_go,uniref90_go,uniref50_ko,uniref90_ko,uniref50_level4ec,uniref90_level4ec,uniref50_pfam,uniref90_pfam,uniref50_eggnog,uniref90_eggnog}, --groups {uniref90_rxn,uniref50_rxn,uniref50_go,uniref90_go,uniref50_ko,uniref90_ko,uniref50_level4ec,uniref90_level4ec,uniref50_pfam,uniref90_pfam,uniref50_eggnog,uniref90_eggnog}
#                         Built-in grouping options
#   -c CUSTOM, --custom CUSTOM
#                         Custom groups file (.tsv or .tsv.gz format)
#   -r, --reversed        Custom groups file is reversed: mapping from features to groups
#   -f {sum,mean}, --function {sum,mean}
#                         How to combine grouped features; default=sum
#   -e PRECISION, --precision PRECISION
#                         Decimal places to round to after applying function; default=Don't round
#   -u {Y,N}, --ungrouped {Y,N}
#                         Include an 'UNGROUPED' group to capture features that did not belong to other groups? default=Y
#   -p {Y,N}, --protected {Y,N}
#                         Carry through protected features, such as 'UNMAPPED'? default=Y
#   -o OUTPUT, --output OUTPUT
#                         Path for modified output table; default=STDOUT

                        
#### 转换基因家族为KO(uniref90_ko)，可选eggNOG(uniref90_eggnog)或酶(uniref90_level4ec)
#ko
for i in `tail -n+2 metadata.txt|cut -f1`;do
humann_regroup_table \
-i temp/humann/${i}_genefamilies.tsv \
-g uniref90_ko \
-o temp/humann/${i}_ko.tsv
done


# 合并，并修正样本名
humann_join_tables \
--input temp/humann/ \
--file_name ko \
--output 03.humann/ko.tsv
sed -i '1s/_Abundance-RPKs//g' 03.humann/ko.tsv
tail 03.humann/ko.tsv

KO合并为高层次L2, L1通路代码
wc -l 03.humann/ko.tsv # 3797 lines
grep -v '|' 03.humann/ko.tsv > 03.humann/ko_clean.tsv
wc -l 03.humann/ko_clean.tsv


# 分层结果：功能和对应物种表(stratified)和功能组成表(unstratified)
humann_split_stratified_table \
  --input 03.humann/ko.tsv \
  --output 03.humann/ 

# KO to level 1/2/3, 也可切换至humann3或qiime2等Python3环境下运行
/mnt/home/huanggy/biosoft/EasyMicrobiome/script/summarizeAbundance.py \
  -i 03.humann/ko_unstratified.tsv \
  -m /mnt/home/huanggy/biosoft/EasyMicrobiome/kegg/KO1-4.txt \
  -c 2,3,4 -s ',+,+,' -n raw \
  -o 03.humann/KEGG


############################################可视化################################################
##metaphlan2 属水平/种水平  丰度热图
##Metaphlan2
sd=/mnt/home/huanggy/biosoft/EasyMicrobiome/script

Rscript $sd/metaphlan_hclust_heatmap.R \
  -i 02.metaphlan3/taxonomy.spf2 \
  -t Family -n 25 \
  -w 183 -e 118 \
  -o 02.metaphlan3/HeatmapFamily
  
  
$ Rscript $sd/metaphlan_hclust_heatmap.R \
>   -i 02.metaphlan3/taxonomy.spf2 \
>   -t Family -n 25 \
>   -w 183 -e 118 \
>   -o 02.metaphlan3/HeatmapFamily
[1] "The input file: 02.metaphlan3/taxonomy.spf2"
[1] "Taxonomy level: Family. Default if Species"
[1] "Number of taxonomy showing: 25"
[1] "Output heatmap filename: 02.metaphlan3/HeatmapFamily"
Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  line 287 did not have 39 elements   ###这一行taxid 也没有
Calls: read.table -> scan



###########################kraken##########################################
### Alpha多样性

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
Rscript $sd/kraken2alpha.R -h    
# -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.txt
Rscript $sd/kraken2alpha.R \
  --input 04.kraken2/tax_count.mpa \
  --depth 0 \
  --species 04.kraken2/tax_count.txt \
  --normalize 04.kraken2/tax_count.norm \
  --output 04.kraken2/tax_count.alpha
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
Rscript $sd/alpha_boxplot.R \
  -i 04.kraken2/tax_count.alpha \
  -a shannon \
  -d metadata \
  -n Cachexia \
  -o 04.kraken2/ \
  -w 89 -e 59

# 批量计算6种指数的箱线图+统计
for i in richness chao1 ACE shannon simpson invsimpson;do
Rscript $sd/alpha_boxplot.R -i 04.kraken2/tax_count.alpha -a ${i} \
  -d metadata -n Cachexia -w 89 -e 59 \
  -o 04.kraken2/
done


########################kraken_bracken#####################################
###heatmap 热图
Rscript ~/bracken.hclust_heatmap.R \
  -i 04.kraken2/bracken.G.count.txt \
  -t Genus -n 25 \
  -w 183 -e 118 \
  -o 04.kraken2/Heatmap

### 堆叠柱状图
	# 以门(P)/种(S)水平为例，结果包括output.sample/group.pdf两个文件
tax=P
Rscript ${sd}/tax_stackplot.R \
  --input 04.kraken2/bracken.${tax}.count.txt --design metadata \
  --group Cachexia --output 04.kraken2/bracken.${tax}.stackplot \
  --legend 10 --width 89 --height 59

### Alpha多样性
    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # 查看帮助
    Rscript $sd/otutab_rare.R -h    
    # -d指定最小样本量，默认0为最小值，抽平文件tax_norm.txt，alpha多样性tax_alpha.
    tax=S
    Rscript $sd/otutab_rare.R \
      --input 04.kraken2/bracken.${tax}.count.txt \
      --depth 0 --seed 1 \
      --normalize 04.kraken2/bracken.${tax}.norm \
      --output 04.kraken2/bracken.${tax}.alpha
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript $sd/alpha_boxplot.R -h # 查看参数
    mkdir -p 04.kraken2/${tax}
    Rscript $sd/alpha_boxplot.R \
      -i 04.kraken2/bracken.${tax}.alpha \
      -a shannon \
      -d metadata \
      -n Cachexia \
      -o 04.kraken2/${tax}/ \
      -w 89 -e 59
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript $sd/alpha_boxplot.R -i 04.kraken2/bracken.${tax}.alpha -a ${i} \
      -d metadata -n Cachexia -w 89 -e 59 \
      -o 04.kraken2/${tax}/
    done
    
### Beta多样性
    # Beta多样性距离矩阵计算
    # Mac 用不了面对的 usearch，使用在线平台 https://www.bic.ac.cn/BIC
    mkdir -p 04.kraken2/beta/
    /mnt/home/huanggy/biosoft/EasyMicrobiome/linux/usearch -beta_div 04.kraken2/bracken.${tax}.norm \
        -filename_prefix 04.kraken2/beta/

    # PCoA分析输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
    # 可选距离有 bray_curtis, euclidean, jaccard, manhattan
    dis=bray_curtis
    Rscript $sd/beta_pcoa.R \
      --input 04.kraken2/beta/${dis}.txt \
      --design metadata \
      --group Cachexia \
      --width 89 --height 59 \
      --output 04.kraken2/pcoa.${dis}.pdf