

###质控,宿主过滤 参考有参部分

######################所有样本一起组装###############################
megahit -t 6 \
        -1 `tail -n+2 metadata.txt|cut -f1 | sed 's/^/filter\//;s/$/.clean.1.gz/' |tr '\n' ','|sed 's/,$//'` \
        -2 `tail -n+2 metadata.txt|cut -f1 | sed 's/^/filter\//;s/$/.clean.2.gz/' |tr '\n' ','|sed 's/,$//'` \
        -o temp/megahit


##sed 's/^/'  匹配开头
######################每个样本单独组装################################
tail -n+2 metadata.txt|cut -f1 | sort | uniq |/mnt/db/tools/rush -j 96 \
"/mnt/db/tools/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -t 16 \
-1 filter/{}.clean.1.gz \
-2 filter/{}.clean.2.gz \
-o temp/{} "



#######################组装后过滤 ，再进行基因预测#################################################################################
tail -n+2 metadata.txt|cut -f1 | sort | uniq  | while read id ;do /mnt/project/tools/seqkit/seqkit seq -m 500 temp/${id}/final.contigs.fa > temp/${id}/${id}.contigs500.fa;done


projdir=/mnt/home/huanggy/project/ganju
sampleID=BE_1
mkdir -p 3.GenePrediction/${sampleID}
/mnt/project/tools/MetaGeneMark/mgm/gmhmmp \
            -a -d -f G \
            -m /mnt/project/tools/MetaGeneMark/mgm/MetaGeneMark_v1.mod \
            -A ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.fa \
            -D ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.fa \
            -o ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.contigs.fa.gff \
            temp/${sampleID}/${sampleID}.contigs500.fa


#对预测结果名称进行修改
sed "/>/s/gene/${sampleID}_gene/" ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.fa \
    | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \
    > ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.rename.fa


#####################################MetaGeneMark循环 预测 ##选择它不选择prodigal 是因为MetaGeneMark 快#############################
for sampleID in `tail -n+2 metadata.txt|cut -f1 | sort | uniq`;
do
            mkdir -p 3.GenePrediction/${sampleID};
            /mnt/project/tools/MetaGeneMark/mgm/gmhmmp \
            -a -d -f G \
            -m /mnt/project/tools/MetaGeneMark/mgm/MetaGeneMark_v1.mod \
            -A ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.fa \
            -D ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.fa \
            -o ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.contigs.fa.gff \
            temp/${sampleID}/${sampleID}.contigs500.fa;
      
done

### 对nucleotide 和protei 进行改名和删除空白行
for sampleID in `tail -n+2 metadata.txt|cut -f1 | sort | uniq`;
do
            ###核酸
            sed "/>/s/gene/${sampleID}_gene/" ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.fa \
            | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \
            > ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.rename.fa;    
            sed '/^\s*$/d' ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.rename.fa \
            > ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.delblank.fa;
            ##蛋白
            sed "/>/s/gene/${sampleID}_gene/" ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.fa \
            | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \
            > ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.rename.fa;    
            sed '/^\s*$/d' ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.rename.fa \
            > ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.delblank.fa;
done


########################################################合并#################################################################
mkdir -p 3.GenePrediction/Cluster
for sampleID in `tail -n+2 metadata.txt|cut -f1 | sort | uniq`;
do
cat ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.nucleotide.delblank.fa >> ${projdir}/3.GenePrediction/Cluster/total.nucleotide.fa;
cat ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.protein.delblank.fa >> ${projdir}/3.GenePrediction/Cluster/total.protein.fa;     
done

##############################################cd-hit去冗余 构建非冗余基因集#################################################
/mnt/project/tools/cd-hit-v4.8.1/cd-hit-est \
    -c 0.95 \
    -aS 0.9 \
    -G 0 \
    -g 1 \
    -d 0 \
    -M 6000 \
    -i ${projdir}/3.GenePrediction/Cluster/total.nucleotide.fa \
    -o ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa 

###发现基因数量200 万以上，很多基因的长度还是小于100 需要过滤
/mnt/project/tools/seqkit/seqkit seq -m 500 ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \
> ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide500.fa
# 翻译成蛋白序列并构建索引
/mnt/project/tools/seqkit/seqkit translate \
${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide500.fa \
> ${projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.fa 
/mnt/project/tools/bowtie2-2.4.4/bowtie2-build -f ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide500.fa \
${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide



##bowtie2 比对
# tail -n+2 metadata.txt|cut -f1 | sort | uniq |/mnt/project/tools/rush -j 96 \
# echo "/mnt/project/tools/bowtie2-2.4.4/bowtie2 -p 16 \
#             -x ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide \
#             -1 ${projdir}/filter/{}.R1.fq.gz \
#             -2 ${projdir}/filter/{}.R2.fq.gz \
#             --end-to-end --sensitive -I 100 -X 200 \
#             -S ${projdir}/3.GenePrediction/{}/{}.abundance.sam"
tail -n+2 metadata.txt|cut -f1 | sort | uniq |/mnt/project/tools/rush -j 96 \
"/mnt/project/tools/bowtie2-2.4.4/bowtie2 -p 16 \
            -x ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide \
            -1 ${projdir}/filter/{}.R1.fq.gz \
            -2 ${projdir}/filter/{}.R2.fq.gz \
            --end-to-end --sensitive -I 100 -X 200 \
            -S ${projdir}/3.GenePrediction/{}/{}.abundance.sam"
# gene支持数
for sampleID in `tail -n+2 metadata.txt|cut -f1 | sort | uniq`;
do
grep -v '^@'  ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.abundance.sam  \
| cut -f3 | awk '{{if($1 != "*")print $0}}'  \
| sort | uniq -c |  awk 'BEGIN{{FS=" "; OFS=","}}{{print $2, $1}}' \
| awk 'BEGIN{{FS=",";OFS=","}}{{if($2>2)print $1, $2;else print $1,"0"}}' \
> ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.gene.count.csv;
done

##基因长度
/mnt/project/tools/seqkit/seqkit fx2tab -j 30 -l  -n -i -H ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide500.fa \
> ${projdir}/3.GenePrediction/Cluster/gene.length.csv

##丰度计算
for sampleID in `tail -n+2 metadata.txt|cut -f1 | sort | uniq`;
do
python get_abudance_table.py --gene_count ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.gene.count.csv \
--gene_length ${projdir}/3.GenePrediction/Cluster/gene.length.csv \
--abundance_table  ${projdir}/3.GenePrediction/${sampleID}/${sampleID}.abundance.table
done

##合并丰度表与count 表
python merge_abudance_table.py --project_dir /mnt/home/huanggy/project/ganju --merge_abundance merge_abundance.txt --merge_count merge_count.txt
##TPM 或者fpkm ,merge count 表转化为 merge tpm 或者 FPKM
python exp_units_convert.py -exp_matrix /mnt/home/huanggy/project/ganju/merge_count.txt -gene_length /mnt/home/huanggy/project/ganju/3.GenePrediction/Cluster/gene.length.csv

##建库
/mnt/home/huanggy/biosoft/Ublastx_stageone2.3/bin/diamondv0.9.24 makedb --threads 96 --in /mnt/data/NCBI/NR/fasta/nr -d /mnt/data/nr/nr

#diamond NR数据库比对; 输出格式为daa, 向后兼容megan
mkdir -p 3.GenePrediction/anno
/mnt/home/huanggy/biosoft/Ublastx_stageone2.3/bin/diamondv0.9.24 \
	blastx \
	-d /mnt/data/nr/nr \
	-q ${projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide500.fa \
	-o ${projdir}/3.GenePrediction/anno/NonRundant.total.nucleotide.daa \
	--outfmt 100 \
	--evalue 0.001

## megan 安装
conda create -n megan megan
数据库下载
nohup wget -c --no-check-certificate https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db.zip &
unzip megan-map-Feb2022.db.zip
conda activate megan

mkdir -p 4.GeneAnnotation/anno
daa2rma \
-i ${projdir}/3.GenePrediction/anno/NonRundant.total.nucleotide.daa \
-t 48 \
-ms 50 \
-me 0.01 \
-top 50 \
-mdb /mnt/data/megan-map-Feb2022.db \
-o ${projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.rma &&\
rma2info \
 -i ${projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.rma \
 -r2c Taxonomy -n true -p true -v \
 > ${projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy

#################################################物种注释及结果整理###################################################################
python3 handle_taxonomy_file.py ${projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy ${projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy.simple

##基因2物种,样本2物种
python ./get_taxonomy_gene_abundantable.py  --abundance_table merge_abundance.txt --taxonomy_table ${projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy.simple \
--out ${projdir}/4.GeneAnnotation/table


###############################################
基因2功能??????