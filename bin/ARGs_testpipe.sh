#!/bin/bash
  
sample=$1
outdir=$2
input_single=$3
if [ $# -lt 3 ];then
    echo "Usage :sh $0 <sample> <outdir> <fq> "
    exit 1
fi

fastp=/mnt/db/software/fastp/fastp
bowtie2_bin=/mnt/project/tools/bowtie2-2.4.4/bowtie2
salmon_card_index=/mnt/home/huanggy/CARD/localDB/card.salmon
human_database_index=/mnt/db/kneaddata/human_genome/hg37dec_v0.1


mkdir -p ${outdir}/result/${sample}/filter
mkdir -p ${outdir}/result/${sample}/amr

##############################filter low quality & replicate############################################
${fastp} -i ${input_single} --dedup  --length_required 50 -q 15 -u 40 --thread 8 --json  ${outdir}/result/${sample}/filter/${sample}.json --html  ${outdir}/result/${sample}/filter/${sample}.html -o ${outdir}/result/${sample}/filter/${sample}.clean.fq.gz
echo $?
echo 'filter finished'

##############################remove host####################################################
${bowtie2_bin} --mm --very-sensitive \
    -x ${human_database_index} \
    -U ${outdir}/result/${sample}/filter/${sample}.clean.fq.gz \
    --threads 16 \
    --un ${outdir}/result/${sample}/filter/${sample}.filter.fq \
    --no-unal --no-head -k 1 -S ${outdir}/result/${sample}/filter/${sample}.host.sam \
    1> ${outdir}/result/${sample}/filter/${sample}.bowtie.log 2>&1
rm -f ${outdir}/result/${sample}/filter/${sample}.host.sam
echo $?
echo 'hostremove finished'
####################################基于耐药基因比对############################################
/mnt/project/tools/seqtk/seqtk seq -a ${outdir}/result/${sample}/filter/${sample}.filter.fq > ${outdir}/result/${sample}/filter/${sample}.filter.fa
/mnt/project/tools/ncbi-blast-2.12.0+/bin/blastn -query ${outdir}/result/${sample}/filter/${sample}.filter.fa -out ${outdir}/result/${sample}/filter/${sample}.blastm6.xls -outfmt 6 -num_threads 20  -evalue 1e-10 -db ~/CARD/blastDB/card.nucl
python  ./blastm6_stat.py -b  ${outdir}/result/${sample}/filter/${sample}.blastm6.xls  -p ${sample}
python  ./blastm6_anno.py  ../DB/card_3.2.3.gff3  ${outdir}/result/${sample}/filter/${sample}.blastm6.stat.xls
