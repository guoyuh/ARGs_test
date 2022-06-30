
from collections import defaultdict
import pandas as pd
import os,sys

def run(card_gff,blast_out):
    model_d = defaultdict(str)
    gene_d = defaultdict(str)
    taxonomy_name_d = defaultdict(str)
    tb = pd.read_table(card_gff,header=0,sep = "\t")
    for ind ,row in tb.iterrows():
        seqid = row['seqid']
        model_type = row['model_type']
        gene = row['attributes']
        taxonomy_name = row['taxonomy_name']

        model_d[seqid] = model_type
        gene_d[seqid] = gene
        taxonomy_name_d[seqid] = taxonomy_name

    tb_blast = pd.read_table(blast_out,header=0,sep = "\t")
    tb_blast["model_type"] = tb_blast["subseqid"].map(model_d)
    tb_blast["gene"] = tb_blast["subseqid"].map(gene_d)
    tb_blast["taxonomy_name"] = tb_blast["subseqid"].map(taxonomy_name_d)
    #df = tb_blast.sort_values(["num_reads", "evalue", "cov", ], ascending=[False, True, False])
    df = tb_blast.sort_values(["num_reads","cov", ], ascending=[False,False])
    df = df[["subseqid","gene","taxonomy_name","model_type","pident","evalue","bitscore","cov","num_reads"]]
    out = os.path.splitext(blast_out)[0] + ".anno.xls"
    df.to_csv(out,sep='\t',index=False,header=True)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <card_gff> <blastm6_stat_out>")
    else:
        run(card_gff=sys.argv[1],blast_out=sys.argv[2])