
import pandas as pd
from collections import defaultdict
import argparse,os

blast_out = "/mnt/home/huanggy/project/20220620_amr/result/21JS944075/filter/21JS944075.blastm6.xls"
def parser_length(x):
    end = int(x.split("_")[-1])
    start = int(x.split("_")[-2])
    return  end - start + 1

def seq_id_counts(seri):
    counts  = {}
    #assert isinstance(seri,pd.Series)
    for i in seri.to_list():
        if i not in counts:
            counts[i] = 1
        else:
            counts[i] +=1
    return counts

# def parser_subseqid_breadth(_df):
#     subseqid_interver = defaultdict(list)
#     subseqid_breadth = defaultdict(int)
#     for ind, row in _df.iterrows():
#         sseq_id = row['subseqid']
#         sstart = row["sstart"]
#         send = row["send"]
#         subseqid_interver[sseq_id].append((sstart, send) if send > sstart else (send, sstart))
#
#     for i in subseqid_interver.items():
#         sid = i[0]
#         interver_list = i[1]
#         if isinstance(interver_list, list):
#             subseqid_breadth[sid] = 0
#             for (start, end) in interver_list:
#                 subseqid_breadth[sid] += end - start
#
#     return subseqid_interver, subseqid_breadth

def merge(subseqid_interver_dic):
    m = {}
    for i in subseqid_interver_dic:
        intervals = subseqid_interver_dic[i]
        intervals.sort(key=lambda x: x[0])
        merged = []
        for interval in intervals:
            # 如果列表为空，或者当前区间与上一区间不重合，直接添加
            if not merged or merged[-1][1] < interval[0]:
                li = list(interval)
                merged.append(li)
            else:
                # 否则的话，我们就可以与上一区间进行合并
                merged[-1][1] = max(merged[-1][1], interval[1])
        m.update({i:merged})
    return m

def getLength(subseqid_interver_merged_dic):
    newlen = {}
    for i in subseqid_interver_merged_dic:
        length = 0
        for intervals in subseqid_interver_merged_dic[i]:
            length += (intervals[1] - intervals[0] + 1)
        newlen.update({i:length})
    return newlen

def parser_subseqid_breadth_2(_df):
    subseqid_interver = defaultdict(list)
    for ind, row in _df.iterrows():
        sseq_id = row['subseqid']
        sstart = row["sstart"]
        send = row["send"]
        subseqid_interver[sseq_id].append((sstart, send) if send > sstart else (send, sstart))

    m = merge(subseqid_interver)
    subseqid_breadth = getLength(m)

    return subseqid_interver, subseqid_breadth



def parser_subseqid_len(_df):
    subseqid_len = defaultdict(int)
    for ind, row in _df.iterrows():
        sseq_id = row['subseqid']

        a = int(sseq_id.split("_")[-1]) - int(sseq_id.split("_")[-2]) + 1
        subseqid_len[sseq_id] = a
    return subseqid_len




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--blast',type= str,help='blast out m6 format')
    parser.add_argument('-p','--prefix',type = str,help='')
    args = parser.parse_args()


    _cols = ["qseqid", "subseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
             "evalue", "bitscore"]
    df = pd.read_csv(args.blast, sep="\t", header=None, names=_cols)

    identitymatch = 95
    evaluematch = 1e-10
    df = df[(df["pident"] >= identitymatch) & (df["evalue"] <= evaluematch)]

    subseqid_interver, subseqid_breadth = parser_subseqid_breadth_2(df)
    subseqid_len = parser_subseqid_len(df)

    df["cov"] = df["subseqid"].map(pd.Series(subseqid_breadth) / pd.Series(subseqid_len))
    df["cov"] = df["cov"].apply(lambda x:format(x,".3%"))
    # df3 = df.groupby("subseqid").agg({"qseqid":['count']})
    df["num_reads"] = df["subseqid"].map(seq_id_counts(df["subseqid"]))
    df2 = df.drop_duplicates(subset='subseqid', keep="first")


    out = os.path.join(os.path.dirname(args.blast) + "/" +  args.prefix + ".blastm6.stat.xls")
    df2.to_csv(out,sep='\t',index= False,header=True)


if __name__ == '__main__':
    main()

