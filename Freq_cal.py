#!/usr/bin/python
from Bio import AlignIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
import pandas as pd
import sys
from pandas import ExcelWriter
import argparse


parser = argparse.ArgumentParser(
    description=f"Calculate the amino acid or nucleotide frequency of multiple alignment sequences,\
        For example: python {sys.argv[0]} -i input.fasta -o output -model 'row/col' -type 'nuc/pro' -ignore_gap 'Ture/False'")
parser.add_argument("-i", type=str, metavar='',
                    help="must be mulit-alignment fasta")
parser.add_argument("-o", type=str, metavar='',
                    help="the output result of amino acid or nucleotide frequency")
parser.add_argument("-model", type=str,metavar='',
                    help="Choose which calculation mode. you can calculate the pro(nuc) frequency of each sequence or site by setting '-model' to 'row' or 'col',Default: row")
parser.add_argument("-type", type=str,metavar='',
                    help="Choose which sequence type.'nuc'(nucleotide) or 'pro'(protein),Default:'pro'")
parser.add_argument("-ignore_gap",type=str,metavar='',
                    help="whether or not ignore gap in multiple alignments,Default = 'False'")
args = parser.parse_args()


# 读取多序列比对结果文件
alignment = AlignIO.read(args.i, "fasta")
# alignment = AlignIO.read("MDV_vaccine_all_pp38.align.fasta", "fasta")
# args.type = "nuc"
# args.model = "col"
# args.ignore_gap = "True"

def calculation_aa_freq(alignment, model,ignore_gap):
    '氨基酸频率计算'
    # 获取序列长度和序列个数
    seq_len = alignment.get_alignment_length()-1  # -1是为了排除最后一列全部的终止密码子的情况。
    seq_count = len(alignment)

    # 初始化频率矩阵
    freq_matrix = {}

    if model == "col":
        # 遍历每个位点，计算氨基酸频率
        if ignore_gap == "True":
            for i in range(seq_len):
                column = alignment[:, i].replace('-','')
                aa_freq = ProteinAnalysis(str(column)).get_amino_acids_percent()
                freq_matrix[i+1] = aa_freq
        else:
            for i in range(seq_len):
                column = alignment[:, i]
                aa_freq = ProteinAnalysis(str(column)).get_amino_acids_percent()
                freq_matrix[i+1] = aa_freq           
    else:
        if ignore_gap =="True":
            for i in range(seq_count):
                row = alignment[i, :].seq.replace('*', '').replace('-','')
                aa_freq = ProteinAnalysis(str(row)).get_amino_acids_percent()
                freq_matrix[alignment[i].id] = aa_freq
        else:
            for i in range(seq_count):
                row = alignment[i, :].seq.replace('*', '')
                aa_freq = ProteinAnalysis(str(row)).get_amino_acids_percent()
                freq_matrix[alignment[i].id] = aa_freq

    return(freq_matrix)


def calculation_nuc_fre(alignment, model,ignore_gap):
    '核苷酸频率的计算'
    seq_len = alignment.get_alignment_length()
    seq_count = len(alignment)
    freq_matrix = {}

    if model == "col":
        nuc_freq = {}
        if ignore_gap == "True":
            for i in range(seq_len):
                column = alignment[:, i].replace('-','')
                seq_no_gap_count = len(column)               
                for j in ["a", "t", "c", "g"]:
                    nuc_freq[j] = Seq(str(column)).lower().count(j) / seq_no_gap_count
                freq_matrix[i+1] = {'A': nuc_freq['a'], 'T': nuc_freq['t'],
                                    'C': nuc_freq['c'], 'G': nuc_freq['g']}
        else:
            for i in range(seq_len):
                column = alignment[:, i]
                for j in ["a", "t", "c", "g"]:
                    nuc_freq[j] = Seq(str(column)).lower().count(j) / seq_count
                freq_matrix[i+1] = {'A': nuc_freq['a'], 'T': nuc_freq['t'],
                                    'C': nuc_freq['c'], 'G': nuc_freq['g']}
    else:
        nuc_freq = {}
        if ignore_gap == "True":
            for i in range(seq_count):
                row = alignment[i, :].seq.replace('-','')
                seq_no_gap_count = len(row)
                for j in ["a", "t", "c", "g"]:
                    nuc_freq[j] = Seq(str(row)).lower().count(j) / seq_no_gap_count
                freq_matrix[i+1] = {'A': nuc_freq['a'], 'T': nuc_freq['t'],
                                    'C': nuc_freq['c'], "G": nuc_freq['g']}
        else:
            for i in range(seq_count):
                row = alignment[i, :].seq
                for j in ["a", "t", "c", "g"]:
                    nuc_freq[j] = Seq(str(row)).lower().count(j) / seq_len
                freq_matrix[i+1] = {'A': nuc_freq['a'], 'T': nuc_freq['t'],
                                    'C': nuc_freq['c'], "G": nuc_freq['g']}            
    return(freq_matrix)



def main():
    if args.type == "nuc":
        freq_matrix = calculation_nuc_fre(alignment, model=args.model,ignore_gap=args.ignore_gap)
    else:
        freq_matrix = calculation_aa_freq(alignment, model=args.model,ignore_gap=args.ignore_gap)

    # 将频率矩阵转换为数据框
    # 将字典freq_matrix转换为Pandas DataFrame格式。其中，字典的键会被转换为DataFrame的列名，字典的值会被转换为DataFrame的一列数据。所以后续需要转置。
    df = pd.DataFrame(freq_matrix).T
    df.index.name = "position"  # df.index是Pandas DataFrame对象的行索引
    df.to_excel(args.o+".xlsx", index=True)


if __name__ == "__main__":
    main()
