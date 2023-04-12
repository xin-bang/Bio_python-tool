#!/usr/bin/python
# -*- coding: UTF-8 -*-
from Bio import SeqIO
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description=f"Get what you want in your genebank file,\
        For example: python {sys.argv[0]} -i input -o output -seek 'pro/cds/metadata'")
parser.add_argument("-i", type=str, metavar='',
                    help="genebank file")
parser.add_argument("-o", type=str, metavar='',
                    help="the output result")
parser.add_argument("-seek", type=str, metavar='',
                    help="Choose which part you want to seek. Here we offer three modes of seeking.\n"
                    "'pro': Get Protein sequence from Genebank.\n"
                    "'cds': Get CDS sequence from Genebank.\n"
                    "'metadata': Get meta-information from Genebank and output non Excel tables. \n"
                    "Meta-information include accession number,organism,cds length,pro length,strain,country,collection_date,host,Defintion\n"
                    "Of course,this information you can download XML file from NCBI database")
args = parser.parse_args()


# file1 = open(args.i,'r')
# file1 = open("sequence.gb","r")
file1 = "M_tuberculosis.gbff"
# file1 = args.i


def get_pro(file1):
    with open(file1, "r") as gb_file:
        for record in SeqIO.parse(gb_file, "genbank"):
            # SeqRecord 对象有一个 annotations 属性，它是一个 Python 字典，包含了与该记录相关的注释信息。这些注释信息可以是关于该序列来源、采集日期、文献引用、分类学信息等等
            # collection_date = record.features[0].qualifiers.get("date", None)
            for feature in record.features:
                if feature.type == "CDS":
                    # feature.qualifiers 是一个字典，它包含了一个 SeqFeature 对象的所有限定符信息，即所有的genebank中的feature框的信息
                    if "translation" in feature.qualifiers:
                        protein_seq = feature.qualifiers["translation"][0]
                        protein_id = feature.qualifiers['protein_id'][0]
                        product = feature.qualifiers["product"][0]
                        gene_id = str(feature.qualifiers.get('gene', None)).replace(
                            '\'', '').replace('[', '').replace(']', '')
                        result = (
                            f">{protein_id}_{record.id}_{gene_id}_({product})\n{protein_seq}\n")
                        with open(args.o, 'a') as f:  # 使用a代替w模式，这样这样每次写入的内容将会被追加到文件末尾，而不是覆盖之前的内容
                            f.write(result)



def get_cds (file1):
    with open(file1, "r") as gb_file:
        for record in SeqIO.parse(gb_file, "genbank"):
            sub_seq = 0  #这里设置sub_seq为全局变量，这样在第一个if语句中得到的sub_seq可以在全局使用
            for feature in record.features:
                if feature.type == "gene":
                    seq_start = feature.location.nofuzzy_start
                    seq_end = feature.location.nofuzzy_end
                    strand = feature.location.strand
                    if strand == 1:
                        sub_seq = record.seq[seq_start:seq_end]
                    else:
                        sub_seq = record.seq[seq_start:seq_end].reverse_complement()
            for feature in record.features:           
                if feature.type == "CDS":
                    protein_id = feature.qualifiers['protein_id'][0]
                    product = feature.qualifiers["product"][0]
                    gene_id = str(feature.qualifiers.get('gene', None)).replace(
                        '\'', '').replace('[', '').replace(']', '')                                
                    result = (
                        f">{protein_id}_{record.id}_{gene_id}_({product})\n{sub_seq}\n")
                    with open(args.o,'a') as f:
                        f.write(result)



def get_metadata(file1):
    with open(file1, "r") as gb_file:
        meta_data = {}
        for record in SeqIO.parse(gb_file, "genbank"):
            definition = record.description
            accession = record.id
            organism = record.annotations.get("organism", None)
            # cds_len = len(record.seq) ##添加每个cds以及蛋白序列的长度
            # OrderedDict（即qualifiers），对象是一个字典,可以使用get() 函数来获取改字典中某个键对应值。
            collection_date = str(record.features[0].qualifiers.get(
                "collection_date", None)).replace('\'', '').replace('[', '').replace(']', '')
            host = str(record.features[0].qualifiers.get("host", None)).replace(
                '\'', '').replace('[', '').replace(']', '')
            country = str(record.features[0].qualifiers.get("country", None)).replace(
                '\'', '').replace('[', '').replace(']', '')
            strain = str(record.features[0].qualifiers.get("strain", None)).replace(
                '\'', '').replace('[', '').replace(']', '')
            meta_data[accession] = {'collection_date': collection_date, 'host': host,
                                    'country': country, 'strain': strain, 'definition': {definition}, 'organism': {organism}}
        df = pd.DataFrame(meta_data).T
        df.index.name = "accession number"
        df.to_excel(args.o+".xlsx", index=True)


def main():
    if args.seek == "pro":
        get_pro(args.i)
    elif args.seek == "cds":
        get_cds(args.i)
    else:
        get_metadata(args.i)


if __name__ == "__main__":
    main()
