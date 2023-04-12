# Bio_python-tool

使用python去处理日常的生物信息学问题：

主要依赖的python包：Biopython；Panda

## Fre_cal.py 用于计算氨基酸或核苷酸频率

主要功能：可以计算多重比对序列后每个位点或每条序列的氨基酸以及核苷酸频率，并以excel表格的形式返回结果。

```python
usage: Freq_cal.py [-h] [-i] [-o] [-model] [-type] [-ignore_gap]

Calculate the amino acid or nucleotide frequency of multiple alignment sequences, For example: python Freq_cal.py -i input.fasta
-o output -model 'row/col' -type 'nuc/pro' -ignore_gap 'Ture/False'

options:
  -h, --help    show this help message and exit
  -i            must be mulit-alignment fasta
  -o            the output result of amino acid or nucleotide frequency
  -model        Choose which calculation mode. you can calculate the pro(nuc) frequency of each sequence or site by setting
                '-model' to 'row' or 'col',Default: row
  -type         Choose which sequence type.'nuc'(nucleotide) or 'pro'(protein),Default:'pro'
  -ignore_gap   whether or not ignore gap in multiple alignments,Default = 'False'
```

* Example:

```python
#计算多重比对的每个位点的氨基酸频率
python Freq_cal.py -i example_data/H9N2_HA_gene_pro.fas -o example_data/H9N2_HA_gene_pro_col -model col -type pro -ignore_gap False

#计算多重比对的每条序列的氨基酸频率
python Freq_cal.py -i example_data/H9N2_HA_gene_pro.fas -o example_data/H9N2_HA_gene_pro_row -model row -type pro -ignore_gap False

```



## 用于Genebank文件中获取相应的序列及信息

主要功能：用于从NCBI数据库下载的GenebanK格式的文件中索引CDS序列、Protein序列、Metadata信息等。

```python
usage: gbk_handle.py [-h] [-i] [-o] [-seek]

Get what you want in your genebank file, For example: python gbk_handle.py -i input -o output -seek 'pro/cds/metadata'

options:
  -h, --help  show this help message and exit
  -i          genebank file
  -o          the output result
  -seek       Choose which part you want to seek. Here we offer three modes of seeking. 'pro': Get Protein sequence from Genebank.    
              'cds': Get CDS sequence from Genebank. 'metadata': Get meta-information from Genebank and output non Excel tables. Meta-  
              information include accession number,organism,cds length,pro length,strain,country,collection_date,host,Defintion Of    
              course,this information you can download XML file from NCBI database
```

* Eaxmple:

  ```python
  #从Genebank文件中获得CDS序列：
  python gbk_handle.py -i example_data/M_tuberculosis.gbff -o example_data/M_tuberculosis.cds.fasta -seek cds

  #从Genebank文件中获得protein序列：
  python gbk_handle.py -i example_data/M_tuberculosis.gbff -o example_data/M_tuberculosis.pro.fasta -seek pro

  #从Genebank文件中获得相应的metadata信息：
  python gbk_handle.py -i example_data/H9N2_HA_gene.gb -o example_data/H9N2_HA_gene.metadata -seek metadata
  ```
