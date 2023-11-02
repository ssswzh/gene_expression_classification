#!/usr/bin/python3
# -*- coding:utf-8 -*-
'''Calculat one sample and other samples' aggregate's correlation.'''

import os
import sys
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter

def main():
    '''Parameter setting'''
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=RawTextHelpFormatter)
    parser.add_argument('-e','--expression',dest='expression',help='RNA gene expression,rows for gene and columns for sample',type=str,required=True)
    parser.add_argument('-s','--samlist',dest='samlist',help='sample list',type=str,required=True)
    parser.add_argument('-m','--method',dest='method',help='correlation method',choices=['pearson','spearman'],type=str,required=True)
    parser.add_argument('-l2on','--log2ornot',dest='log2ornot',help='whether proceed log2(data+1) or not',type=str,required=True)
    parser.add_argument('-am','--aggregatemethod',dest='aggregatemethod',help='aggregate method,such as sum ,mean and other',choices=['sum','mean','median'],type=str,required=True)
    parser.add_argument('-o','--odir',dest='odir',help='output path',type=str,required=True)
    parser.add_argument('-p','--prefix',dest='prefix',help='output prefix',type=str,required=True)
    argv = vars(parser.parse_args())

    expression = os.path.abspath(argv['expression'].strip())
    samlist = os.path.abspath(argv['samlist'].strip())
    method = argv['method']
    log2ornot = argv['log2ornot']
    aggregatemethod = argv['aggregatemethod']
    odir = os.path.abspath(argv['odir'].strip())
    prefix = argv['prefix'].strip()

    expression_df = pd.read_table(expression,header=0,index_col=0,sep='\t')
    sam_list = pd.read_table(samlist,header=None)[0].tolist()
    sam_list = [ str(sam) for sam in sam_list ]
    sam_aggregate_dict = {}
    for sam in sam_list:
        sam_expression_df = expression_df[[sam]]
        othersam_list = sam_list.copy()
        othersam_list.remove(sam)
        othersam_expression_df = expression_df[othersam_list]
        if aggregatemethod == 'sum':
            othersam_aggregate_expression_df = othersam_expression_df.apply(lambda x:x.sum(),axis=1)
        elif aggregatemethod == 'mean':
            othersam_aggregate_expression_df = othersam_expression_df.apply(lambda x:x.mean(),axis=1)
        elif aggregatemethod == 'median':
            othersam_aggregate_expression_df = othersam_expression_df.apply(lambda x:x.median(),axis=1)
        othersam_aggregate_expression_df = pd.DataFrame(othersam_aggregate_expression_df,columns=['othersam'])
        sam_othersam_expression_df = pd.concat([sam_expression_df,othersam_aggregate_expression_df],axis=1)
        if log2ornot == 'YES':
            sam_othersam_expression_df = np.log2(sam_othersam_expression_df+1)
        cor = sam_othersam_expression_df.corr(method = method).loc[sam,'othersam']
        sam_aggregate_dict[sam] = [cor]

    am_aggregate_df = pd.DataFrame(sam_aggregate_dict).T
    am_aggregate_df.columns = ['R']
    am_aggregate_df.index.name = 'Sample'
    am_aggregate_df.to_csv('{0}/{1}.leavel-one-out.cor.tsv'.format(odir,prefix),header=True,index=True,sep='\t')

if __name__ == '__main__':
    main()
