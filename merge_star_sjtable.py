#!/usr/bin/env python3
# coding: utf-8
# author: zhang.siwen
# data: 2023.03.14


import argparse
import os
import sys
import pandas as pd
import warnings
warnings.filterwarnings('ignore')


def GetArgs():
    parser = argparse.ArgumentParser(description='Merge STAR.SJ.out.tab files to a matrix', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--files', help='input STAR.SJ.out.tab files', action='store', dest='files', required=True, nargs='+')
    required.add_argument('--out', help='output file name with suffix and path', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--reads', help='count for uniquely-mapped reads OR multi-mapped reads OR both, DEFAULT "%(default)s"', default='unique', choices=['unique', 'multi', 'both'], action='store', dest='reads', required=False)
    optional.add_argument('--feature', help='row-wise sum( feature>=X ) >= Y*samples, keep feature,\nuse "none" for not filtering, \nformat: X,Y, separated by comma ",", DEFAULT "%(default)s"', default='2,0.2', action='store', dest='feature', required=False)
    args = parser.parse_args()
    return args


sjtable_colnames = ['chrom', 'start', 'end', 'strand', 'intron_motif', 'annotation', 'uniquely_mapping_reads', 'multi_mapping_reads', 'max_overhang']


def ReadSJTable(file, reads='unique'):
    sjtable = pd.read_csv(file, sep='\t', low_memory=False, index_col=False, names=sjtable_colnames)
    sjtable.index = sjtable[['chrom', 'start', 'end']].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
    if reads == 'unique':
        sjtable[['reads']] = sjtable[['uniquely_mapping_reads']]
    elif reads == 'multi':
        sjtable[['reads']] = sjtable[['multi_mapping_reads']]
    elif reads == 'both':
        sjtable[['reads']] = sjtable[['uniquely_mapping_reads']] + sjtable[['multi_mapping_reads']]
    return sjtable[['reads']]


def ReadFiles(files, reads='unique'):
    file_names = []
    file_dir = []
    file_path = []
    for f in files:
        f = os.path.abspath(f)
        file_names.append(os.path.basename(f).split('.')[0])
        file_dir.append(os.path.dirname(f).split('/')[-1])
        file_path.append(f)
    # build a dict for file id : file path
    if len(set(file_dir)) == len(file_path):
        file_result = dict(zip(file_dir, file_path))
    elif len(set(file_names)) == len(file_path):
        file_result = dict(zip(file_names, file_path))
    else:
        sys.exit("None of file folder name or file name is unique.")
    # read tables
    table_content = {}
    for name,file in file_result.items():
        table = ReadSJTable(file)
        table.columns = [name]
        table_content[name] = table
    return table_content


def MergeSJTables(table_content, out="output", feature='none'):
    tables = [i for i in table_content.values()]
    # merge all read count from all samples, fill NAs as 0, remove sites that have no counts in all samples
    tables_merged = pd.concat(tables, axis=1, join='outer').fillna(0)
    tables_merged = tables_merged[ tables_merged.sum(axis=1) != 0 ]
    if feature != 'none':
        feature_count, feature_ratio = feature.strip().split(',')
        keep_features = tables_merged[ tables_merged.ge(int(feature_count)).sum(axis=1) >= float(feature_ratio)*tables_merged.shape[1] ].index.to_list()
        tables_merged = tables_merged.loc[keep_features]
    # output
    outdir = os.path.dirname(os.path.abspath(out))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    tables_merged.index.name = 'name'
    tables_merged = tables_merged.astype(int)
    tables_merged.to_csv(out, sep='\t', index=True, mode='w')
    return None


def Main():
    args = GetArgs()
    table_content = ReadFiles(args.files, args.reads)
    MergeSJTables(table_content, args.out, args.feature)


if __name__ == '__main__':
    Main()
