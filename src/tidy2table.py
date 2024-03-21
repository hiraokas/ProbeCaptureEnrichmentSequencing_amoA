#!/usr/bin/env python
#coding: UTF-8

import sys 
import re
import os
import pandas as pd

def _help():
    print("""
==================================================================================================================================================================================
DESCRIPTION:
    Satoshi Hiraoka
    hiraokas@jamstec.go.jp
    created: 20230904
    History: 20230904
    - This script is for converting tidy data to table using pandas
INPUT
    sample  OTU  count
USAGE:
    this.sh tidy.tsv
OUTPUT FORMAT:
    table.tsv
==================================================================================================================================================================================
 """)


if __name__=='__main__':
    
    argvs = sys.argv
    argc  = len(argvs)
    if argc<2 or argvs[1]=="-h":
        _help()
        exit()

    filename = argvs[1]

    print("Road {}".format(filename))
    df = pd.read_csv(filename, sep='\t')

    #df.groupby(['Sample', 'OTU']).mean()
    print("convert")
    df_table = pd.pivot_table(df, index='OTU', columns='sample', values='count', fill_value=0)

    #output
    outputfile=filename+".table"
    df_table.to_csv(outputfile, sep='\t', index=True)
    print("output: {}".format(outputfile))


    print("")
    print("All done.")