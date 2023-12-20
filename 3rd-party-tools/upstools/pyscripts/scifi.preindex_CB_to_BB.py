#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='preindex+CB to BB. Preindex are at end of readname')
parser.add_argument('--in', type=str, dest="inbam", help='input bam file')

args = parser.parse_args()
inbamf = args.inbam

import os
import numpy as np
import pandas as pd
import pysam
from time import perf_counter as pc

def run():
    """ Run standard NMF on rank """
    start_time = pc()
    """ init input files """
    print("Add barcodes to field CU...")
    generate_bams(inbamf)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bams(inbamf):
    inb = pysam.AlignmentFile(inbamf, "rb")
    outb = pysam.AlignmentFile(inbamf + ".BB.bam", "wb", template=inb)
    for read in inb.fetch(until_eof=True):
        try:
            sp = read.get_tag("CB", with_value_type=False)
            #sp2 = read.get_tag("UB", with_value_type=False)
        except KeyError:
            continue
        if len(sp) != 18:
            print("Barcodes length in CB is not 18bp. 10X output is 16bp + '-1'.\n") ### since I will modifiy some files as my old habits..
        sp_bc = sp[0:16]
        pindex = read.qname.split(":")[-1]
        cu = pindex + sp_bc
        read.tags = read.tags + [("BB",cu)]
        
        outb.write(read)
    inb.close()
    outb.close()
    # os.system("mv " + inbamf + ".tmp " + str(inbamf))

if __name__ == "__main__":
    run()


