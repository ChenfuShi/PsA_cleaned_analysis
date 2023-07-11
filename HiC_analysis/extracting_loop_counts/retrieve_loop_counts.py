########################################
# extra loops counts from a list of loops

# requires /mnt/jw01-aruk-home01/projects/functional_genomics/bin/new_basic_software

# the way this works is it will actually take the central pixel but also the pixel before and after. 
# so if the loop is 2.5kb and res is 2.5kb you will get a 3x3 matrix which then gets summed to 1 number
# if you have a 5k loop it would give you a 4x4 matrix. if it's aligned it should work fine

########################################

import os
import argparse
import pandas as pd
import numpy as np
import hicstraw 
from multiprocessing import Pool
from functools import partial
import math
def myfloor(x, base=2500):
    return base * math.floor(x/base)
def myceil(x, base=2500):
    return base * math.ceil(x/base)

def get_loop(x, mzd, res):
    A_start = myfloor(x["A_start"], res)
    A_end = myceil(x["A_end"], res)
    B_start = myfloor(x["B_start"], res)
    B_end = myceil(x["B_end"], res)
    
    return mzd.getRecordsAsMatrix(A_start-res, A_end, B_start-res, B_end).sum()

def retrieve_loops_chr(chrom, loops, hic_file,resolution):
    hic = hicstraw.HiCFile(hic_file)
    mzd = hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", resolution)
    loops_chr = loops[loops["chrA"] == chrom].copy()
    __get_loop = partial(get_loop, mzd=mzd, res=resolution)
    loops_chr["extracted_counts"] = loops_chr.apply(__get_loop, axis = 1)
    return loops_chr
    
def retrieve_loops_from_list(hic_file, f_loops, resolution = 2500):
    loops = pd.read_csv(f_loops, sep = "\t", dtype = {"BIN1_CHR":str,"BIN2_CHROMOSOME":str})
    loops.columns = ["chrA","A_start","A_end","chrB","B_start","B_end"] + list(loops.columns[6:])
    __retrieve_loops_chr = partial(retrieve_loops_chr, loops=loops,hic_file=hic_file,resolution=resolution)
    with Pool(4) as p:
        res = p.map(__retrieve_loops_chr, ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"])
    return pd.concat(res)


def main():
    parser = argparse.ArgumentParser(description='extract counts from Hi-C map given list of loops')

    parser.add_argument("-i",'--input', dest='in_loops', action='store', required=True,
                        help='input loops file')
    parser.add_argument("-o",'--output', dest='out_loops', action='store', required=True,
                        help='output loops file with extra column')
    parser.add_argument("-f",'--hic_file', dest='hic_file', action='store', required=True,
                        help='hic file to extra loops from')
    parser.add_argument("-r",'--resolution', dest='res', action='store', required=False, default = 2500,
                        help='resolution of the map to use')

    args = parser.parse_args()

    output_loops = retrieve_loops_from_list(args.hic_file, args.in_loops, int(args.res))
    output_loops.to_csv(args.out_loops, header = True, index = False, sep="\t")

if __name__=="__main__":
    main()