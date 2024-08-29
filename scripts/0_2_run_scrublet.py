#!/usr/bin/env python3
import sys,re,io
from collections import defaultdict
import os
import pandas as pd
import scrublet as scr
import scipy.io
import numpy as np

def main():
    counts, expected_doublets = getArgs()

    sampArr = counts.split("/")

    sampArrLen = len(sampArr) - 1
    sampDir ='/'.join(sampArr[:(sampArrLen)])
    sampName = counts.split('/')[-4]
    #sampName = sampArr[0] + '_' + sampArr[1]
    #print(sampDir)
    #print(sampName)

    counts = scipy.io.mmread(counts).T.tocsc()
    exp_doublets = getExpectedDoublets(expected_doublets = expected_doublets, sampName = sampName)
    scrub = scr.Scrublet(counts, expected_doublet_rate = exp_doublets, sim_doublet_ratio = 2)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85, n_prin_comps = 30)
    scrub.call_doublets(threshold = 0.2)
    outDf = pd.DataFrame({
        'scrublet_score_obs': scrub.doublet_scores_obs_,
        'scrublet_doublet': scrub.predicted_doublets_
    })
    simDf = pd.DataFrame({
        'scrublet_score_sim': scrub.doublet_scores_sim_
    })
    ### change here for filt/unfilt data
    outFile = sampDir + '/' + sampName + "_scrublet_res.tsv"
    bars = pd.read_table("barcodes.tsv", header = None)
    bars = bars.rename(columns = {0:"barcode"})
    bars['barcode'] = sampName + "_" + bars['barcode'].astype(str)
    outDf['barcode'] = bars['barcode']
    outDf.to_csv(outFile, sep = "\t", index = False)
    simDf.to_csv((sampDir + '/' + sampName + "_scrublet_sim.tsv"), sep = "\t", index = False)


def getArgs():
    import argparse
    parser = argparse.ArgumentParser(description = "run scrublet on a sparse matrix")
    parser.add_argument('-counts', action = 'store', type = str, help = 'sparse matrix containing counts data')
    parser.add_argument('-expected_doublets', action = 'store', type = str, help = 'table containing expected doublet frequencies')
    args = parser.parse_args()
    counts, expected_doublets = args.counts, args.expected_doublets
    return(counts, expected_doublets)

def getExpectedDoublets(expected_doublets, sampName):
    dubs = pd.read_table(expected_doublets, delimiter = "\t")
    expDubs = dubs[dubs['sample'] == sampName]['doublets'].iloc[0]
    return(expDubs)

if __name__ == "__main__":
	main()
