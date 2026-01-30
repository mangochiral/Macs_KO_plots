#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:50:35 2026

@author: chandrima.modak
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 15:31:31 2025

@author: chandrima.modak
"""
import re
import sys
import os
import pandas as pd
import argparse
from  preprocess import *
import multiprocessing as mp
import preprocess_macs as pprm

def _run_one(run):
    """Worker for a single (lane, sample_name) job."""
    datadir, output_dir, exp, lane, sample_name = run
    try:
        path_sample_dir = os.path.join(
            datadir, sample_name, "processed/cellranger/outs/per_sample_outs", lane
        )

        xdata = os.path.join(path_sample_dir, "count", "sample_filtered_feature_bc_matrix.h5")

        if not os.path.exists(xdata):
            print(f"Missing input: {xdata}")
            return (lane, sample_name)

        gex_a, crispr_a, pre_filter_counts, post_filter_counts = pprm.process_cellranger_h5_macs(
            xdata, exp, sample_name, lane
        )

        sample_outdir = os.path.join(output_dir, f"{sample_name}_{lane}")
        os.makedirs(sample_outdir, exist_ok=True)

        gex_out = os.path.join(sample_outdir, f"{sample_name}_gex_preprocessed.h5ad")
        crispr_out = os.path.join(sample_outdir, f"{sample_name}_crispr_preprocessed.h5ad")

        gex_a.write_h5ad(gex_out)
        crispr_a.write_h5ad(crispr_out)

        with open(os.path.join(sample_outdir, "stats_report.txt"), "w", encoding="utf-8") as file:
            file.write(f"Pre filter cells count: {pre_filter_counts}\n")
            file.write(f"Post filter cells count: {post_filter_counts}\n")

        return (lane, sample_name)

    except Exception as e:
        print(f"Error for {sample_name},{lane}: {e}")
        return (lane, sample_name)

    

def main():
    parser = argparse.ArgumentParser(description='Processing CRISPR experiment')
    parser.add_argument('--datadir',type=str,required=True,help='Path to directory that contains CRISPRia_Cellanome_<lane> folders')
    parser.add_argument('--exp',type=str,default='crispr',help="Experiment type; for Perturb-seq keep as 'crispr'")
    parser.add_argument('--exp_info',type=str,default='experiment_info.csv',help="Experiment metadata csv")
    parser.add_argument('--output_dir',type=str,required=True,help='Path to directory that processed data to be saved')
    parser.add_argument('--nprocs', type=int, help='Number of worker processes')
    args = parser.parse_args()
    
   
    exp_info = pd.read_csv((os.path.join(args.datadir, args.exp_info)))
    # Make sure output_dir exists
    os.makedirs(args.output_dir, exist_ok=True)
    num_cores = min(args.nprocs, mp.cpu_count())
    
    # Running lanes as one job
    runs = [
    (args.datadir, args.output_dir, args.exp, lane, sample_name)
    for lane in exp_info.columns
    for sample_name in exp_info[lane].dropna().astype(str)
    if os.path.isdir(os.path.join(args.datadir, sample_name))
    ]
    ctx =  mp.get_context('fork' if sys.platform != 'win32' else 'spawn')
    with ctx.Pool(processes= num_cores)as pool:
        for lane, sample_name in pool.imap(_run_one, runs):
            print(f'completed {sample_name},{lane}')

if __name__ == "__main__":
    main()

