"""
Code: cfMEX_tools (Version 0.4.0)
By: Songphon Sutthitthasakul (Moon)
Date: 30-05-2025
"""

from __init__ import *
import os
import argparse
from time import time
from multiprocessing import Pool

base_dir = os.path.dirname(os.path.abspath(__file__))

def main():
    parser = argparse.ArgumentParser(description="Extract multi-feature from cfDNA Fragment file.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Fragment size
    parser_FS = subparsers.add_parser("fragsize", help="Extract End motif")
    parser_FS.add_argument("--input", '-i', type=str, help="Fragment file path", required=True)
    parser_FS.add_argument("--output", '-o', type=str, help="Output file name", required=False, default='FS_output')
    parser_FS.add_argument("--id", type=str, help="Sample ID", required=False)

    # End motif
    parser_EM = subparsers.add_parser("endmotif", help="Extract End motif")
    parser_EM.add_argument("--input", '-i', type=str, help="Fragment file path", required=True)
    parser_EM.add_argument("--output", '-o', type=str, help="Output file name", required=False, default='EM_output')
    parser_EM.add_argument("--id", type=str, help="Sample ID", required=False)
    parser_EM.add_argument("--thread", '-t', type=int, help="Number of Threads (default: 1)", required=False, default=1)

    # SpL Ratio
    parser_SR = subparsers.add_parser("splratio", help="Extract SpL ratio")
    parser_SR.add_argument("--input", '-i', type=str, help="Fragment file path", required=True)
    parser_SR.add_argument("--output", '-o', type=str, help="Output file name", required=False, default='SR_output')
    parser_SR.add_argument("--id", type=str, help="Sample ID", required=False)
    parser_SR.add_argument("--thread", '-t', type=int, help="Number of Threads (default: 1)", required=False, default=1)

    # BAM-2-Fragment
    parser_BF = subparsers.add_parser("bam2frag", help="Convert BAM to Fragment")
    parser_BF.add_argument("--input", '-i', type=str, help="BAM file path", required=True)
    
    args = parser.parse_args()

    if args.command == "fragsize":
        start_time = time()

        FS_Run(args.input, args.id, args.output)
        Time_Stamp(start_time)
        
    elif args.command == "endmotif":
        start_time = time()
        
        if args.thread == 0:
            args.thread = 1

        if args.thread > 24:
            args.thread = 24
            
        if args.thread == 1:
            location = singleThreads_Location(f"{base_dir}/{EM_window}")
            EM_Run(location, args.input, args.id)
            cfMex_EM().save_csv(args.id, args.output)
            os.remove(f'{args.id}_Metadata.csv')
            Time_Stamp(start_time)

        elif args.thread > 1:
            location = multiThreads_Location(f"{base_dir}/{EM_window}", args.thread)
            with Pool(args.thread) as p:
                p.starmap(EM_Run, zip(location, [args.input]*args.thread, [args.id]*args.thread))
            cfMex_EM().save_csv(args.id, args.output)
            os.remove(f'{args.id}_Metadata.csv')
            Time_Stamp(start_time)

    elif args.command == "splratio":
        start_time = time()
        
        if args.thread == 0:
            args.thread = 1
            
        if args.thread == 1:
            location = singleThreads_Location(f"{base_dir}/{Num_1Mb_window}")
            EM_Run(location, args.input, args.id)
            cfMex_SR().save_csv(args.id, args.output)
            os.remove(f'{args.id}_Metadata.csv')
            Time_Stamp(start_time)

        elif args.thread > 1:
            location = multiThreads_Location(f"{base_dir}/{Num_1Mb_window}", args.thread)
            with Pool(args.thread) as p:
                p.starmap(SR_Run, zip(location, [args.input]*args.thread, [args.id]*args.thread))
            cfMex_SR().save_csv(args.id, args.output)
            os.remove(f'{args.id}_Metadata.csv')
            Time_Stamp(start_time)
    
    elif args.command == "bam2frag":
        start_time = time()
        BAM2FRAG(args.input).Convert()
        Time_Stamp(start_time)
        

if __name__ == "__main__":
    main()
