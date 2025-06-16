"""
Code: cfMEX_tools (Version 0.4.2)
By: Songphon Sutthitthasakul (Moon)
Date: 15-06-2025
"""

import argparse
import os
from multiprocessing import Pool
from time import time

from __init__ import *

base_dir = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(
        description="Extract multi-feature from cfDNA Fragment file."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Fragment size
    parser_FS = subparsers.add_parser("fragsize", help="Extract Fragment size")
    parser_FS.add_argument(
        "--input", "-i", type=str, help="Fragment file path", required=True
    )
    parser_FS.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output file name",
        required=False,
        default="FS_output",
    )
    parser_FS.add_argument("--id", type=str, help="Sample ID", required=False)

    # End motif
    parser_EM = subparsers.add_parser("endmotif", help="Extract End motif")
    parser_EM.add_argument(
        "--input", "-i", type=str, help="Fragment file path", required=True
    )
    parser_EM.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output file name",
        required=False,
        default="EM_output",
    )
    parser_EM.add_argument("--id", type=str, help="Sample ID", required=False)
    parser_EM.add_argument(
        "--location", "-l", type=str, help="Spcific location file", required=False, default=f"{base_dir}/dataset/EM_Window.csv"
    )
    parser_EM.add_argument(
        "--thread",
        "-t",
        type=int,
        help="Number of Threads (default: 1)",
        required=False,
        default=1,
    )

    # SpL Ratio
    parser_SR = subparsers.add_parser("splratio", help="Extract SpL ratio")
    parser_SR.add_argument(
        "--input", "-i", type=str, help="Fragment file path", required=True, default=f"{base_dir}/dataset/1Mb_Window_Num.csv"
    )
    parser_SR.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output file name",
        required=False,
        default="SR_output",
    )
    parser_SR.add_argument("--id", type=str, help="Sample ID", required=False)
    parser_SR.add_argument(
        "--location", "-l", type=str, help="Spcific location file", required=False
    )
    parser_SR.add_argument(
        "--thread",
        "-t",
        type=int,
        help="Number of Threads (default: 1)",
        required=False,
        default=1,
    )

    # BAM-2-Fragment
    parser_BF = subparsers.add_parser("bam2frag", help="Convert BAM to Fragment")
    parser_BF.add_argument(
        "--input", "-i", type=str, help="BAM file path", required=True
    )

    # Split Fragment
    parser_SF = subparsers.add_parser("splitfrag", help="Split Fragment by Range")
    parser_SF.add_argument(
        "--input", "-i", type=str, help="Fragment file path", required=True
    )
    parser_SF.add_argument(
        "--range", "-r", type=str, help="Spcific fragment range", required=False, default=f"{base_dir}/dataset/FS_Range.csv"
    )
    parser_SF.add_argument(
        "--thread",
        "-t",
        type=int,
        help="Number of Threads (default: 1)",
        required=False,
        default=1,
    )

    args = parser.parse_args()

    if args.command == "fragsize":
        start_time = time()

        c_format = check_n_Columns(args.input)
        FS_Run(args.input, args.id, args.output, c_format)
        Time_Stamp(start_time)

    elif args.command == "endmotif":
        start_time = time()

        if args.thread == 0:
            args.thread = 1

        if args.thread > 24:
            args.thread = 24

        if args.thread == 1:
            loc_n, location = singleThreads_Location(EM_window)
            c_format = check_n_Columns(args.input)
            EM_Run(location, args.input, args.id, c_format)
            cfMex_EM().save_csv(args.id, args.output, c_format)
            os.remove(f"{args.id}_Metadata_ref.csv")
            if c_format == "c6":
                os.remove(f"{args.id}_Metadata_ori.csv")
            Time_Stamp(start_time)

        elif args.thread > 1:
            loc_n, location = multiThreads_Location(EM_window, args.thread)
            c_format = check_n_Columns(args.input)
            with Pool(args.thread) as p:
                p.starmap(
                    EM_Run,
                    zip(
                        location,
                        [args.input] * args.thread,
                        [args.id] * args.thread,
                        [c_format] * args.thread,
                    ),
                )
            cfMex_EM().save_csv(args.id, args.output, c_format)
            os.remove(f"{args.id}_Metadata_ref.csv")
            if c_format == "c6":
                os.remove(f"{args.id}_Metadata_ori.csv")
            Time_Stamp(start_time)

    elif args.command == "splratio":
        start_time = time()

        if args.thread == 0:
            args.thread = 1

        if args.thread == 1:
            loc_n, location = singleThreads_Location(Num_1Mb_window)
            c_format = check_n_Columns(args.input)
            SR_Run(location, args.input, args.id, c_format)
            cfMex_SR().save_csv(args.id, args.output, loc_n)
            os.remove(f"{args.id}_Metadata.csv")
            Time_Stamp(start_time)

        elif args.thread > 1:
            loc_n, location = multiThreads_Location(Num_1Mb_window, args.thread)
            c_format = check_n_Columns(args.input)
            with Pool(args.thread) as p:
                p.starmap(
                    SR_Run,
                    zip(
                        location,
                        [args.input] * args.thread,
                        [args.id] * args.thread,
                        [c_format] * args.thread,
                    ),
                )
            cfMex_SR().save_csv(args.id, args.output, loc_n)
            os.remove(f"{args.id}_Metadata.csv")
            Time_Stamp(start_time)

    elif args.command == "bam2frag":
        start_time = time()
        BAM2FRAG(args.input).Convert()
        Time_Stamp(start_time)

    elif args.command == "splitfrag":
        start_time = time()

        if args.thread == 0:
            args.thread = 1

        if args.thread == 1:
            loc_n, location = singleThreads_Range(FS_range)
            c_format = check_n_Columns(args.input)
            Range_Run(location, args.input, c_format)
            Time_Stamp(start_time)

        elif args.thread > 1:
            loc_n, location = multiThreads_Range(FS_range, args.thread)
            c_format = check_n_Columns(args.input)
            with Pool(args.thread) as p:
                p.starmap(
                    Range_Run,
                    zip(location, [args.input] * args.thread, [c_format] * args.thread),
                )
            Time_Stamp(start_time)


if __name__ == "__main__":
    main()
