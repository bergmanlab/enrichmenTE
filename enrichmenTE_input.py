#!/usr/bin/env python3

import argparse
import os


def get_args():
    parser = argparse.ArgumentParser(
        description="""Script to detect reference and non-reference TEs from TE-NGS data
    
    enrichmenTE consists of two major steps:
    - DETECT detects TE insertions
    - CLUSTER add given non-reference TE profiles from samples to a fixed phylogeny
    """
    )
    subparsers = parser.add_subparsers(help="modes", dest="sub")

    parser_detect = subparsers.add_parser(
        "detect",
        help="Detect TE insertions using TE-NGS data",
    )

    ## required ##
    parser_detect.add_argument(
        "-1",
        "--read1",
        type=str,
        help="read 1 fastq.gz file or file list",
        required=True,
    )
    parser_detect.add_argument(
        "-2",
        "--read2",
        type=str,
        help="read 2 fastq.gz file or file list",
        required=True,
    )
    parser_detect.add_argument(
        "-r",
        "--reference",
        type=str,
        help="masked augmented reference genome",
        required=True,
    )
    parser_detect.add_argument(
        "-f",
        "--filter_region",
        type=str,
        help="filter out region in bed format",
        required=True,
    )
    parser_detect.add_argument(
        "-g",
        "--gff",
        type=str,
        help="reference TE annotation in bed file",
        required=True,
    )
    parser_detect.add_argument(
        "-d",
        "--depth_config",
        type=str,
        help="depth cutoff config file for identifying TE clusters",
        required=True,
    )
    ## optional
    parser_detect.add_argument(
        "-p",
        "--prefix",
        help="Prefix for output files",
        required=False,
    )
    parser_detect.add_argument(
        "-w",
        "--window",
        type=int,
        help="merge window for identifying TE clusters (default = 100bp) ",
        required=False,
    )
    parser_detect.add_argument(
        "--tsd_max", type=int, help="maximum TSD size (default = 25) ", required=False
    )
    parser_detect.add_argument(
        "--gap_max", type=int, help="maximum gap size (default = 25) ", required=False
    )
    parser_detect.add_argument(
        "-t", "--thread", type=int, help="thread (default = 1) ", required=False
    )
    parser_detect.add_argument(
        "-o", "--out", type=str, help="output dir (default = '.') ", required=False
    )

    # parser_detect = subparsers.add_parser(
    #     "cluster",
    #     help="Detect TE insertions using TE-NGS data",
    # )

    args = parser.parse_args()

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.thread is None:
        args.thread = 1

    if args.window is None:
        args.window = 100

    if args.tsd_max is None:
        args.tsd_max = 25

    if args.gap_max is None:
        args.gap_max = 25

    return args
