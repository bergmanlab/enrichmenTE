#!/usr/bin/env python3

import sys
import os
import time
from shutil import copyfile
import datetime
import logging
from enrichmenTE_utility import mkdir, get_lines, format_time
from enrichmenTE_input import get_args
from enrichmenTE_detect import detect
from enrichmenTE_cluster import cluster

# from enrichmenTE_cluster import cluster

"""
Author: Shunhua Han <shhan@uga.edu>
"""


def main():
    # Fetch command-line options
    args = get_args()

    # create output directory
    out_dir = args.out
    mkdir(out_dir)

    # logging config
    formatstr = "%(asctime)s: %(levelname)s: %(message)s"
    datestr = "%m/%d/%Y %H:%M:%S"
    # log_filename = ".".join(["enrichmenTE", prefix, "log"])
    now = datetime.datetime.now()
    log_filename = "enrichmenTE_{}.log".format(now.strftime("%d_%m_%Y_%T"))
    logging.basicConfig(
        level=logging.DEBUG,
        filename=os.path.join(out_dir, log_filename),
        filemode="w",
        format=formatstr,
        datefmt=datestr,
    )
    logging.info("CMD: " + " ".join(sys.argv))
    start_time = time.time()

    if not args.sub:
        print(
            "Please choose one of the two modes ('detect' or 'cluster'). See --help for more information."
        )
        return

    logging.info("WORKING DIR: {0}".format(os.path.abspath(out_dir)))
    logging.info("Start enrichmenTE...")

    # detect
    if args.sub == "detect":
        logging.info("MODE: DETECT")
        logging.info("Read 1: {0}".format(os.path.abspath(args.read1)))
        logging.info("Read 2: {0}".format(os.path.abspath(args.read2)))
        logging.info("Reference genome: {0}".format(os.path.abspath(args.reference)))

        # set up default parameters
        if args.thread is None:
            args.thread = 1

        if args.window is None:
            args.window = 100

        if args.tsd_max is None:
            args.tsd_max = 25

        if args.gap_max is None:
            args.gap_max = 25

        # get prefix
        if args.prefix is None:
            print(args.read1)
            r1 = (args.read1).split(",")
            prefix = os.path.basename(r1[0]).replace(".fastq.gz", "")
        else:
            prefix = args.prefix

        # create output directory
        out_dir = os.path.join(out_dir, prefix)
        mkdir(out_dir)

        # create module output directory
        detect_out_dir = os.path.join(out_dir, "detect")
        mkdir(detect_out_dir)

        # create directory for intermediate files
        tmp_dir = os.path.join(detect_out_dir, "detect_intermediate_files")
        mkdir(tmp_dir)

        nonref_pred = detect(
            prefix=prefix,
            read1=args.read1,
            read2=args.read2,
            reference=args.reference,
            outdir=tmp_dir,
            depth_config=args.depth_config,
            ref_te_bed=args.gff,
            window=args.window,
            tsd_max=args.tsd_max,
            gap_max=args.gap_max,
            filter_region=args.filter_region,
            thread=args.thread,
        )

        nonref_output = os.path.join(detect_out_dir, prefix + ".nonref.bed")
        copyfile(nonref_pred, nonref_output)
        num_nonref = get_lines(nonref_output)

        proc_time_all = time.time() - start_time
        logging.info("enrichmenTE DETECT finished in " + format_time(proc_time_all))
        logging.info("Number of non-reference TEs: " + str(num_nonref))
    elif args.sub == "cluster":
        logging.info("MODE: CLUSTER")

        # get prefix
        if args.prefix is None:
            prefix = "clustergram"
        else:
            prefix = args.prefix

        # create output directory
        out_dir = os.path.join(out_dir, prefix)
        mkdir(out_dir)

        # create module output directory
        cluster_out_dir = os.path.join(out_dir, "cluster")
        mkdir(cluster_out_dir)

        tree_file, pdf_file = cluster(
            prefix=prefix,
            enrichmente_out_dirs=args.enrichmente_out_dirs,
            filter_region=args.filter_region,
            outgroup=args.outgroup,
            out=cluster_out_dir,
            include_families=args.include_families,
            exclude_families=args.exclude_families,
            exclude_samples=args.exclude_samples,
        )

        proc_time_all = time.time() - start_time
        logging.info("enrichmenTE CLUSTER finished in " + format_time(proc_time_all))
        logging.info("Clustergram file in newick format: " + tree_file)
        logging.info("Clustergram file in PDF format: " + pdf_file)


main()
