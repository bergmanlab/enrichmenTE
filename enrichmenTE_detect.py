#!/usr/bin/env python3

import sys
import argparse
import os
import subprocess
from glob import glob
import re
from enrichmenTE_utility import *

"""
this script takes family and primer set multiplexed pair end illumina data (if data is primer set/family non-multiplexed, a list of R1/R2 can be provided) and output reference/non-reference TEs for six focal TE families: 'copia', '1731', 'roo', 'mdg3', 'mdg1', '297'.
"""


def detect(
    prefix,
    read1,
    read2,
    reference,
    outdir,
    depth_config,
    ref_te_bed,
    window,
    tsd_max,
    gap_max,
    filter_region,
    thread,
):

    families = ["1731", "297", "copia", "mdg1", "mdg3", "roo"]

    # process input read list
    read1_list = read1.split(",")
    read2_list = read2.split(",")

    input_files = read1 + "," + read2
    input_files = input_files.replace(" ", "")

    # unzip and merge input files for R1 and R2
    fq1 = outdir + "/" + prefix + ".R1.fastq"
    with open(fq1, "w") as output:
        for read1 in read1_list:
            subprocess.call(["gunzip", "-c", read1], stdout=output)
    fq2 = outdir + "/" + prefix + ".R2.fastq"
    with open(fq2, "w") as output:
        for read2 in read2_list:
            subprocess.call(["gunzip", "-c", read2], stdout=output)

    # align R2 to masked augmented reference genome
    bam_r2 = outdir + "/" + prefix + ".R2.bam"
    make_bam(fq2, reference, str(thread), bam_r2)

    # align R1 to masked augmented reference genome
    bam_r1 = outdir + "/" + prefix + ".R1.bam"
    make_bam(fq1, reference, str(thread), bam_r1)

    # clean fastq files
    os.remove(fq1)
    os.remove(fq2)

    # demultiplex BAM file by family and primer set
    sets = ["set1", "set2"]
    for family in families:
        family_dir = os.path.join(outdir, family)
        mkdir(family_dir)
        for set in sets:
            group_name = ".".join([prefix, family, set])
            read_list = family_dir + "/" + group_name + ".txt"
            extract_reads(bam_r2, family, set, read_list)

            # use read ID to filter read 1 bam file by family and set
            if os.path.isfile(read_list) and os.stat(read_list).st_size != 0:
                bam_out = family_dir + "/" + group_name + ".R1.bam"
                filter_bam(bam_r1, read_list, bam_out)
            else:
                print("family with no R2 alignment for " + set + ": " + family)
            os.remove(read_list)

    # for each family, generate reference/non-reference TE predictions and report summary
    for family in families:
        family_dir = os.path.join(outdir, family)
        bam_set1 = family_dir + "/" + prefix + "." + family + ".set1.R1.bam"
        bam_set2 = family_dir + "/" + prefix + "." + family + ".set2.R1.bam"

        bed_set1 = bam_set1.replace("bam", "cluster.bed")
        bed_set2 = bam_set2.replace("bam", "cluster.bed")

        if os.path.isfile(bam_set1):
            get_cluster(
                bam_set1,
                bed_set1,
                config=depth_config,
                window=window,
                set="set1",
                family=family,
            )
            os.remove(bam_set1)
        if os.path.isfile(bam_set2):
            get_cluster(
                bam_set2,
                bed_set2,
                config=depth_config,
                window=window,
                set="set2",
                family=family,
            )
            os.remove(bam_set2)

        if os.path.isfile(bed_set1) and os.path.isfile(bed_set2):
            get_nonref(bed_set1, bed_set2, family_dir, family, tsd_max, gap_max)

    # gather non-reference TE predictions
    nonref_bed = outdir + "/" + prefix + ".nonref.bed"
    pattern = "/*/*.nonref.bed"
    bed_files = glob(outdir + pattern, recursive=True)
    genome = get_genome_file(reference)
    merge_bed(
        bed_in=bed_files, bed_out=nonref_bed, genome=genome, filter_method="overlap"
    )
    return nonref_bed
