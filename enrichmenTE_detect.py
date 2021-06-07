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

    # sample_name = os.path.basename(read1_list[0]).replace(
    #     ".fastq.gz", ""
    # )  # TODO: generalize
    # # sample_name = re.sub("-set1.*", "-set1-and-set2", sample_name)
    # # sample_name = re.sub("-set3.*", "-set3", sample_name)

    # sample_dir = outdir + "/" + prefix
    # mkdir(sample_dir)

    # summary_dir = sample_dir + "/summary"
    # mkdir(summary_dir)

    # detection_summary = summary_dir + "/" + sample_name + ".detection.summary"
    # map_summary = outdir + "/" + prefix + ".alignment.summary"
    # demultiplex_summary = outdir + "/" + prefix + ".demultiplex.summary"
    # input_file_list = summary_dir + "/" + sample_name + ".input.txt"

    # # write header for summary file
    # with open(map_summary, "w") as output:
    #     out_line = "\t".join(
    #         [
    #             "sample_name",
    #             "read_name",
    #             "total",
    #             "mapped",
    #             "unique",
    #             "1731",
    #             "297",
    #             "copia",
    #             "mdg1",
    #             "mdg3",
    #             "roo",
    #             "focus_families",
    #             "other",
    #         ]
    #     )
    #     output.write(out_line + "\n")

    # with open(demultiplex_summary, "w") as output:
    #     out_line = "\t".join(
    #         [
    #             "sample_name",
    #             "demultiplex_family",
    #             "demultiplex_set",
    #             "demultiplex_count",
    #         ]
    #     )
    #     output.write(out_line + "\n")

    # with open(detection_summary, "w") as output:
    #     out_line = "\t".join(
    #         [
    #             "sample_name",
    #             "family",
    #             "te_type",
    #             "whole_genome (set1,set2)",
    #             "norm_region (set1,set2)",
    #         ]
    #     )
    #     output.write(out_line + "\n")

    # with open(input_file_list, "w") as output:
    #     output.write(input_files + "\n")

    # tmp_dir = sample_dir + "/intermediate"
    # mkdir(tmp_dir)

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

    # # report alignment stats for R1 and R2
    # ref_index = reference + ".fai"
    # get_bam_stats(bam_r1, map_summary, ref_index, families, outdir)
    # get_bam_stats(bam_r2, map_summary, ref_index, families, outdir)

    # demultiplex BAM file by family and primer set
    sets = ["set1", "set2"]
    for family in families:
        family_dir = os.path.join(outdir, family)
        mkdir(family_dir)
        for set in sets:
            group_name = ".".join([prefix, family, set])
            read_list = family_dir + "/" + group_name + ".txt"
            extract_reads(bam_r2, family, set, read_list)

            # # output number of reads by set and by family
            # with open(demultiplex_summary, "a") as output:
            #     output.write(prefix + "\t")
            #     output.write(family + "\t")
            #     output.write(set + "\t")
            #     read_count = get_lines(read_list)
            #     output.write(str(read_count) + "\n")

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
            # os.remove(bam_set1)
        if os.path.isfile(bam_set2):
            get_cluster(
                bam_set2,
                bed_set2,
                config=depth_config,
                window=window,
                set="set2",
                family=family,
            )
            # os.remove(bam_set2)

        # if (
        #     ref_te_bed is not None
        #     and os.path.isfile(bed_set1)
        #     and os.path.isfile(bed_set2)
        # ):
        #     get_ref(bed_set1, bed_set2, ref_te_bed, family_dir, family)

        if os.path.isfile(bed_set1) and os.path.isfile(bed_set2):
            get_nonref(bed_set1, bed_set2, family_dir, family, tsd_max, gap_max)

        # get_ref_te(
        #     bed_set1,
        #     bed_set2,
        #     filter_region,
        #     ref_te_bed,
        #     family,
        #     detection_summary,
        #     sample_dir,
        #     tmp_dir,
        #     sample_name,
        # )
        # get_nonref_te(
        #     bed_set1,
        #     bed_set2,
        #     filter_region,
        #     family,
        #     detection_summary,
        #     sample_dir,
        #     tmp_dir,
        #     sample_name,
        # )

    # gather non-reference TE predictions
    nonref_bed = outdir + "/" + prefix + ".nonref.bed"
    pattern = "/*/*.nonref.bed"
    bed_files = glob(outdir + pattern, recursive=True)
    genome = get_genome_file(reference)
    merge_bed(
        bed_in=bed_files, bed_out=nonref_bed, genome=genome, filter_method="overlap"
    )

    # # gather reference TE predictions
    # ref_bed = outdir + "/" + prefix + ".ref.bed"
    # if ref_te_bed is not None:
    #     pattern = "/*/*.reference.bed"
    #     bed_files = glob(outdir + pattern, recursive=True)
    #     merge_bed(
    #         bed_in=bed_files, bed_out=ref_bed, genome=genome, filter_method="duplicate"
    #     )
    # else:
    #     open(ref_bed, "w").close()

    # # merge bed files from all families and sort
    # ## merge non-reference TE bed files
    # merge_nonref_tmp = sample_dir + "/" + sample_name + ".merge.nonref.tmp"
    # with open(merge_nonref_tmp, "w") as output:
    #     for family in families:
    #         bed = sample_dir + "/" + sample_name + "." + family + ".nonref.norm.bed"
    #         if os.path.isfile(bed) and os.stat(bed).st_size != 0:
    #             with open(bed, "r") as input:
    #                 for line in input:
    #                     output.write(line + "\n")
    # merge_nonref = sample_dir + "/" + sample_name + ".merge.nonref.norm.bed"
    # with open(merge_nonref, "w") as output:
    #     subprocess.call(["bedtools", "sort", "-i", merge_nonref_tmp], stdout=output)
    # os.remove(merge_nonref_tmp)

    # ## merge reference TE bed files
    # merge_ref_tmp = sample_dir + "/" + sample_name + ".merge.ref.tmp"
    # with open(merge_ref_tmp, "w") as output:
    #     for family in families:
    #         bed = sample_dir + "/" + sample_name + "." + family + ".ref.norm.bed"
    #         if os.path.isfile(bed) and os.stat(bed).st_size != 0:
    #             with open(bed, "r") as input:
    #                 for line in input:
    #                     output.write(line + "\n")
    # merge_ref = sample_dir + "/" + sample_name + ".merge.ref.norm.bed"
    # with open(merge_ref, "w") as output:
    #     subprocess.call(["bedtools", "sort", "-i", merge_ref_tmp], stdout=output)
    # os.remove(merge_ref_tmp)

    return nonref_bed
