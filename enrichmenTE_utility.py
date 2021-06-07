#!/usr/bin/env python3

import sys
import argparse
import os
import subprocess
from glob import glob
from datetime import datetime, timedelta
import re

"""
This script provides utility functions that are used by the main TE detection pipeline for multiplexed TE-NGS data
"""


def get_bam_stats(bam, stat, ref_index, families, dir):
    # this script generates alignment stats
    bed_unique = dir + "/" + "unique.bed"
    bed_none = dir + "/" + "none.bed"
    with open(ref_index, "r") as input, open(bed_unique, "w") as unique, open(
        bed_none, "w"
    ) as none:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            out_line = "\t".join([entry[0], "1", entry[1]])
            if "chr" in line:
                unique.write(out_line + "\n")
            elif all(x not in line for x in families):
                none.write(out_line + "\n")

    bam_name = os.path.basename(bam).replace(".bam", "")
    sample_name = bam_name.split(".")[0]
    read_name = bam_name.split(".")[1]
    with open(stat, "a") as output:
        output.write(sample_name + "\t")
        output.write(read_name + "\t")
        num = subprocess.Popen(["samtools", "view", "-c", bam], stdout=subprocess.PIPE)
        total_reads = num.stdout.read().decode().replace("\n", "")
        output.write(total_reads + "\t")
        num = subprocess.Popen(
            ["samtools", "view", "-c", "-F", "260", bam], stdout=subprocess.PIPE
        )
        mapped_reads = num.stdout.read().decode().replace("\n", "")
        pt = "{:.1%}".format(int(mapped_reads) / int(total_reads))
        output.write(pt + "\t")
        # unique region
        readc, pt = count_reads(bam, bed_unique, total_reads)
        output.write("{:.1%}".format(pt) + "\t")
        pt_tes = 0
        # per family
        for family in families:
            bed = dir + "/" + family + ".bed"
            with open(ref_index, "r") as INDEX, open(bed, "w") as BED:
                for line in INDEX:
                    entry = line.replace("\n", "").split("\t")
                    out_line = "\t".join([entry[0], "1", entry[1]])
                    if family in line:
                        BED.write(out_line + "\n")
            readc, pt = count_reads(bam, bed, total_reads)
            pt_tes = pt_tes + pt
            output.write("{:.1%}".format(pt) + "\t")
            os.remove(bed)
        # focal TE regions
        output.write("{:.1%}".format(pt_tes) + "\t")
        # none region
        readc, pt = count_reads(bam, bed_none, total_reads)
        output.write("{:.1%}".format(pt) + "\t")
        output.write("\n")
        os.remove(bed_none)
        os.remove(bed_unique)


def format_time(time):
    d = datetime(1, 1, 1) + timedelta(seconds=time)
    if d.hour == 0 and d.minute == 0:
        return "%d seconds" % (d.second)
    elif d.hour == 0 and d.minute != 0:
        return "%d minutes %d seconds" % (d.minute, d.second)
    else:
        return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)


def get_lines(path):
    if os.path.isfile(path) == False or os.stat(path).st_size == 0:
        count = 0
    else:
        count = len(open(path).readlines())
    return count


def bed_rm_overlap(bed_in, bed_out):
    bed_merge = bed_in + ".redundant"
    with open(bed_merge, "w") as output:
        command = 'bedtools merge -d 0 -o collapse -c 2,3,4,5,6 -delim "," -i ' + bed_in
        subprocess.call(command, shell=True, stdout=output)

    with open(bed_merge, "r") as input, open(bed_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if "," not in entry[3]:
                chromosome = entry[0]
                start = entry[3]
                end = entry[4]
                info = entry[5]
                score = entry[6]
                strand = entry[7]
                out_line = "\t".join([chromosome, start, end, info, score, strand])
                output.write(out_line + "\n")
    os.remove(bed_merge)


def merge_bed(bed_in, bed_out, genome, filter_method="overlap"):
    # merge bed files from all families, check overlap, merge or remove entries if necessary
    bed_merge = bed_out + ".merge.tmp"
    with open(bed_merge, "w") as output:
        for bed in bed_in:
            if os.path.isfile(bed) and os.stat(bed).st_size != 0:
                with open(bed, "r") as input:
                    for line in input:
                        output.write(line)

    if get_lines(bed_merge) != 0:
        # sort bed files
        bed_sort = bed_out + ".merge.sort.tmp"
        with open(bed_sort, "w") as output:
            subprocess.call(
                ["bedtools", "sort", "-i", bed_merge, "-g", genome], stdout=output
            )
        os.remove(bed_merge)

        # remove entries if overlap with multiple families
        if filter_method == "overlap":
            bed_rm_overlap(bed_sort, bed_out)
        else:
            bed_rm_duplicate(bed_sort, bed_out)
        os.remove(bed_sort)

    else:
        os.rename(bed_merge, bed_out)


def bed_rm_duplicate(bed_in, bed_out):
    with open(bed_out, "w") as output:
        command = "cat " + bed_in + " | sort | uniq"
        subprocess.call(command, shell=True, stdout=output)


def get_genome_file(ref):
    subprocess.call(["samtools", "faidx", ref])
    ref_index = ref + ".fai"
    genome = ref + ".genome"
    with open(ref_index, "r") as input, open(genome, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            out_line = "\t".join([entry[0], entry[1]])
            output.write(out_line + "\n")
    return genome


def get_cluster(bam, bed, config, window, set, family):
    # parse cutoff file
    cutoff = 30
    with open(config, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if family == entry[0] and set == entry[1]:
                cutoff = int(entry[2])

    # generate potential TE enriched clusters based on depth profile
    depth = bam + ".depth"
    with open(depth, "w") as output:
        subprocess.call(["samtools", "depth", bam, "-d", "0", "-Q", "1"], stdout=output)

    if os.path.isfile(depth) == False or os.stat(depth).st_size == 0:
        print("No depth info for " + family + "\n")
        return None

    depth_filter = depth + ".filter"
    with open(depth, "r") as input, open(depth_filter, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if int(entry[2]) > cutoff:
                out_line = "\t".join(
                    [entry[0], str(entry[1]), str(int(entry[1]) + 1), str(entry[2])]
                )
                output.write(out_line + "\n")

    if os.path.isfile(depth_filter) == False or os.stat(depth_filter).st_size == 0:
        print("No depth info for " + family + "\n")
        return None

    bed_tmp = bed + ".tmp"
    with open(bed_tmp, "w") as output:
        subprocess.call(
            [
                "bedtools",
                "merge",
                "-d",
                str(window),
                "-c",
                "4",
                "-o",
                "mean",
                "-i",
                depth_filter,
            ],
            stdout=output,
        )

    # filter out non-chr entries
    with open(bed, "w") as output, open(bed_tmp, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            contig = entry[0]
            if contig != family:
                output.write(line)
                # out_line = "\t".join(
                #     [entry[0], entry[1], str(int(entry[2]) + 1), entry[3]]
                # )
            # if "chr" in entry[0]:

    os.remove(depth)
    os.remove(depth_filter)
    os.remove(bed_tmp)
    return None


def count_reads(bam, bed, total_reads):
    # cout number of reads mapped to a specific region
    command = "bedtools intersect -abam " + bam + " -b " + bed + " | samtools view -c"
    num = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    unique_reads = num.stdout.read().decode().replace("\n", "")
    pt = int(unique_reads) / int(total_reads)
    return unique_reads, pt


def extract_reads(bam, family, set, out):
    # extract read list from R2 BAM file that belong to set1/set2
    if set == "set2":
        command = (
            "samtools view -F 0x10 " + bam + " " + family + " | cut -f1 | sort | uniq"
        )
    else:
        command = (
            "samtools view -f 0x10 " + bam + " " + family + " | cut -f1 | sort | uniq"
        )
    with open(out, "w") as output:
        subprocess.call(command, stdout=output, shell=True)


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        print("Creation of the directory %s failed" % dir)
    else:
        print("Successfully created the directory %s " % dir)


def make_bam(fq, ref, thread, bam):
    # alignment and generate sorted bam file
    sam = bam + ".sam"
    with open(sam, "w") as output:
        subprocess.call(["bwa", "mem", "-v", "0", "-t", thread, ref, fq], stdout=output)

    command = (
        "samtools view -Sb -t "
        + ref
        + " "
        + sam
        + " | "
        + "samtools sort -@ "
        + thread
        + " -o "
        + bam
    )
    subprocess.call(command, shell=True)
    subprocess.call(["samtools", "index", bam])
    os.remove(sam)


def filter_bam(bam_in, read_list, bam_out):
    # filter bam file based on provided read list
    command = (
        "picard FilterSamReads I="
        + bam_in
        + " O="
        + bam_out
        + " READ_LIST_FILE="
        + read_list
        + " FILTER=includeReadList"
    )
    with open(bam_out, "w") as output:
        subprocess.call(command, stdout=output, shell=True)
    subprocess.call(["samtools", "index", bam_out])


def get_nonref_te(bed1, bed2, norm, family, stats, sample_dir, tmp_dir, sample_name):
    overlap = tmp_dir + "/" + sample_name + ".overlap.tsv"
    if os.path.isfile(bed1) and os.path.isfile(bed2):
        # generate non-reference TE predictions based on cluster signals from set1 and set2 R1 alignments
        with open(overlap, "w") as output:
            subprocess.call(
                ["bedtools", "intersect", "-wao", "-a", bed1, "-b", bed2], stdout=output
            )

    # parse overlap
    if os.path.isfile(overlap):
        nonref = tmp_dir + "/" + sample_name + "." + family + ".nonref.bed"
        with open(overlap, "r") as input, open(nonref, "w") as output:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                if entry[4] != "." and int(entry[8]) <= 12:
                    chr = entry[0]
                    if abs(int(entry[1]) - int(entry[6])) < abs(
                        int(entry[2]) - int(entry[5])
                    ):
                        start = entry[1]
                        end = entry[6]
                    else:
                        start = entry[5]
                        end = entry[2]
                    score = (float(entry[3]) + float(entry[7])) / 2
                    score = "{:.2f}".format(score)
                    strand = "."
                    # family = "copia"
                    out_line = "\t".join(
                        [chr, str(start), str(end), family, str(score), strand]
                    )
                    output.write(out_line + "\n")
        nonref_count = get_lines(nonref)
        os.remove(overlap)

        if os.path.isfile(nonref):
            nonref_norm = (
                sample_dir + "/" + sample_name + "." + family + ".nonref.norm.bed"
            )
            with open(nonref_norm, "w") as output:
                subprocess.call(
                    ["bedtools", "intersect", "-a", nonref, "-b", norm, "-u"],
                    stdout=output,
                )
            nonref_norm_count = get_lines(nonref_norm)
        else:
            nonref_norm_count = 0
    else:
        nonref_count = 0
        nonref_norm_count = 0

    with open(stats, "a") as output:
        output.write(sample_name + "\t")
        output.write(family + "\t")
        output.write("non-ref\t")
        output.write(str(nonref_count) + "\t")
        output.write(str(nonref_norm_count) + "\n")


def get_nonref(bed1, bed2, outdir, family, tsd_max, gap_max):
    overlap = outdir + "/" + family + ".overlap.tsv"
    if os.path.isfile(bed1) and os.path.isfile(bed2):
        with open(overlap, "w") as output:
            # subprocess.call(
            #     ["bedtools", "intersect", "-wao", "-a", bed1, "-b", bed2], stdout=output
            # )
            subprocess.call(
                ["bedtools", "window", "-w", str(gap_max), "-a", bed1, "-b", bed2],
                stdout=output,
            )

    # parse overlap
    if os.path.isfile(overlap) and os.stat(overlap).st_size != 0:
        nonref = outdir + "/" + family + ".nonref.bed"
        with open(overlap, "r") as input, open(nonref, "w") as output:
            for line in input:
                ins_pass = True
                entry = line.replace("\n", "").split("\t")
                chromosome = entry[0]
                if (
                    (int(entry[1]) - int(entry[5])) > 0
                    and (int(entry[2]) - int(entry[6])) > 0
                ) or (
                    (int(entry[5]) - int(entry[1])) > 0
                    and (int(entry[6]) - int(entry[2])) > 0
                ):  # get rid of entries if one is within another
                    if abs(int(entry[1]) - int(entry[6])) < abs(
                        int(entry[2]) - int(entry[5])
                    ):
                        if int(entry[1]) < int(entry[6]):
                            start = entry[1]
                            end = entry[6]
                        else:
                            start = entry[6]
                            end = entry[1]
                        strand = "-"
                        dist = int(entry[1]) - int(
                            entry[6]
                        )  # positive: gap; negative: overlap
                    else:
                        if int(entry[2]) < int(entry[5]):
                            start = entry[2]
                            end = entry[5]
                        else:
                            start = entry[5]
                            end = entry[2]
                        strand = "+"
                        dist = int(entry[5]) - int(
                            entry[2]
                        )  # positive: gap; negative: overlap
                    if dist < 0:
                        if -dist > tsd_max:
                            ins_pass = False
                    else:
                        if dist > gap_max:
                            ins_pass = False

                    if ins_pass:
                        score = "."
                        family_info = "|".join([family, str(dist)])
                        out_line = "\t".join(
                            [
                                chromosome,
                                str(start),
                                str(end),
                                family,  # TODO: save dist
                                str(score),
                                strand,
                            ]
                        )
                        output.write(out_line + "\n")
        # os.remove(overlap)


def get_ref(bed1, bed2, rm_bed, out_dir, family, window=50):
    # calculate clusters that jointly support ref TEs (all, norm) with a percentage
    ref_rm = out_dir + "/" + family + ".ref_rm.bed"
    family_ref_count = 0
    with open(ref_rm, "w") as output, open(rm_bed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if entry[3] == family:
                family_ref_count = family_ref_count + 1
                output.write(line)

    if family_ref_count > 0:
        # calculate clusters that jointly support ref TEs (all, norm) with a percentage
        ref_sm = bed1.replace(".bed", ".ref.bed")
        if os.path.isfile(bed1):
            with open(ref_sm, "w") as output:
                subprocess.call(
                    [
                        "bedtools",
                        "window",
                        "-w",
                        str(window),
                        "-a",
                        ref_rm,
                        "-b",
                        bed1,
                        "-u",
                    ],
                    stdout=output,
                )
        ref_ms = bed2.replace(".bed", ".ref.bed")
        if os.path.isfile(bed2):
            with open(ref_ms, "w") as output:
                subprocess.call(
                    [
                        "bedtools",
                        "window",
                        "-w",
                        str(window),
                        "-a",
                        ref_rm,
                        "-b",
                        bed2,
                        "-u",
                    ],
                    stdout=output,
                )
        if os.path.isfile(ref_sm) and os.path.isfile(ref_ms):
            ref_both = out_dir + "/" + family + ".reference.bed"
            with open(ref_both, "w") as output:
                subprocess.call(
                    ["bedtools", "intersect", "-a", ref_sm, "-b", ref_ms, "-u"],
                    stdout=output,
                )
        # if os.path.isfile(ref_sm):
        #     os.remove(ref_sm)
        # if os.path.isfile(ref_ms):
        #     os.remove(ref_ms)


def get_ref_te(
    bed1, bed2, norm, gff, family, stat, sample_dir, tmp_dir, sample_name, window=100
):
    # get reference count
    ref_te = tmp_dir + "/" + family + ".ref.bed"
    with open(ref_te, "w") as output, open(gff, "r") as input:
        for line in input:
            check_family = "Name=" + family + "{}"
            if check_family in line:
                entry = line.replace("\n", "").split("\t")
                family = entry[8].split(";")[1]
                family = re.sub("Name=", "", family)
                family = re.sub("{}.*", "", family)
                out_line = "\t".join(
                    [entry[0], entry[3], entry[4], family, ".", entry[6]]
                )
                output.write(out_line + "\n")
    ref_count = get_lines(ref_te)
    ref_te_norm = tmp_dir + "/" + family + ".ref.norm.bed"
    with open(ref_te_norm, "w") as output:
        subprocess.call(
            ["bedtools", "intersect", "-a", ref_te, "-b", norm, "-u"], stdout=output
        )

    # calculate clusters that jointly support ref TEs (all, norm) with a percentage
    ref_set1 = bed1.replace(".bed", ".ref.bed")
    ref_set1_norm = bed1.replace(".bed", ".ref.norm.bed")
    if os.path.isfile(bed1):
        with open(ref_set1, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te,
                    "-b",
                    bed1,
                    "-u",
                ],
                stdout=output,
            )
        ref_set1_count = get_lines(ref_set1)

        with open(ref_set1_norm, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te_norm,
                    "-b",
                    bed1,
                    "-u",
                ],
                stdout=output,
            )
        ref_set1_norm_count = get_lines(ref_set1_norm)
    else:
        ref_set1_count = 0
        ref_set1_norm_count = 0

    ref_set2 = bed2.replace(".bed", ".ref.bed")
    ref_set2_norm = bed2.replace(".bed", ".ref.norm.bed")
    if os.path.isfile(bed2):
        with open(ref_set2, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te,
                    "-b",
                    bed2,
                    "-u",
                ],
                stdout=output,
            )
        ref_set2_count = get_lines(ref_set2)

        with open(ref_set2_norm, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te_norm,
                    "-b",
                    bed2,
                    "-u",
                ],
                stdout=output,
            )
        ref_set2_norm_count = get_lines(ref_set2_norm)
    else:
        ref_set2_count = 0
        ref_set2_norm_count = 0

    # joinly support ref
    if os.path.isfile(ref_set1) and os.path.isfile(ref_set2):
        ref_both = tmp_dir + "/" + sample_name + "." + family + ".ref.bed"
        with open(ref_both, "w") as output:
            subprocess.call(
                ["bedtools", "intersect", "-a", ref_set1, "-b", ref_set2, "-u"],
                stdout=output,
            )
        ref_both_count = get_lines(ref_both)
    else:
        ref_both_count = 0

    if os.path.isfile(ref_set1_norm) and os.path.isfile(ref_set2_norm):
        ref_both_norm = sample_dir + "/" + sample_name + "." + family + ".ref.norm.bed"
        with open(ref_both_norm, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "intersect",
                    "-a",
                    ref_set1_norm,
                    "-b",
                    ref_set2_norm,
                    "-u",
                ],
                stdout=output,
            )
        ref_both_norm_count = get_lines(ref_both_norm)
    else:
        ref_both_norm_count = 0

    # write stats to summary
    with open(stat, "a") as output:
        output.write(sample_name + "\t")
        output.write(family + "\t")
        output.write("ref\t")
        out_num = (
            str(ref_both_count)
            + "("
            + str(ref_set1_count)
            + ","
            + str(ref_set2_count)
            + ")"
        )
        output.write(out_num + "\t")
        out_num = (
            str(ref_both_norm_count)
            + "("
            + str(ref_set1_norm_count)
            + ","
            + str(ref_set2_norm_count)
            + ")"
        )
        output.write(out_num + "\n")

    # clean tmp files
    os.remove(ref_te)
    os.remove(ref_te_norm)
