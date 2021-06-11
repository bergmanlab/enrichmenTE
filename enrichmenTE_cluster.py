#!/usr/bin/env python3
import sys
import argparse
import os
import subprocess
import glob
import time
import pandas as pd
from multiprocessing import Process, Pool
from Bio import Phylo

# import dendropy


# def get_args():
#     parser = argparse.ArgumentParser(
#         description="Script to generate clustering"
#     )

#     # required
#     parser.add_argument(
#         "--prefix",
#         type=str,
#         help="output prefix",
#         required=True,
#     )
#     parser.add_argument(
#         "--enrichmente_out_dirs",
#         type=str,
#         help="list of enrichmenTE output directories",
#         nargs="+",
#         required=True,
#     )
#     parser.add_argument(
#         "--filter_region",
#         type=str,
#         help="filter region in bed format, separated by comma",
#         required=True,
#     )
#     parser.add_argument("--outgroup", type=str, help="outgroup", required=False)
#     parser.add_argument(
#         "-o",
#         "--out",
#         type=str,
#         help="directory to output results)",
#         required=True,
#     )
#     parser.add_argument(
#         "--include_families",
#         type=str,
#         help="TE families to use in the phylogeny (separated by comma)",
#         required=False,
#     )
#     parser.add_argument(
#         "--exclude_families",
#         type=str,
#         help="TE families to exclude in the phylogeny (separated by comma)",
#         required=False,
#     )
#     parser.add_argument(
#         "--exclude_samples",
#         type=str,
#         help="sample to exclude in the phylogeny (separated by comma)",
#         required=False,
#     )
#     parser.add_argument(
#         "--thread",
#         type=int,
#         help="max cpu threads to use (default = '1')",
#         required=False,
#     )

#     args = parser.parse_args()

#     # sets up default value for optional variable
#     if thread is None:
#         thread = 1

#     return args


def mkdir(dir):
    if os.path.isdir(dir):
        print("Directory %s exists" % dir)
        return
    try:
        os.mkdir(dir)
    except OSError:
        print("Creation of the directory %s failed" % dir)
    else:
        print("Successfully created the directory %s " % dir)


def merge_bed(bed_list, bed_out):
    with open(bed_out, "w") as output:
        for bed in bed_list:
            sample_name = os.path.basename(bed).replace(".nonref.bed", "")
            with open(bed, "r") as input:
                for line in input:
                    entry = line.replace("\n", "").split("\t")
                    chrom = entry[0]
                    start = entry[1]
                    end = entry[2]
                    family = entry[3]
                    score = entry[4]
                    strand = entry[5]
                    info = "|".join([family, sample_name])
                    out_line = "\t".join([chrom, start, end, info, score, strand])
                    output.write(out_line + "\n")


def filter_region_bed(bed_in, region_filter, bed_out):
    with open(bed_out, "w") as output:
        subprocess.call(
            ["bedtools", "intersect", "-a", bed_in, "-b", region_filter, "-u"],
            stdout=output,
        )


def filter_family_bed(bed_in, family_filter, bed_out, method):
    families = set(family_filter.split(","))
    with open(bed_out, "w") as output, open(bed_in, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            family = entry[3].split("|")[0]
            if method == "include":
                if family in families:
                    output.write(line)
            else:
                if family not in families:
                    output.write(line)


def sort_bed(bed_in, bed_out):
    with open(bed_out, "w") as output:
        subprocess.call(["bedtools", "sort", "-i", bed_in], stdout=output)


def cluster_bed(bed_in, bed_out):
    window = 0
    with open(bed_out, "w") as output:
        subprocess.call(
            ["bedtools", "cluster", "-s", "-d", str(window), "-i", bed_in],
            stdout=output,
        )


def get_te_info(meta):
    te_info_dict = dict()
    with open(meta, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            te_info_dict[entry[0]] = entry[1]
    return te_info_dict


def get_method_bed(
    te_out_dirs,
    filter_region,
    outdir,
    include_families,
    exclude_families,
    exclude_samples,
    prefix,
):
    output_prefix = prefix
    pattern = "/**/detect/*nonref.bed"
    bed_list = []
    for te_out_dir in te_out_dirs:
        bed_files = glob.glob(te_out_dir + pattern, recursive=True)
        bed_list = bed_list + bed_files

    # exclude samples
    if exclude_samples is not None:
        exclude_sample_list = exclude_samples.replace(" ", "").split(",")
        new_bed_list = []
        for bed in bed_list:
            include = True
            for exclude_sample in exclude_sample_list:
                if exclude_sample in bed:
                    include = False
            if include:
                new_bed_list.append(bed)
        bed_list = new_bed_list

    # keep only nonref, merge per method
    bed_merged = os.path.join(outdir, output_prefix + ".merge.bed")
    merge_bed(bed_list=bed_list, bed_out=bed_merged)

    # filter by region
    bed_filtered = os.path.join(outdir, output_prefix + ".filter.bed")
    filter_region_bed(
        bed_in=bed_merged,
        region_filter=filter_region,
        bed_out=bed_filtered,
    )

    # filter by family
    if include_families is not None:
        bed_filtered_tmp = bed_filtered + ".tmp"
        filter_family_bed(
            bed_in=bed_filtered,
            family_filter=include_families,
            bed_out=bed_filtered_tmp,
            method="include",
        )
        os.rename(bed_filtered_tmp, bed_filtered)

    if exclude_families is not None:
        bed_filtered_tmp = bed_filtered + ".tmp"
        filter_family_bed(
            bed_in=bed_filtered,
            family_filter=exclude_families,
            bed_out=bed_filtered_tmp,
            method="exclude",
        )
        os.rename(bed_filtered_tmp, bed_filtered)

    # sort and cluster, use method specific window
    bed_sort = os.path.join(outdir, output_prefix + ".sort.bed")
    sort_bed(bed_in=bed_filtered, bed_out=bed_sort)

    bed_cluster = os.path.join(outdir, output_prefix + ".cluster.bed")
    cluster_bed(bed_in=bed_sort, bed_out=bed_cluster)

    # os.remove(bed_merged)
    # os.remove(bed_filtered)
    # os.remove(bed_sort)

    return bed_cluster


def bed2matrix(bed, outgroup):
    # from clustered bed to binary data matrix
    header = [
        "chr",
        "start",
        "end",
        "info",
        "score",
        "strand",
        "cluster",
    ]
    df = pd.read_csv(bed, delimiter="\t", names=header)
    info = df["info"].str.split("|", expand=True)
    df["family"] = info[0]
    df["sample"] = info[1]
    df.drop(["info"], inplace=True, axis=1)

    ## filter out clusters with where a sample appears > 1 time
    df.drop_duplicates(subset=["sample", "cluster"], keep=False, inplace=True)

    ## filter out clusters with more than two TE families
    df = df.groupby("cluster").filter(lambda x: x["family"].nunique() == 1)

    # convert to matrix
    df["value"] = 1
    matrix = df.pivot_table(
        index="sample", columns="cluster", values="value", fill_value=0
    )
    if outgroup is not None and outgroup != "":
        root_row = pd.Series(0, index=matrix.columns)
        root_row.name = outgroup
        matrix = matrix.append(root_row)
    return matrix


def get_cluster(input, outdir, prefix):
    tree_file = os.path.join(outdir, prefix + ".nwk")
    pdf_file = os.path.join(outdir, prefix + ".pdf")
    subprocess.call(
        ["Rscript", "--vanilla", "get_cluster.R", input, tree_file, pdf_file]
    )

    return tree_file, pdf_file


def cluster(
    prefix,
    enrichmente_out_dirs,
    filter_region,
    outgroup,
    out,
    include_families,
    exclude_families,
    exclude_samples,
):
    # args = get_args()

    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    # get clustered bed from each matrix
    bed_cluster = get_method_bed(
        enrichmente_out_dirs,
        filter_region,
        out,
        include_families,
        exclude_families,
        exclude_samples,
        prefix,
    )
    # generate data matrix
    matrix = bed2matrix(bed_cluster, outgroup=outgroup)
    
    # write matrix to csv file
    matrix_file = os.path.join(out, prefix + ".matrix.csv")
    matrix.to_csv(
        matrix_file,
        sep=",",
        index=True,
        header=True,
    )

    # generate NJ tree
    tree_file, pdf_file = get_cluster(matrix_file, out, prefix)

    return tree_file, pdf_file


# main()
