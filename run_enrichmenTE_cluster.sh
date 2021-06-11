#!/bin/bash
#SBATCH --job-name=enrichmenTE_test
#SBATCH --partition=bergman_p
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sh60271@uga.edu
#SBATCH --output=/scratch/sh60271/enrichmenTE/test_0611.out
#SBATCH --export=NONE

source activate enrichmenTE

# prepare input data
enrichmenTE_dir="/home/sh60271/git/enrichmenTE"
norm_region="/home/sh60271/git/DGRC/data/annotation/recomb_border_r6.bed"
run_dir="/scratch/sh60271/enrichmenTE/test_0611"

# set up directory for the job
mkdir -p $run_dir
data_dir=$run_dir/data
out_dir=$run_dir/out
script_dir=$run_dir/script
mkdir -p $data_dir
mkdir -p $out_dir
mkdir -p $script_dir
cd $out_dir

# ######## cluster ########
python3 $enrichmenTE_dir/enrichmenTE.py cluster --prefix "clustergram" --enrichmente_out_dirs $out_dir --filter_region $norm_region --out $out_dir --thread $SLURM_NTASKS --exclude_families "mdg3"