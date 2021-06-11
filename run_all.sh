#!/bin/bash
#SBATCH --job-name=enrichmenTE_test
#SBATCH --partition=bergman_p
#SBATCH --nodes=1
#SBATCH --tasks-per-node=28
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sh60271@uga.edu
#SBATCH --output=/scratch/sh60271/enrichmenTE/test_0611.out
#SBATCH --export=NONE

source activate enrichmenTE

# prepare input data
enrichmenTE_dir="/home/sh60271/git/enrichmenTE"
read_dir="/work/cmblab/cbergman/dgrc"
norm_region="/home/sh60271/git/DGRC/data/annotation/recomb_border_r6.bed"
curated_annotation="/home/sh60271/git/DGRC/data/annotation/dm6_te.gff"
mask_aug_ref="/work/cmblab/shunhua/dm6_masked_augmented/dm6.cov.fasta.masked.aug"
run_dir="/scratch/sh60271/enrichmenTE/test_0611"
depth_config=$enrichmenTE_dir/utility/cutoff.config

# set up directory for the job
mkdir -p $run_dir
data_dir=$run_dir/data
out_dir=$run_dir/out
script_dir=$run_dir/script
mkdir -p $data_dir
mkdir -p $out_dir
mkdir -p $script_dir
cd $out_dir

# prepare index for modified dm6
if [ ! -f $data_dir/dm6.masked.aug.fasta ]; then
	cp $mask_aug_ref $data_dir/dm6.masked.aug.fasta
	chmod u+x $data_dir/dm6.masked.aug.fasta
	bwa index $data_dir/dm6.masked.aug.fasta
fi
ref_masked_aug=$data_dir/dm6.masked.aug.fasta

if [[ ! -f $ref_masked_aug".fai" ]]; then
	samtools faidx $ref_masked_aug
fi
if [[ ! -f $ref_masked_aug".bwt" ]]; then
	bwa index $ref_masked_aug
fi

# run
cd $script_dir
echo "#!/bin/bash" >jobscript.sh
echo "#SBATCH --partition=batch" >>jobscript.sh
echo "#SBATCH --nodes=1" >>jobscript.sh
echo "#SBATCH --tasks-per-node=24" >>jobscript.sh
echo "#SBATCH --mem=50G" >>jobscript.sh
echo "#SBATCH --time=10:00:00" >>jobscript.sh
echo "#SBATCH --export=None" >>jobscript.sh
echo "source activate enrichmenTE" >>jobscript.sh
echo "python3 \$enrichmenTE_dir/enrichmenTE.py detect -1 \$read1 -2 \$read2 -o \$out_dir -f \$filter_region -r \$ref_masked_aug -g \$curated_annotation -d \$depth_config -t \$SLURM_NTASKS --prefix \$prefix" >>jobscript.sh

for file in $(ls -d -1 /work/cmblab/shunhua/dgrc/*_R1_001.fastq.gz); do
	sample_name=$(basename "$file" _R1_001.fastq.gz)
	fq1=$file
	fq2=$(echo ${file} | sed -E "s/_R1_/_R2_/g")
	log_output=$out_dir/$sample_name.log
	sbatch --job-name=${prefix} --export=enrichmenTE_dir="${enrichmenTE_dir}",read1="${fq1}",read2="${fq2}",out_dir="${out_dir}",filter_region="${norm_region}",ref_masked_aug="${ref_masked_aug}",curated_annotation="${curated_annotation}",depth_config="${depth_config}",prefix="${sample_name}" --output="${log_output}" --error="${log_output}" jobscript.sh
done

# # ######## cluster ########
# tree_dir=$out_dir/tree2
# mkdir -p $tree_dir
# # python3 $enrichmenTE_dir/enrichmenTE_cluster.py --prefix "test_cluster" --enrichmente_out_dirs $out_dir --filter_region $norm_region --outgroup "ISO1" --out $tree_dir --thread $SLURM_NTASKS --exclude_families "mdg3"
# python3 $enrichmenTE_dir/enrichmenTE.py cluster --prefix "cluster" --enrichmente_out_dirs $out_dir --filter_region $norm_region --out $out_dir --thread $SLURM_NTASKS --exclude_families "mdg3"

# source activate te_ngs

# utility_dir="/home/sh60271/git/DGRC/src/utility"
# run_dir="/scratch/sh60271/primer/read_0514_batch4"

# # merge and cluster the bed files
# bed_dir="/scratch/sh60271/primer/read_0514_batch4/out"
# out_dir=$run_dir/out/merge_cluster
# mkdir -p $out_dir
# python3 $utility_dir/merge_cluster.py -d $bed_dir -o $out_dir

# bed="/scratch/sh60271/enrichmenTE/test_0607/out/tree/test_cluster.cluster.bed"
# # produce neighbor joining tree
# png=$out_dir/merge.cluster.png
# tree=$out_dir/merge.cluster.nw
# Rscript --vanilla $enrichmenTE_dir/nj.R $bed $png $tree
