# Working directory
work_dir="<path_to_working_directory>/"
bam_L1="${work_dir}count_RNA/RNA_K562_L1/outs/possorted_genome_bam.bam"
bam_L2="${work_dir}count_RNA/RNA_K562_L2/outs/possorted_genome_bam.bam"
out_dir="${work_dir}count_RNA/subsetted_bam"
ref="<path_to_cellranger_reference>"



cd $work_dir
# Load samtools
module load samtools/1.19
module load R/4.3.2
module load bedtools/2.31.1
module load cellranger/8.0.0

# Create output directory
mkdir -p "$out_dir"

#######################################
# Step 1: filter T7 reads from bam    #
#######################################

# Getting a list of targeted CRE
# regions that have been expanded by 
# 1Kb on each side
Rscript ./scripts/blacklist_t7_ranges.R

echo Filtering T7 reads from lane 1
bedtools intersect -abam $bam_L1 -b blacklist_t7_ranges.bed -v > ${out_dir}/t7_filtered.L1.possorted_genome.bam
samtools index ${out_dir}/t7_filtered.L1.possorted_genome.bam
echo
echo Filtering T7 reads from lane 2
bedtools intersect -abam $bam_L2 -b blacklist_t7_ranges.bed -v > ${out_dir}/t7_filtered.L2.possorted_genome.bam
samtools index ${out_dir}/t7_filtered.L2.possorted_genome.bam
echo Finished

#######################################
# Step 2: convert filtered bam back   #
# to fastq files to run cellranger    #
# count on                            #
#######################################

# Output FASTQ dirs
fastq_L1_dir="${out_dir}/t7-filtered-fastq_L1"
fastq_L2_dir="${out_dir}/t7-filtered-fastq_L2"

echo "Converting filtered BAM (L1) -> FASTQ..."
# --relaxed tolerates unpaired/duplicate reads after filtering
# Increase threads as appropriate for your node
cellranger bamtofastq \
  --nthreads=16 \
  --relaxed \
  "${out_dir}/t7_filtered.L1.possorted_genome.bam" \
  "$fastq_L1_dir"

echo
echo "Converting filtered BAM (L2) -> FASTQ..."
cellranger bamtofastq \
  --nthreads=16 \
  --relaxed \
  "${out_dir}/t7_filtered.L2.possorted_genome.bam" \
  "$fastq_L2_dir"

echo
echo "bamtofastq complete."
echo "FASTQs are in:"
echo "  $fastq_L1_dir"
echo "  $fastq_L2_dir"

#######################################
# Step 3: re-running cellranger count #
# one job per lane, two flowcells     #
#######################################

# Output FASTQ dirs
fastq_L1_dir="${out_dir}/t7-filtered-fastq_L1"
fastq_L2_dir="${out_dir}/t7-filtered-fastq_L2"

recount_out="${work_dir}count_RNA/t7_filtered_counts"
mkdir -p "$recount_out"
cd "$recount_out"

# Sample name is taken from FASTQ prefix 'bamtofastq_*'
sample_name="bamtofastq"

# Lane 1: two flowcell dirs (NO SPACES after comma)
fastqs_L1="${fastq_L1_dir}/RNA_K562_L1_0_1_AACYC32HV,${fastq_L1_dir}/RNA_K562_L1_0_1_AAHJTC2M5"

echo "=== cellranger count: t7-filtered_RNA_L1 ==="
cellranger count \
  --id=t7-filtered_RNA_L1 \
  --transcriptome="$ref" \
  --fastqs="$fastqs_L1" \
  --sample="$sample_name" \
  --create-bam=true \
  --chemistry=auto \
  --nosecondary \
  --expect-cells=20000 \
  --localcores=24 \
  --localmem=495

# Lane 2: two flowcell dirs (NO SPACES after comma)
fastqs_L2="${fastq_L2_dir}/RNA_K562_L2_0_1_AACYC32HV,${fastq_L2_dir}/RNA_K562_L2_0_1_AAHJTC2M5"

echo
echo "=== cellranger count: t7-filtered_RNA_L2 ==="
cellranger count \
  --id=t7-filtered_RNA_L2 \
  --transcriptome="$ref" \
  --fastqs="$fastqs_L2" \
  --sample="$sample_name" \
  --create-bam=true \
  --chemistry=auto \
  --nosecondary \
  --expect-cells=20000 \
  --localcores=24 \
  --localmem=495

echo
echo "Done. Outputs in: $recount_out"
