# --- Set these paths for your environment before running ---
work_dir=<path_to_fastq_directory>
bam_dir=<path_to_bam_output_directory>
env=<path_to_conda_environment>
index=<path_to_bwa_genome_index>

# Activate the conda environment
source activate $env

# Load the samtools module
module load samtools/1.19

# Change to the working directory
cd $work_dir

# prefixes.fofn: plain text file, one fastq filename prefix per line (one per sample/replicate),
# matching the fastq files as deposited on GEO for this dataset.
# e.g. Top100v2_T7_K562_3N_3d_PEmax_rep1
for prefix in $(cat prefixes.fofn); do

    # NOTE: _R1_001 / _R2_001 suffixes below match Illumina bcl2fastq output naming used internally.
    # GEO-deposited fastqs are named <prefix>_R1.fastq.gz / <prefix>_R2.fastq.gz (no "_001" suffix).
    # Confirm/adjust these suffixes to match actual filenames accordingly.
    R1_file=${prefix}_R1_001.fastq.gz
    R2_file=${prefix}_R2_001.fastq.gz

    echo "Aligning $prefix to hg38"
    bwa mem -t 16 -Y $index $R1_file $R2_file | samtools view -b -o ${bam_dir}${prefix}.bam
    
    echo "Sorting and indexing $prefix"
    samtools sort -o ${bam_dir}${prefix}.sorted.bam ${bam_dir}${prefix}.bam
    samtools index ${bam_dir}${prefix}.sorted.bam

    echo "Finished $prefix"
    echo
done


