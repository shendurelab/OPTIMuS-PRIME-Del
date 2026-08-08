work_dir=<path_to_working_directory>/
fastq_dir=${work_dir}250410_fastq/
bam_dir=${work_dir}250410_alignments/

# CRE-deletions.fa has a reference for each of the two
# barcodes that were designed for any given unique pair
# of pegRNAs
cp ${work_dir}CRE-deletions.fa ${bam_dir}CRE-deletions.fa
index=${bam_dir}CRE-deletions.fa 

env=<path_to_conda_environment>/
# Activate the conda environment
source activate $env

# Load the samtools module
module load samtools/1.19

# Indexing reference
echo Indexing reference
cd $bam_dir
bwa index $index
echo 

# Change to the working directory
cd $fastq_dir

# Loop through the prefixes listed in the file
for prefix in $(cat prefixes.fofn); do

    echo "Aligning reads for sample $prefix"
    bwa mem -t 4 $index ${prefix}_R1_001.fastq.gz ${prefix}_R2_001.fastq.gz | samtools view -o ${bam_dir}${prefix}.bam
    
    echo "Sorting and indexing .bam file for sample $prefix"
    samtools sort -o ${bam_dir}${prefix}.sorted.bam ${bam_dir}${prefix}.bam
    samtools index ${bam_dir}${prefix}.sorted.bam

    echo "Finished $prefix"
    echo
done

