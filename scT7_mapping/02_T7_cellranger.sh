#########################################
# Running cellranger count on T7 IST+PAP 
# capture data for the purpose of cell
# barcode and UMI correction
#########################################

# Definining paths
ref=<path_to_cellranger_reference>
work_dir=<path_to_working_directory>
fastq1=${work_dir}/<path_to_first_sequencing_run_fastq_directory>
fastq2=${work_dir}/<path_to_second_sequencing_run_fastq_directory>
outdir=${work_dir}/<path_to_count_output_directory>

# Loading cellranger
module load cellranger/8.0.0

cd $outdir

# prefixes.fofn: plain text file (located alongside the first run's fastqs), listing
# one sample prefix per line, e.g. t7_K562_L1_S1
# This includes the "_S#" sample index assigned during demultiplexing; the sed
# command below strips it to get the sample_id used for both --id and --sample.
#
# NOTE: both sequencing runs (250825 and 250902) were demultiplexed with the same
# sample sheet, so the same "_S#" index is preserved for a given sample across both
# run dates internally (e.g. t7_K562_L1_S1 appears in both run directories).
# GEO-deposited fastqs are named differently: they drop the "_S#" index and "_001"
# suffix entirely and instead encode the run date, e.g.
# t7_K562_L1_S1_R1_001.fastq.gz (internal, either run) ->
# Top100v2_scT7capture_K562_L1_run250825_R1.fastq.gz /
# Top100v2_scT7capture_K562_L1_run250902_R1.fastq.gz (GEO, per run)
# Confirm/adjust naming accordingly if not using the original internal fastq
# directory structure.
for sample in $(cat ${fastq1}/prefixes.fofn); do
  echo "Processing sample: $sample"

  sample_id=$(echo $sample | sed 's/_S[0-9]\+//')

  cellranger count \
    --id ${sample_id} \
    --transcriptome $ref \
    --fastqs ${fastq1},${fastq2} \
    --sample ${sample_id} \
    --create-bam=true \
    --chemistry=auto  \
    --nosecondary \
    --localcores 16 \
    --localmem   240

done


