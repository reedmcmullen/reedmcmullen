#!/bin/bash
#$ -V
#$ -S /bin/bash                                                 # Specifies shell for job script as /bin/bash
#$ -e /wynton/scratch/rmcmullen/NEMP25_demux_tags		            # Directs error output to specified directory.
#$ -o /wynton/scratch/rmcmullen/NEMP25_demux_tags		            # Directs output to the specified directory.
#$ -j y                                                         # Joins error and output into a single file.
#$ -cwd                                                         # Executes job from the current working dir.
#$ -l mem_free=50G                                              # Requests n Gigabytes of free memory.
#$ -l h_rt=12:00:00                                             # Requests a hard runtime limit.
#$ -l scratch=50G                                               # Requests n Gigabytes of scratch space.
#$ -pe smp 7                                                    # Requests a parallel env. with n slots.
#$ -t 1-14                                                      # Define the array job range.

# Define the required variables.
FASTQ_DIR="/wynton/group/pollen/reads/RMAP05"  					                        # Directory containing the fastq files.
PROGRAM_DIR="/wynton/home/pollenlab/reedmcmullen/utils/cellbouncer"  		        # Directory containing the demux_species program.
WHITELIST="/wynton/home/pollenlab/reedmcmullen/libraries/3M-february-2018.txt"  # File of 10X Genomics' whitelist for scRNA-seq GEX barcodes.
LIBS_FILE="/wynton/home/pollenlab/reedmcmullen/libraries/NEMP25_MS_samples.txt" # File of MULTI-seq samples, one per line.
SEQS_FILE="/wynton/home/pollenlab/reedmcmullen/libraries/NEMP25_MS_seqs.txt"	  # File of MULTI-seq sequences.
NAMES_FILE="/wynton/home/pollenlab/reedmcmullen/libraries/NEMP25_MS_names.txt"  # File of MULTI-seq names.
OUT_DIR="/wynton/scratch/rmcmullen/NEMP25_demux_tags/Outs"

# Activate the conda environment.
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cellbouncer

# Make the output directory and set the working directory to the program directory.
mkdir -p "$OUT_DIR"
cd "$PROGRAM_DIR"

# Find paired read1 and read2 files for each unique sample.
SAMPLE=$(head -n $SGE_TASK_ID "$LIBS_FILE" | tail -1)
SAMPLE_FILES=$(find "$FASTQ_DIR" -type f -name "${SAMPLE}*R[12]_001*" | sort)
READS_ARGS="${OUT_DIR}/${SAMPLE}_reads.txt"
> "$READS_ARGS"  # Clear or create the file

# Process each file and append to the output file with the correct prefix.
for file in $SAMPLE_FILES; do
    if [[ "$file" == *"R1_001"* ]]; then
        echo "--read1 $file" >> "$READS_ARGS"
    elif [[ "$file" == *"R2_001"* ]]; then
        echo "--read2 $file" >> "$READS_ARGS"
    fi
done

# Output the content of the generated reads file
cat "$READS_ARGS"

# Run the demux_tags program.
./demux_tags \
--output "${OUT_DIR}/${SAMPLE}" \
$(cat "$READS_ARGS") \
--whitelist "$WHITELIST" \
--names "$NAMES_FILE" \

# Deactivate the conda environment.
conda deactivate
