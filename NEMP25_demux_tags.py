import os

#Define variables
samples = ['NEMP25_GEX_t1_1', 'NEMP25_GEX_t1_2', 'NEMP25_GEX_t1_3', 'NEMP25_GEX_t1_4', 'NEMP25_GEX_t2_1', 'NEMP25_GEX_t2_2', 'NEMP25_GEX_t2_3', 'NEMP25_GEX_t2_4', 'NEMP25_GEX_t3_1', 'NEMP25_GEX_t3_2', 'NEMP25_GEX_t3_3', 'NEMP25_GEX_t3_4', 'NEMP25_GEX_t3_5', 'NEMP25_GEX_t3_6']
fastq_dir = "/wynton/group/pollen/reads/RMAP05" 
whitelist = "/wynton/home/pollenlab/reedmcmullen/libraries/3M-february-2018.txt"
names = "/wynton/home/pollenlab/reedmcmullen/libraries/NEMP25_MS_names.txt"
seqs = "/wynton/home/pollenlab/reedmcmullen/libraries/NEMP25_MS_seqs.txt"
program_dir = "/wynton/home/pollenlab/reedmcmullen/utils/cellbouncer"   
out_dir = "/wynton/scratch/rmcmullen/NEMP25_demux_tags/Outs"

#Set working directory to the program directory.
os.chdir(program_dir)

#Define a function to find MULTI-seq fastq files for a given sample and prepend --read1 or --read2, as necessary.
def find_fastq_files(directory, sample):
    # List to hold the matched file names
    matched_files = []

    # Traverse the directory
    for filename in os.listdir(directory):
        # Check if the file name contains 'sample' and '_MS'
        if sample in filename and '_MS' in filename:
            # Check for R1_001 or R2_001
            if 'R1_001' in filename:
                matched_files.append(f'--read1 {filename}')
            elif 'R2_001' in filename:
                matched_files.append(f'--read2 {filename}')
    matched_files.sort()
    return matched_files
  
for sample in samples:
  find_fastq_files(fastq_dir, sample)
  ./demux_tags -o f'{out_dir}/{sample}' -w whitelist -N names -s seqs -m 2 -u 12 matched_files
