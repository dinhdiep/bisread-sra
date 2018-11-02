# Bisread-sra


## View the pipeline

```bash
snakemake --forceall --dag | dot -Tpng > dag.png
```

## Preliminary steps

Make sure that [Miniconda 3] is installed. When prompted, answer yes to add conda to the PATH environment variable.
```bash
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
```

Download bisulfite-seq pipeline
```bash
   git clone https://github.com/hdinhdp/bisulfite-seq.git
```

Go to the bisread-sra directory
```bash
   cd bisulfite-seq/bisread-sra
```

Use the yaml file to install all required software into an isolated Conda environment with the name [bisread-preprocess]
```bash
   conda env create --name bisread-sra --file environment.yml
```

Activate the environment
```bash
   source activate bisread-sra
```

Note: to exit the environment, just execute the following command. Don't do this until you are done working.
```bash
   source deactivate
```

## Generating samples.json and configuring config.yml

Make the samples.json file with all the available fastq files in data/fastq folder. Note that ```make_samples.py``` might need to be modified in order to recognize your fastq files and guess the sample names. Keep sample names that corresponds to the fastq files for now as they should ideally be merged in the calling methylation frequency stage. This script assumes that read 1 all have ```_R1_``` in their names and read 2 all have ```_R2_``` in their names. The sample names are the first substring of the file name separated by a ```_```.

```bash
   python3 ../scripts/make_samples.py
```
The current ```config.yml``` file may not be compatible with the library type. Make sure to modify this file using a text editor to modify the following lines.

```yaml
# Set the trimming mode to be BSPP, WGBS, or RRBS for this run
TRIM_MODE: 'BSPP'

# Set the UMI sequence length to be extracted from read 1
UMI_LENGTH: 10

# Set the read 1 adapter sequence, use default 'AGATCGGAAGAGC' if you don't know/have any. Do not leave empty.
ADAPTER1: 'AGATCGGAAGAGC'

# Set the read 2 adapter sequence, use default 'AGATCGGAAGAGC' if you don't know/have any. Do not leave empty.
ADAPTER2: 'AGATCGGAAGAGC'

# Set the quality base for Illumina reads 64 or 33
QUAL_BASE: 64

```

## Run the snakemake

It is possible to speedup the process by allowing snakemake to spawn multiple jobs in parallel whenever possible by specifying the job number.
```bash 
   snakemake --jobs 4
```

