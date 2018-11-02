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

Download bisread-sra pipeline
```bash
   git clone https://github.com/hdinhdp/bisread-sra.git
```

Go to the bisread-sra directory
```bash
   cd bisread-sra
```

Use the yaml file to install all required software into an isolated Conda environment with the name **bisread-sra**
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

Make the samples.json file with all the SRA ids.

The current ```config.yml``` file may not be compatible with the library type. Make sure to modify this file using a text editor to your specific
library parameters if this is not a WGBS data.


## Run the snakemake

It is possible to speedup the process by allowing snakemake to spawn multiple jobs in parallel whenever possible by specifying the job number.
```bash 
   snakemake --jobs 4
```

