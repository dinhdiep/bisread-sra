'''
Snakefile
bisread-sra

Analysis of bisulfite reads from sra Run Ids

2018-10-30 Copyright Dinh Diep

'''
import json
from os.path import join, basename, dirname
from subprocess import check_output
from itertools import chain

# Globals ---------------------------------------------------------------------

# Configurations for this run is indicated in this file.
configfile: 'config.yml'

# Full path to a folder where intermediate output files will be created.
FASTQ_DIR = config['FASTQ_DIR']

# Full path to a folder where the intermediate BAM files will be created.
BAM_DIR = config['BAM_DIR']

# Full path to the output methylation directory.
METHYL_DIR = config['METHYL_DIR']

# Full path to FASTA file with all the chromosome sequences.
GENOME = config['GENOME']

quality_base = config['QUAL_BASE']

min_map_score = config['MAP_SCORE']

# Trim mode set
TRIM_MODE = config['TRIM_MODE']

# Adapter 1 set
ADAPTER1 = config['ADAPTER1']

# Adapter 2 set
ADAPTER2 = config['ADAPTER2']

# Samples and their corresponding filenames.
FILES = json.load(open(config['SAMPLES_JSON']))
SAMPLES = sorted(FILES.keys())

# Cutadapt parameters
error_rate = config['ERROR_RATE']
quality_3prime = config['QUAL_CUTOFF']
minimum_length = config['MINIMUM_LENGTH']
overlap_length = config['REQUIRED_OVERLAP']
quality_base = config['QUAL_BASE']

cutadapt1 = '-f fastq -e ' \
            + str(error_rate) \
            + ' -q ' + str(quality_3prime) \
            + ' -m ' + str(minimum_length) \
            + ' --overlap=' + str(overlap_length) \
            + ' --quality-base=' \
            + str(quality_base) \
            + ' -a ' + str(ADAPTER1)

cutadapt2 = '-f fastq -e ' \
            + str(error_rate) \
            + ' -q ' + str(quality_3prime) \
            + ' -m ' + str(minimum_length) \
            + ' --overlap=' + str(overlap_length) \
            + ' --quality-base=' \
            + str(quality_base) \
            + ' -a ' + str(ADAPTER2)

# Do we need to extract the UMI from read 1?
UMI_LENGTH = config['UMI_LENGTH']

if TRIM_MODE == 'BSPP':
    minimum_length = minimum_length + 27
else:
    minimum_length = minimum_length + 5

# Functions -------------------------------------------------------------------

def get_mode(file1, file2):
    mode = 'paired'
    if file2[0] == "":
        mode = 'single'
    return mode


# Rules -----------------------------------------------------------------------

rule all:
    input:
        'mapping_report.txt'


rule genome_prep:
    input:
        g = GENOME
    output:
        fai = GENOME + '.fai',
        ct = GENOME + '.bis.CT',
        ga = GENOME + '.bis.GA',
        cpg = GENOME + '.cpg.positions.txt'
    run:
        # verify the GENOME file first
        shell('samtools faidx {GENOME}')
        shell('bin/genomePrep.pl {input.g} convert=yes context=cg')
        shell('cat {input.g}.chr*cpositions.txt > {output.cpg}')


rule minimap_index:
    input:
        ct = GENOME + '.bis.CT',
        ga = GENOME + '.bis.GA'
    output:
        ct = GENOME + '.bis.CT.mmi',
        ga = GENOME + '.bis.GA.mmi'
    log:
        join(dirname(GENOME), 'minimap2_index.log')
    run:
        shell('minimap2 -x sr -d {output.ct} {input.ct} >> {log} 2>&1')
        shell('minimap2 -x sr -d {output.ga} {input.ga} >> {log} 2>&1')


rule samtools_faidx:
    input:
        ct = GENOME + '.bis.CT',
        ga = GENOME + '.bis.GA'
    output:
        GENOME + '.bis.CT.fai',
        GENOME + '.bis.GA.fai'
    run:
        shell('samtools faidx {input.ct}')
        shell('samtools faidx {input.ga}')


rule fastq_dump:
    output: 
        temp(join(FASTQ_DIR, '{sample}.fastq.lst'))
    log:
        join('log', '{sample}.fqdump.log')
    run:
        r1 = FILES[wildcards.sample]['R1']
	shell('fastq-dump --split-files {r1} --outdir {FASTQ_DIR} > {log} 2>&1')
        shell('ls {FASTQ_DIR}/{r1}* > {output}')

        
rule trim_umi:
    input:
        join(FASTQ_DIR, '{sample}.fastq.lst')
    output:
        temp(join(FASTQ_DIR, '{sample}.umi.lst'))
    run:
        r1 = FILES[wildcards.sample]['R1']
        read_files = [line.rstrip('\n') for line in open(input[0])]
        mode = get_mode(read_files[0], read_files[1])
        umi_ex = FASTQ_DIR + '/' + wildcards.sample + '.' + mode + '.umi'
        read_files = [line.rstrip('\n') for line in open(input[0])]
        if mode == 'paired':
            if UMI_LENGTH > 0:
                shell('bin/getUMI.pl {read_files[0]} {read_files[1]} {umi_ex} {UMI_LENGTH}')
                shell('rm {read_files[0]} {read_files[1]}');
                shell('echo {umi_ex}.R1.fq > {output}')
                shell('echo {umi_ex}.R2.fq >> {output}')
            else:
                shell('echo {read_files[0]} > {output}')
                shell('echo {read_files[1]} >> {output}')
        else:
            if UMI_LENGTH > 0:
                shell('bin/getUMI.pl {read_files[0]} {umi_ex} {UMI_LENGTH}')
                shell('rm {read_files[0]}');
                shell('echo {umi_ex}.R1.fq > {output}')
            else:
                shell('echo {read_files[0]} > {output}')


rule trim_cutadapt:
    input: 
        join(FASTQ_DIR, '{sample}.umi.lst')
    output:
        temp(join(FASTQ_DIR, '{sample}.cutadapt.lst'))
    log:
        l1 = join('log', '{sample}_cutadapt_R1.log'),
        l2 = join('log', '{sample}_cutadapt_R2.log')
    run:
        r1 = FILES[wildcards.sample]['R1']
        read_files = [line.rstrip('\n') for line in open(input[0])]
        mode = get_mode(read_files[0], read_files[1])
        trimmed = FASTQ_DIR + '/' + wildcards.sample + '.' + mode + '.trimmed'
        if mode == 'paired':
            shell('cutadapt {cutadapt1} -o {trimmed}.R1.fq {read_files[0]} > {log.l1} 2>&1')
            shell('cutadapt {cutadapt2} -o {trimmed}.R2.fq {read_files[1]} > {log.l2} 2>&1')
            shell('rm {read_files[0]} {read_files[1]}');
            shell('echo {trimmed}.R1.fq {log.l1} > {output}')
            shell('echo {trimmed}.R2.fq {log.l2} >> {output}')
        else:
            shell('cutadapt {cutadapt1} -o {trimmed}.R1.fq {read_files[0]} > {log.l1} 2>&1')
            shell('rm {read_files[0]}');
            shell('echo {trimmed}.R1.fq {log.l1} > {output}')


rule encode:
    input:
        join(FASTQ_DIR, '{sample}.cutadapt.lst')
    output:
        temp(join(FASTQ_DIR, '{sample}.encoded.lst'))
    run:
        r1 = FILES[wildcards.sample]['R1']
        read_files = [line.rstrip('\n').split()[0] for line in open(input[0])]
        mode = get_mode(read_files[0], read_files[1])
        encoded = FASTQ_DIR + '/' + wildcards.sample + '.' + mode + '.encoded'
        proper = FASTQ_DIR + '/' + wildcards.sample + '.' + mode + '.proper'
        if mode == 'paired':
            shell('cat {read_files[0]} {read_files[1]} | paste - - - - | sort -k1,1 -T {FASTQ_DIR} | bin/identify_pairs.pl {proper}') 
            shell('bin/encodeFastq.pl {proper}.R1.fq {encoded}.R1.fq {TRIM_MODE}')
            shell('bin/encodeFastq.pl {proper}.R2.fq {encoded}.R2.fq {TRIM_MODE}')
            shell('rm {proper}.* {read_files[0]} {read_files[1]}')
            shell('echo {encoded}.R1.fq > {output}')
            shell('echo {encoded}.R2.fq >> {output}')
        else:
            shell('bin/encodeFastq.pl {read_files[0]} {encoded}.R1.fq {TRIM_MODE}')
            shell('rm {read_files[0]}')
            shell('echo {encoded}.R1.fq > {output}')
     

rule minimap2:
    input:
        expand(GENOME + '.bis{suffix}', suffix = ['.CT.mmi', '.GA.mmi', '.CT.fai', '.GA.fai']),
	lst = join(FASTQ_DIR, '{sample}.encoded.lst')
    output:
        sam_ct = temp(join(BAM_DIR, '{sample}.raw.CT.sam')),
        sam_ga = temp(join(BAM_DIR, '{sample}.raw.GA.sam')),
        o1 = temp(join('temp', '{sample}.fastq'))
    threads: 24
    benchmark:
        join('benchmark', '{sample}.minimap.benchmark.txt')
    log:
        join('log', '{sample}.minimap2.log')
    params:
        rg = '@RG\\tID:{sample}\\tSM:{sample}'
    run:
        read_files = expand('{file}', file = [line.rstrip('\n') for line in open(input.lst)])
        shell('cat {read_files} > {output.o1}')
        shell("minimap2 --secondary=no -a -x sr -R '{params.rg}' -t {threads} {GENOME}.bis.CT.mmi {output.o1} > {output.sam_ct} 2> {log}")
        shell("minimap2 --secondary=no -a -x sr -R '{params.rg}' -t {threads} {GENOME}.bis.GA.mmi {output.o1} > {output.sam_ga} 2>> {log}")
        
        #Uncomment this after testing:
        #shell('rm {read_files}')


rule filter_hits:
    input:
        sam_ct = join(BAM_DIR, '{sample}.raw.CT.sam'),
        sam_ga = join(BAM_DIR, '{sample}.raw.GA.sam')
    output:
        temp(join(BAM_DIR, '{sample}.filtered.sam'))
    threads: 1
    run:
        shell('samtools view -h -q {min_map_score} {input.sam_ct} > {output}')
        shell('samtools view -q {min_map_score} {input.sam_ga} >> {output}')


rule process_hits:
    input:
        fai = GENOME + '.fai',
        sam = join(BAM_DIR, '{sample}.filtered.sam')
    output:
        temp(join(BAM_DIR, '{sample}.processed.sam'))
    log:
        join('log', '{sample}.process_hits.log')
    threads: 1
    run:
        shell('sort -k1,1rd -T temp -S 20G {input.sam} | bin/process_sam_hits.pl {input.fai} {output} > {log} 2>&1')


rule fix_mate_info:
    input:
        fai = GENOME + '.fai',
        sam = join(BAM_DIR, '{sample}.processed.sam')
    output:
        temp(join(BAM_DIR, '{sample}.fixed.mates.bam'))
    threads: 1
    run:
        shell('bin/fixMateInfo.pl < {input.sam} | samtools view -bSt {input.fai} - | samtools sort -o {output} - ')


rule clip_overlap:
    input:
        join(BAM_DIR, '{sample}.fixed.mates.bam')
    output:
        temp(join(BAM_DIR, '{sample}.clipped.sorted.bam'))
    log:
        join('log', '{sample}.clip_overlap.log')
    threads: 1
    run:
        shell('bam clipOverlap --stats --poolSize 4000000 --in {input} --out {output} > {log} 2>&1')


rule rmdup:
    input:
        fai = GENOME + '.fai',
        bam = join(BAM_DIR, '{sample}.clipped.sorted.bam')
    output:
        temp(join(BAM_DIR, '{sample}.clipped.sorted.rmdup.bam'))
    log:
        join('log', '{sample}.rmdup.log')
    threads: 1
    run:
        if UMI_LENGTH > 0:
            shell('samtools view {input.bam} | bin/sam_UMI_filter_Poisson.pl | samtools view -bSt {input.fai} - > {output}')
        else:
            shell('samtools rmdup -S {input.bam} {output} > {log} 2>&1')


rule make_hapinfo:
    input:
        cpg = GENOME + '.cpg.positions.txt',
        bam = join(BAM_DIR, '{sample}.clipped.sorted.rmdup.bam')
    output:
        protected(join(METHYL_DIR, '{sample}.hapinfo.txt'))
    log:
        join('log', '{sample}.get_hapinfo.log')
    threads: 1
    run:
        shell('bin/getHaplo_PE_cgOnly.pl {input.bam} {input.cpg} > {output} 2> {log}')


rule mapping_report:
    input:
        e = expand(join(METHYL_DIR, '{sample}.hapinfo.txt'), sample = SAMPLES)
    output:
        protected('mapping_report.txt')
    run:
        shell('echo {input} > {output}')


rule preprocess_report:
    input:
        e = expand(join(FASTQ_DIR, '{sample}.cutadapt.lst'), sample = SAMPLES)
    output:
        protected('preprocess_report.txt')
    run:
        b = lambda x: bytes(x, 'UTF8')

        with open(output[0], 'wb') as out:
            out.write(b("""
Bisulfite-seq fastq preprocessing workflow\n
==========================================\n\n

R1 = represents the first read\n
R2 = represents the second read of paired-end/mate-pairs\n

Reads were decompressed using UNIX less command.\n

Then if UMI length was specified, the corresponding base pairs were trimmed \n
from the start of R1 to represent molecular barcodes (UMI). \n
        
Trimming of adapters from the 3' end was performed using cutadapt using the specified settings:\n

for R1: {}\n
for R2: {}\n\n

The reads were trimmed at 5' using the trimming mode: {}.\n

All Cytosine calls were masked and the original reads was appended to the read IDs.\n
""".format(cutadapt1, cutadapt2, TRIM_MODE)))

            out.write(b("""
The following are starting reads count and trimmed reads count for each sample:\n"""))
            for i in input.e:
                read_files = [line.rstrip('\n') for line in open(i)]
                for f in read_files:
                    n_written = 0
                    n_reads = 0
                    with open(f.split()[1]) as log:
                        for l in log:
                            if 'Total reads' in l:
                                n_reads = l.split()[3].rstrip('\n')
                            if 'Reads written' in l:
                                n_written = l.split()[4].rstrip('\n')
                    out.write(b('{} {} {} {} {}\n'.format(i, f.split()[1], n_reads, n_written, float(n_written.replace(",", ""))/float(n_reads.replace(",", "")))))


