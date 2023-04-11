import os
import pandas as pd
from collections import Counter
DIRECTION=["1","2"]
GROUPS=set()
sample_sheet=pd.read_csv("samples.tsv", sep="\t",dtype=object)
genome = '/mnt/d/2023/EPICS/EPIC_ITD_781_PRS/ITD_1030_PRS_Panel_quality/02_analyses/00_ref/human_g1k_v37_decoy.fasta'


def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

SAMPLES=list(sample_sheet["ID"])


for b in set(SAMPLES):
    os.makedirs("../01_raw/" +b+ "/fastqc", exist_ok=True)
    check_symlink(sample_sheet[sample_sheet["ID"]==b]["forward reads"].values[0], "../01_raw/"+b +"/"+ b +"_1P.fastq.gz")
    check_symlink(sample_sheet[sample_sheet["ID"]==b]["reverse reads"].values[0], "../01_raw/"+b +"/"+ b +"_2P.fastq.gz")
    os.makedirs("../02_trimmed/" + b +"/fastqc", exist_ok=True)
    os.makedirs("../03_mapped/"+b, exist_ok=True)


rule all:
    input:
        expand('../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html',direction=DIRECTION, sample=SAMPLES),
        expand("../03_mapped/{sample}/{sample}.mkdup.bam", sample=SAMPLES),

rule fastqc1:
    input:
        r = '../01_raw/{sample}/{sample}_{direction}P.fastq.gz',
    threads: 2
    priority: 50
    output:
        '../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html'
    conda:
        "envs/fastqc.yaml"
    shell:
        'fastqc -o ../01_raw/{wildcards.sample}/fastqc -t {threads} --extract {input.r}'

rule trimming:
    input:
        r1= '../01_raw/{sample}/{sample}_1P.fastq.gz',
        r2= '../01_raw/{sample}/{sample}_2P.fastq.gz',
    output:
        p1="../02_trimmed/{sample}/{sample}_trimmed_1P.fastq.gz",
        u1="../02_trimmed/{sample}/{sample}_trimmed_1U.fastq.gz",
        p2="../02_trimmed/{sample}/{sample}_trimmed_2P.fastq.gz",
        u2="../02_trimmed/{sample}/{sample}_trimmed_2U.fastq.gz",
    conda:
        "envs/trimmomatic.yaml"
    threads: 4
    priority: 50
    shell:
        'trimmomatic PE -threads {threads} \
        {input.r1} {input.r2} \
        {output.p1} {output.u1} \
        {output.p2} {output.u2} \
        ILLUMINACLIP:data/TruSeq3-PE.fa:2:30:10:2:true \
        LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36'

rule Minimap2:
    input:
        r1 = "../02_trimmed/{sample}/{sample}_trimmed_1P.fastq.gz",
        r2 = "../02_trimmed/{sample}/{sample}_trimmed_2P.fastq.gz",
        #r1= '../01_raw/{sample}/{sample}_1P.fastq.gz',
        #r2= '../01_raw/{sample}/{sample}_2P.fastq.gz',
        #sampleID="{sample}",
    output:
        ancient("../03_mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam"),
    conda:
        "envs/minimap.yaml"
    threads: 8
    priority: 50
    shell:
        "minimap2 -ax sr -t {threads} {genome} {input.r1} {input.r2} |\
        samtools view -b | samtools sort -o ../03_mapped/{wildcards.sample}/{wildcards.sample}Aligned.sortedByCoord.out.bam"


rule MarkDuplicates:
    input:
        bamOrig = ancient("../03_mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam")

    output:
        bamMkdup=  ancient("../03_mapped/{sample}/{sample}.mkdup.bam")

    threads: 2

    conda:
        "envs/samtools.yaml"

    shell:
        "samtools collate -o /dev/stdout {input.bamOrig} | samtools fixmate -m /dev/stdin /dev/stdout | samtools sort | samtools markdup /dev/stdin {output.bamMkdup}; samtools index {output.bamMkdup}"