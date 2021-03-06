"""
Process SNARE-seq ATAC data
Requires bwa-mem, samtools, sinto, bgzip, and tabix
"""

IND, = glob_wildcards("replicates/{rep}.txt")


rule all:
    input: "mapped/fragments.sort.bed.gz.tbi"

rule get_genome:
    output:
        "genome/mm10.fa.gz"
    threads: 1
    message: "Download genome and build bwa index"
    shell:
        """
        cd genome
        wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
        """

rule bwa_build:
    input:
        "genome/mm10.fa.gz"
    output:
        "genome/mm10.fa"
    threads: 1
    message: "Build index for bwa mem"
    shell:
        """
        gzip -d genome/mm10.fa.gz
        bwa index genome/mm10.fa
        """

rule download_attach:
    input:
        "replicates/{rep}.txt"
    output:
        "fastq/{rep}/done.txt"
    threads: 1
    message:
        """
        Download fastq files for each replicate
        Decompress fastq files
        Add barcodes to read 1 and read 2
        """
    shell:
        """
        wget -i {input} -P fastq/{wildcards.rep}
        cd fastq/{wildcards.rep}
        gzip -S .gz.1 -d *
        sinto barcode -b 12 --barcode_fastq *R1_001.fastq --read1 *R2_001.fastq --read2 *R3_001.fastq --prefix {wildcards.rep}_
        mv *.barcoded.fastq ..
        rm *.fastq
        touch done.txt
        """

rule cat_fastq:
    input:
        expand("fastq/{rep}/done.txt", rep=IND)
    output:
        "fastq/read1.fastq"
    threads: 1
    message: "Concatenate FASTQ files from different replicates"
    shell:
        """
        cat fastq/*R2_001.barcoded.fastq > fastq/read1.fastq
        cat fastq/*R3_001.barcoded.fastq > fastq/read2.fastq
        rm fastq/*R2_001.barcoded.fastq fastq/*R3_001.barcoded.fastq
        """

rule map_reads:
    input:
        read = "fastq/read1.fastq",
        idx = "genome/mm10.fa"
    output:
        "mapped/aln.bam"
    threads: 8
    message: "Map reads to genome"
    shell:
        """
        bwa mem -t {threads} {input.idx} fastq/read1.fastq fastq/read2.fastq \
            | samtools view -b - > {output}
        """

rule sort_bam:
    input:
        "mapped/aln.bam"
    output:
        "mapped/aln.sort.bam.bai"
    threads: 8
    message: "Sort and index bam file"
    shell:
        """
        cd mapped
        samtools sort -@ {threads} aln.bam -o aln.sort.bam
        rm aln.bam
        samtools index -@ {threads} aln.sort.bam
        """

rule create_fragments:
    input:
        "mapped/aln.sort.bam.bai"
    output:
        "mapped/fragments.bed"
    threads: 8
    message: "Create fragment file from BAM file"
    shell:
        """
        sinto fragments -b mapped/aln.sort.bam -p {threads} -f {output} --barcode_regex "[^:]*"
        """

rule sort_fragments:
    input:
        "mapped/fragments.bed"
    output:
        "mapped/fragments.sort.bed.gz.tbi"
    threads: 8
    message: "Sort and index fragment file"
    shell:
        """
        sort -k1,1 -k2,2n {input} > mapped/fragments.sort.bed
        bgzip -@ {threads} mapped/fragments.sort.bed
        tabix -p bed mapped/fragments.sort.bed.gz
        rm {input}
        """