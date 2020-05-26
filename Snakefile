"""
Process SNARE-seq ATAC data
Requires bwa-mem, samtools, sinto, pigz, bgzip, and tabix
"""

IDS, = glob_wildcards("replicates/{rep}.txt")


rule all:
    input: "mapped/fragments.sort.bed.gz.tbi"

rule get_genome:
    """Download genome and build bwa index"""
    output:
        "genome/mm10.fa.gz"
    shell:
        """
        cd genome
        wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
        """

rule bwa_build:
    """Build index for bwa mem"""
    input:
        "genome/mm10.fa.gz"
    output:
        "genome/mm10.fa"
    shell:
        """
        gzip -d mm10.fa.gz
        bwa mem index {input} -p genome/mm10
        """

rule download_attach:
    """
    Download fastq files for each replicate
    Decompress fastq files
    Add barcodes to read 1 and read 2
    Compress fastqs (parallel using pigz)
    """
    input:
        expand("replicates/{rep}.txt", rep=IDS)
    output:
        "fastq/{wildcards.rep}/done.txt"
    threads: 8
    shell:
        """
        mkdir fastq
        cd fastq
        wget -i {input}
        gzip -S .gz.1 -d *
        sinto barcode -b 12 --barcode_fastq *R1_001.fastq --read1 *R2_001.fastq --read2 *R3_001.fastq --prefix {{rep}}_
        pigz -p {threads} *.fastq
        touch done.txt
        """

rule cat_fastq:
    """Concatenate fastq files from different reps"""
    input:
        "fastq/{wildcards.rep}/done.txt"
    output:
        "fastq/read1.fastq.gz"
    shell:
        """
        cat fastq/{{rep}}/*R2_001.barcoded.fastq.gz > fastq/read1.fastq.gz
        cat fastq/{{rep}}/*R3_001.barcoded.fastq.gz > fastq/read2.fastq.gz
        """

rule map_reads:
    """Map reads to genome"""
    input:
        "fastq/read1.fastq.gz",
        "genome/mm10.fa"
    output:
        "mapped/aln.bam"
    threads: 8
    shell:
        """
        bwa mem -t {threads} genome/mm10 fastq/read1.barcoded.fastq.gz fastq/read2.barcoded.fastq.gz \
            | samtools view -b - > {output}
        """

rule sort_bam:
    """Sort and index bam file"""
    input:
        "mapped/aln.bam"
    output:
        "mapped/aln.sort.bam.bai"
    threads: 8
    shell:
        """
        cd mapped
        samtools sort -@ {threads} aln.bam -o aln.sort.bam
        rm aln.bam
        samtools index -@ {threads} aln.sort.bam
        """

rule create_fragments:
    """Create fragment file from BAM file"""
    input:
        "mapped/aln.sort.bam.bai"
    output:
        "mapped/fragments.bed"
    threads: 8
    shell:
        """
        sinto fragments -b mapped/aln.bam -p {threads} -f {output} --barcode_regex "[^:]*"
        """

rule sort_fragments:
    """Sort and index fragment file"""
    input:
        "mapped/fragments.bed"
    output:
        "mapped/fragments.sort.bed.gz.tbi"
    threads: 8
    shell:
        """
        sort -k1,1 -k2,2n {input} > mapped/fragments.sort.bed
        bgzip -@ {threads} mapped/fragments.sort.bed
        tabix -p bed mapped/fragments.sort.bed.gz
        rm {input}
        """