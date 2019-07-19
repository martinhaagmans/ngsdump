REF = "/home/mahaagmans/Documents/referentie/hg19.fa"
REFMM2 = "/home/mahaagmans/Documents/referentie/hg19.omni"
FAST5 = "/media/sf_E_DRIVE/Nanopore/Methylering/1_2_3/20190514_1219_MN29521_FAK43976_42bf5829/fast5_pass"

INTERVALS = [
    "chr16:86541923-86543923",
    "chr6:2841139-2843139",
    "chr12:14925986-14927986",
    "chr1:97990020-97992020",
    "chr1:38226521-38228521",
    "chr9:96622674-96624674",
    "chr5:148519669-148521669"
             ]

SAMPLES = ["BC13", "BC14", "BC15"]

rule all:
    input:
        expand("{sample}.sorted.bam", sample=SAMPLES),
        expand("{sample}.fastq.index.readdb", sample=SAMPLES),
        expand("{sample}.methylation_called.txt", sample=SAMPLES),

rule preprocessing:
    input:
        fast5 = FAST5,
        fastq = "{sample}.fastq"
    output:    
        index = "{sample}.fastq.index",
        fai = "{sample}.fastq.index.fai",
        gzi = "{sample}.fastq.index.gzi",
        readdb = "{sample}.fastq.index.readdb" 
    shell:
        "nanopolish index -d {input.fast5} {input.fastq}"

rule mapreads:
    input:
        fastq = "{sample}.fastq"
    output:
        bam = "{sample}.sorted.bam"
    shell:
        "minimap2 -a -x map-ont {REFMM2} {input.fastq} | samtools sort -T tmp -o {output.bam} && samtools index {output.bam}"

rule call_methylation:
    input:
        fastq = "{sample}.fastq",
        bam = rules.mapreads.output.bam,
        readdb = rules.preprocessing.output.readdb,
         
    output:
        "{sample}.methylation_called.txt"

    run:
        for interval in INTERVALS:
            interval_name = interval.replace(':', '-')
            interval_name = interval_name.replace('-', '_')
            shell("""nanopolish call-methylation -r {input.fastq} -b {input.bam} -g {REF} -w "{interval}" > {wildcards.sample}_{interval_name}_methylation_calls.tsv""")
        shell("touch {output}")



# scripts/calculate_methylation_frequency.py -i methylation_calls.tsv > methylation_frequency.tsv

