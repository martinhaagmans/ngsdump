PICARD = "/home/mahaagmans/programs/jre1.8.0_144/bin/java -jar /home/mahaagmans/programs/picard.jar"
REF = "/home/mahaagmans/Documents/referentie/hg19.fa"
TI = "HyperWES_target.interval_list"
BI = "HyperWES_tiled.interval_list"

SAMPLES = ["18D8473", "18D12942", "CEPH1", "CEPH2", "CEPH3", "NA12878", "NA12891", "NA12892"]


rule all:
    input:
        expand("{sample}.DM.bam", sample=SAMPLES),
        expand("output/{sample}.HSMetrics.txt", sample=SAMPLES),
        expand("output/{sample}.InsertSize.txt", sample=SAMPLES),
        expand("output/{sample}.InsertSize.pdf", sample=SAMPLES),
        expand("output/{sample}.AlignmentMetrics.txt", sample=SAMPLES)


rule mapreads:
    input:
        "reads/{sample}.R1.fastq.gz",
        "reads/{sample}.R2.fastq.gz"
    output:
        temp("output/{sample}.sorted.bam")
    threads:
        2
    log:
        "logfiles/{sample}.BWAlignment.log"
    message:
        "Aligning reads with bwa mem"
    params:
        rg = "@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tPU:{sample}\\tSM:{sample}"
    shell:
        '''(bwa mem  -R '{params.rg}' -t 1 -M {REF} {input} |\
        samtools view -Shu - |\
        samtools sort -T {wildcards.sample}.tmp -O bam - > {output}) > {log}  2>&1
        '''


rule markduplicates:
    input:
        rules.mapreads.output
    output:
        bam = "output/{sample}.DM.bam",
        metrics = "output/{sample}.dupmark.txt"
    log:
        "logfiles/{sample}.MarkDuplicates.log"
    message:
        "Marking duplicates with picard"
    shell:
        '''{PICARD} MarkDuplicates  I={input} O={output.bam}  \
        M={output.metrics} REMOVE_DUPLICATES=FALSE CREATE_INDEX=true > {log}  2>&1
        '''


rule hsmetrics:
    input:
        rules.markduplicates.output.bam
    output:
        "output/{sample}.HSMetrics.txt"
    log:
        "logfiles/{sample}.HSmetrics.log"
    message:
        "Calculating HS-metrics with picard"
    run:
        shell('''{PICARD} CalculateHsMetrics R={REF} I={input} \
        O={output[0]} TI={TI} BI={BI} > {log} 2>&1
        ''')


rule insertsizemetrics:
    input:
        rules.markduplicates.output.bam
    output:
        tabel = "output/{sample}.InsertSize.txt",
        pdf = "output/{sample}.InsertSize.pdf"
    log:
        "logfiles/{sample}.InsertSize.log"
    message:
        "Calculating insertsize metrics"
    run:
        shell('''{PICARD} CollectInsertSizeMetrics R={REF} I={input} \
        O={output.tabel} HISTOGRAM_FILE={output.pdf} > {log} 2>&1
        ''')


rule alignmetrics:
    input:
        rules.markduplicates.output.bam
    output:
        "output/{sample}.AlignmentMetrics.txt"
    log:
        "logfiles/{sample}.AlignmentMetrics.log"
    message:
        "Calculating alignment metrics"
    run:
        shell('''{PICARD} CollectAlignmentSummaryMetrics R={REF} I={input} \
        O={output[0]}  MAX_INSERT_SIZE=500 > {log} 2>&1
        ''')


