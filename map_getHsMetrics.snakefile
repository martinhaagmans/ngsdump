REF = "/path/to/hg19.fa"
PICARD = "/path/to/java -jar /path/to/picard.jar"
SAMPLES = []


rule all:
    input:
        expand("output/{sample}.DM.bam", sample=SAMPLES),
        expand("output/{sample}.HSMetrics.txt", sample=SAMPLES)

rule mapreads:
    input:
        "{sample}.R1.fastq.gz",
        "{sample}.R2.fastq.gz"
    output:
        temp("tmp/{sample}.sorted.bam")
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
        metrics = temp("tmp/{sample}.dupmark.txt")
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
        bam = rules.markduplicates.output.bam,
        target = "/home/mahaagmans/Documents/ngstargets/captures/DLv4_target.interval_list"
    output:
        "output/{sample}.HSMetrics.txt"
    log:
        "logfiles/{sample}.HSmetrics.log"
    message:
        "Calculating HS-metrics with picard"
    run:
        shell('''{PICARD} CollectHsMetrics R={REF} I={input.bam} \
        O={output[0]} TI={input.target} BI={input.target} > {log} 2>&1
        ''')
