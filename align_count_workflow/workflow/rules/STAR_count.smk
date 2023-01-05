rule STAR:
    input:
        read1 = datadir_fq + "dedub.{sample}_R1_001.fastq.gz.gz",
        read2 = datadir_fq + "dedub.{sample}_R2_001.fastq.gz.gz"
    output:
        "unsorted_bam/{sample}Aligned.out.bam"
    threads:
        24
    message: "Aligning {input} using STAR: {threads} threads"
    params:
        genomedir = STARINDEX,
        gtf = config['ref']['gtf_anno'],
        outprefix = "{sample}",
        starlogs = 'mapped/starlogs' 
    shell:
        """
        STAR --genomeDir {params.genomedir} --sjdbGTFfile {params.gtf} \
        --readFilesIn {input.read1} {input.read2} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --chimSegmentMin 20 \
        --outSAMattributes All \
        --outSAMtype BAM Unsorted \
        --limitBAMsortRAM 10000000000 \
        --limitOutSJcollapsed 10000000 \
        --limitIObufferSize 500000000 \
        --quantMode GeneCounts \
        --outSAMunmapped Within KeepPairs \
        --runThreadN {threads} \
        --outFileNamePrefix {params.outprefix}\
        && mkdir -p {params.starlogs} && mv {params.outprefix}Log.final.out \
        {params.outprefix}Log.out {params.outprefix}Log.progress.out {params.starlogs}
        """

rule sort_by_name:
    input:
        "unsorted_bam/{sample}Aligned.out.bam"
    output: 
        output_directory + "/name_sorted/{sample}.sortedByName.bam"
    threads: 4
    params:
        tmp_dir = config['tmp_dir']
    shell:
        """
        samtools sort -on {input} -@ {threads} -T {params.tmp_dir} -o {output}
        """

rule featurecounts:
    input:
        output_directory + "/name_sorted/{sample}.sortedByName.bam"
    output:
        output_directory + "/featurecounts/{sample}.exon.counts.tsv"
    threads: 16
    params:
        anno = config['ref']['gtf_anno']
    shell:
        """
        featureCounts -t exon -g gene_id -a {params.anno} -T {threads} -o {output} {input} 
        """

