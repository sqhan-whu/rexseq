#!/usr/bin/env runsnakemake
"""
REXseq: RNA Editing Detection Pipeline
A comprehensive Snakemake pipeline for detecting RNA editing sites from RNA-seq data.
"""

configfile: "config.yaml"

# Import configuration
samples = config["samples"]
threads = config["threads"]
path = config["path"]
py_script = config["py_script"]
BED = config["BED"]
genome_fa = config["hg38_fa"]
genome_size = config["hg38_size"]
hisat2_index = config["hg38_hisat2_index"]
ANNOVER = config["ANNOVER"]
anno_db = config["ANNOVER_hg38_DATABASE"]
Alu_bed = config["Alu_hg38_bed"]
trans_id = config["trans_id"]

# Define directories
indir = path
outdir = "results"
bam_dir = outdir + "/bam"
hit_dir = bam_dir + "/hit"
seperate_dir = hit_dir + "/seperate"
sites_dir = outdir + "/final_sites"
anno_dir = outdir + "/annover"

# Parameters for trimming
seqtk_fq1_b = 3
seqtk_fq1_e = 10
seqtk_fq2_b = 10
seqtk_fq2_e = 3

rule all:
    input:
        expand(sites_dir + "/{sample}.annotation.txt", sample=samples),
        expand(anno_dir + "/{sample}.variant_function", sample=samples)

# ==================== STEP 1: Read Mapping ====================
rule trim_galore:
    input:
        indir + "/fq/{sample}_R1.fq.gz",
        indir + "/fq/{sample}_R2.fq.gz"
    output:
        outdir + "/clean_fq/{sample}_R1_val_1.fq",
        outdir + "/clean_fq/{sample}_R2_val_2.fq"
    params:
        outdir + "/clean_fq"
    threads: threads
    shell:
        "trim_galore -j {threads} --quality 20 --phred33 --stringency 3 --length 30 --dont_gzip -o {params} --paired {input}"

rule fastuniq:
    input:
        fq1 = outdir + "/clean_fq/{sample}_R1_val_1.fq",
        fq2 = outdir + "/clean_fq/{sample}_R2_val_2.fq"
    output:
        fq1 = temp(outdir + "/clean_fq/{sample}_uniq_1.fq"),
        fq2 = temp(outdir + "/clean_fq/{sample}_uniq_2.fq")
    shell:
        "perl {py_script}/readCollapse.pl -1 {input.fq1} -2 {input.fq2} -p {output.fq1} -q {output.fq2}"

rule seqtk_trimfq1:
    input:
        outdir + "/clean_fq/{sample}_uniq_1.fq"
    output:
        outdir + "/clean_fq/{sample}_seqtk_1.fq"
    shell:
        "seqtk trimfq -b {seqtk_fq1_b} -e {seqtk_fq1_e} {input} > {output}"

rule seqtk_trimfq2:
    input:
        outdir + "/clean_fq/{sample}_uniq_2.fq"
    output:
        outdir + "/clean_fq/{sample}_seqtk_2.fq"
    shell:
        "seqtk trimfq -b {seqtk_fq2_b} -e {seqtk_fq2_e} {input} > {output}"

rule hisat2_mapping:
    input:
        fq1 = outdir + "/clean_fq/{sample}_seqtk_1.fq",
        fq2 = outdir + "/clean_fq/{sample}_seqtk_2.fq"
    output:
        sam = temp(bam_dir + "/{sample}.sam"),
        log = bam_dir + "/{sample}.mapping_rate.txt"
    threads: threads
    params:
        "-I 0 -X 1000 --ma 0 --mp 4,2 --score-min L,0,-0.3 --rna-strandness RF --no-softclip --max-intronlen 4000 --dta-cufflinks --no-mixed --no-discordant --no-unal --no-temp-splicesite"
    shell:
        "hisat2 -p {threads} {params} -x {hisat2_index} -1 {input.fq1} -2 {input.fq2} -S {output.sam} 2> {output.log}"

rule samtools_sam2bam:
    input:
        bam_dir + "/{sample}.sam"
    output:
        temp(bam_dir + "/{sample}.ori.bam")
    params:
        bam_dir + "/{sample}.tmp"
    threads: threads
    shell:
        "samtools view -q 5 -F 2308 -@ {threads} -bS {input} | samtools sort -@ {threads} -m 1000M -o {output} -T {params}"

rule filter_bam:
    input:
        bam_dir + "/{sample}.ori.bam"
    output:
        bam_dir + "/{sample}.sorted.filter.bam"
    shell:
        "bamtools filter -in {input} -out {output} -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'"

rule bam_index:
    input:
        bam_dir + "/{sample}.sorted.filter.bam"
    output:
        bam_dir + "/{sample}.sorted.filter.bam.bai"
    threads: threads
    shell:
        "samtools index -@ {threads} {input}"

# ==================== STEP 2: Gene Hitting ====================
rule bedtools_intersect_1:
    input:
        bam_dir + "/{sample}.sorted.filter.bam"
    output:
        temp(hit_dir + "/{sample}.mate1.bam")
    shell:
        "bamtools filter -in {input} -out {output} -isFirstMate true"

rule bedtools_intersect_2:
    input:
        bam_dir + "/{sample}.sorted.filter.bam"
    output:
        temp(hit_dir + "/{sample}.mate2.bam")
    shell:
        "bamtools filter -in {input} -out {output} -isSecondMate true"

rule intersect_genebed_mate1:
    input:
        hit_dir + "/{sample}.mate1.bam"
    output:
        temp(hit_dir + "/{sample}.mate1.hit.bam")
    params: BED
    shell:
        "bedtools intersect -a {input} -b {params} -u -S -f 1.0 > {output}"

rule intersect_genebed_mate2:
    input:
        hit_dir + "/{sample}.mate2.bam"
    output:
        temp(hit_dir + "/{sample}.mate2.hit.bam")
    params: BED
    shell:
        "bedtools intersect -a {input} -b {params} -u -s -f 1.0 > {output}"

rule merge_hitbam:
    input:
        mate1 = hit_dir + "/{sample}.mate1.hit.bam",
        mate2 = hit_dir + "/{sample}.mate2.hit.bam"
    output:
        hit_dir + "/{sample}.hit.bam"
    shell:
        "bamtools merge -in {input.mate1} -in {input.mate2} -out {output}"

# ==================== STEP 3: Strand Separation ====================
rule positive_strand:
    input:
        hit_dir + "/{sample}.hit.bam"
    output:
        temp(seperate_dir + "/{sample}.pos.sam")
    shell:
        "samtools view {input} | grep 'XS:A:+' > {output}"

rule negative_strand:
    input:
        hit_dir + "/{sample}.hit.bam"
    output:
        temp(seperate_dir + "/{sample}.neg.sam")
    shell:
        "samtools view {input} | grep 'XS:A:-' > {output}"

rule get_head:
    input:
        hit_dir + "/{sample}.hit.bam"
    output:
        seperate_dir + "/{sample}.header"
    shell:
        "samtools view -H {input} > {output}"

rule merge_head_pos_sam:
    input:
        seperate_dir + "/{sample}.header",
        seperate_dir + "/{sample}.pos.sam"
    output:
        temp(seperate_dir + "/{sample}.pos.0.sam")
    shell:
        "cat {input} > {output}"

rule merge_head_neg_sam:
    input:
        seperate_dir + "/{sample}.header",
        seperate_dir + "/{sample}.neg.sam"
    output:
        temp(seperate_dir + "/{sample}.neg.0.sam")
    shell:
        "cat {input} > {output}"

rule samtools_sam2bam_pos:
    input:
        seperate_dir + "/{sample}.pos.0.sam"
    output:
        temp(seperate_dir + "/{sample}.pos.bam")
    shell:
        "samtools view -bS {input} > {output}"

rule samtools_sam2bam_neg:
    input:
        seperate_dir + "/{sample}.neg.0.sam"
    output:
        temp(seperate_dir + "/{sample}.neg.bam")
    shell:
        "samtools view -bS {input} > {output}"

rule bamtools_sort_pos:
    input:
        seperate_dir + "/{sample}.pos.bam"
    output:
        seperate_dir + "/{sample}.pos.sort.bam"
    shell:
        "bamtools sort -in {input} -out {output}"

rule bamtools_sort_neg:
    input:
        seperate_dir + "/{sample}.neg.bam"
    output:
        seperate_dir + "/{sample}.neg.sort.bam"
    shell:
        "bamtools sort -in {input} -out {output}"

rule samtools_index_pos:
    input:
        seperate_dir + "/{sample}.pos.sort.bam"
    output:
        seperate_dir + "/{sample}.pos.sort.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_index_neg:
    input:
        seperate_dir + "/{sample}.neg.sort.bam"
    output:
        seperate_dir + "/{sample}.neg.sort.bam.bai"
    shell:
        "samtools index {input}"

# ==================== STEP 4: Site Detection ====================
rule pos_strand_sites:
    input:
        seperate_dir + "/{sample}.pos.sort.bam",
        seperate_dir + "/{sample}.pos.sort.bam.bai"
    output:
        sites_dir + "/{sample}.pos.rough.txt"
    params:
        script = py_script + "/AG.py",
        new_bam = sites_dir + "/{sample}.AG.bam"
    shell:
        "python {params.script} {input[0]} {genome_fa} {genome_size} {params.new_bam} > {output}"

rule neg_strand_sites:
    input:
        seperate_dir + "/{sample}.neg.sort.bam",
        seperate_dir + "/{sample}.neg.sort.bam.bai"
    output:
        sites_dir + "/{sample}.neg.rough.txt"
    params:
        script = py_script + "/TC.py",
        new_bam = sites_dir + "/{sample}.TC.bam"
    shell:
        "python {params.script} {input[0]} {genome_fa} {genome_size} {params.new_bam} > {output}"

rule trans_neg2pos:
    input:
        sites_dir + "/{sample}.pos.rough.txt",
        sites_dir + "/{sample}.neg.rough.txt"
    output:
        sites_dir + "/{sample}.process.site.txt"
    params:
        py_script + "/trans_neg2pos.py"
    shell:
        "python {params} {input} > {output}"

rule stat_all_site:
    input:
        sites_dir + "/{sample}.process.site.txt"
    output:
        sites_dir + "/{sample}.process.site.txt_high_sites.txt",
        sites_dir + "/{sample}.process.site.txt_low_sites.txt",
        sites_dir + "/{sample}.process.site.txt_total.sites.txt"
    shell:
        "python {py_script}/find_error2.py {input}"

rule final_sites:
    input:
        sites_dir + "/{sample}.process.site.txt_total.sites.txt"
    output:
        sites_dir + "/{sample}.site.txt"
    params:
        py_script + "/error2sites.py"
    shell:
        "python {params} {input} > {output}"

# ==================== STEP 5: Annotation ====================
rule pre_annover:
    input:
        sites_dir + "/{sample}.site.txt"
    output:
        temp(anno_dir + "/{sample}.tmp.txt")
    params:
        py_script + "/pre_annover.py"
    shell:
        "python {params} {input} {output}"

rule annover:
    input:
        anno_dir + "/{sample}.tmp.txt"
    output:
        anno_dir + "/{sample}.variant_function"
    params:
        anno_db = anno_db,
        build = "hg38",
        annover = ANNOVER,
        name = anno_dir + "/{sample}"
    shell:
        "perl {params.annover} -out {params.name} -build {params.build} -dbtype ensGene {input} {params.anno_db}"

rule annotation:
    input:
        sites_dir + "/{sample}.site.txt",
        anno_dir + "/{sample}.variant_function",
        Alu_bed,
        genome_fa,
        trans_id
    output:
        sites_dir + "/{sample}.annotation.txt"
    params:
        py_script + "/sites_table2.py"
    shell:
        "python {params} {input} {output}"