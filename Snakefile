wildcard_constraints:
    source = r"[MucusEV]+",
    seq_id = r"\d+",
    sample = r"\d+"

rule all:
    input: 
        expand(
            "{source}_HZ_{seqid}_S{sample}_L001_unmapped.bracken", 
            zip, 
            source=["Mucus","Mucus","Mucus","Mucus","Mucus","Mucus","EV","EV","EV","EV","EV","EV"],
            seqid=[1,2,3,4,5,6,7,8,9,10,11,12],
            sample=[11,9,1,5,2,10,7,4,8,12,6,3]),
        expand(
            "{source}_HZ_{seqid}_S{sample}_L001_bhist.txt", 
            zip, 
            source=["Mucus","Mucus","Mucus","Mucus","Mucus","Mucus","EV","EV","EV","EV","EV","EV"],
            seqid=[1,2,3,4,5,6,7,8,9,10,11,12],
            sample=[11,9,1,5,2,10,7,4,8,12,6,3])
            #("Mucus", 1, 11),   ("Mucus", 2, 9),
            #("Mucus", 3, 1),    ("Mucus", 4, 5),
            #("Mucus", 5, 2),    ("Mucus", 6, 10),
            #("EV", 7, 7),       ("EV", 8, 4),
            #("EV", 9, 8),       ("EV", 10, 12),
            #("EV", 11, 6),      ("EV", 12, 3),

# ------------------------ Step 1: preprocess paired-end reads ------------------------

rule histogram_generation:
    input:
        in1="../input/{source}_HZ_{seq_id}_S{sample}_L001_R1_001.fastq.gz",
        in2="../input/{source}_HZ_{seq_id}_S{sample}_L001_R2_001.fastq.gz"

    output:
        bhist="{source}_HZ_{seq_id}_S{sample}_L001_bhist.txt",
        qhist="{source}_HZ_{seq_id}_S{sample}_L001_qhist.txt"

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=highmem",

    shell: """bbduk.sh in1={input.in1} in2={input.in2} \\
            bhist={output.bhist} qhist={output.qhist}"""

rule interleave:
    input: 
        in1="../input/{source}_HZ_{seq_id}_S{sample}_L001_R1_001.fastq.gz",
        in2="../input/{source}_HZ_{seq_id}_S{sample}_L001_R2_001.fastq.gz"

    conda: "envs/ev_analysis.yaml"

    output: temp("{source}_HZ_{seq_id}_S{sample}_L001.fastq.gz")

    shell: """
        reformat.sh in1={input.in1} in2={input.in2} \
        out={output} trd addslash int
    """

rule quality_trimming:
    input: "{source}_HZ_{seq_id}_S{sample}_L001.fastq.gz"

    output: temp("{source}_HZ_{seq_id}_S{sample}_L001_clean.fastq.gz")

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=highmem",
        mem_mb = 128000*workflow.attempt

    shell: "bbduk.sh in={input} out={output} qtrim=r trimq=10"

rule deinterleave:
    input: "{source}_HZ_{seq_id}_S{sample}_L001_clean.fastq.gz"

    conda: "envs/ev_analysis.yaml"

    output: 
        out1 = temp("{source}_HZ_{seq_id}_S{sample}_L001_R1_001_cleaner.fastq.gz"),
        out2 = temp("{source}_HZ_{seq_id}_S{sample}_L001_R2_001_cleaner.fastq.gz")

    shell: """
        reformat.sh in={input} out1={output.out1} out2={output.out2}
    """

# -------------------------- Step 2: align to human reference -------------------------

rule get_human_ref_index:
    output: directory("GRCh38_noalt_as/")

    resources:
        slurm_extra = "--qos=huge-long"

    shell: """
            curl https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip -o human_index.zip
            unzip human_index.zip"""

rule align_trimmed_reads_to_index:
    input:
        in1="{source}_HZ_{seq_id}_S{sample}_L001_R1_001_cleaner.fastq.gz",
        in2="{source}_HZ_{seq_id}_S{sample}_L001_R2_001_cleaner.fastq.gz",
        index=rules.get_human_ref_index.output

    output: temp("{source}_HZ_{seq_id}_S{sample}_L001_aligned.sam"),

    resources:
        slurm_extra = "--qos=huge-long",
        runtime=480*workflow.attempt

    conda: "envs/ev_analysis.yaml"

    threads: 8

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "bowtie2 -x {input.index}/GRCh38_noalt_as -p 8 -1 {input.in1} -2 {input.in2} -S {output}"

# ----------------------------- Step 3: get unmapped reads ----------------------------

rule from_sam_to_bam:
    input: "{source}_HZ_{seq_id}_S{sample}_L001_aligned.sam"

    output: temp("{source}_HZ_{seq_id}_S{sample}_L001_aligned.bam")

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "samtools view -bS {input} > {output}"

rule from_bam_to_sorted:
    input: "{source}_HZ_{seq_id}_S{sample}_L001_aligned.bam"

    output: "{source}_HZ_{seq_id}_S{sample}_L001_aligned.sorted.bam"

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "samtools sort -n {input} -o {output}"

rule output_unmapped_to_fastq:
    input: "{source}_HZ_{seq_id}_S{sample}_L001_aligned.sorted.bam"

    output:
        out1="{source}_HZ_{seq_id}_S{sample}_L001_R1_001_unmapped.fastq",
        out2="{source}_HZ_{seq_id}_S{sample}_L001_R2_001_unmapped.fastq"

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "samtools fastq -1 {output.out1} -2 {output.out2} --rf 0x4 {input}"

# -------------------------- Step 4: classify unmapped reads --------------------------

rule create_std_kraken2_database:
    output: directory("k2db/taxonomy")

    threads: 16

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "curl -s -L https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240112.tar.gz \
    | tar xvz -C k2db"

rule classify_unmapped_reads:
    input:
        in1="{source}_HZ_{seq_id}_S{sample}_L001_R1_001_unmapped.fastq",
        in2="{source}_HZ_{seq_id}_S{sample}_L001_R2_001_unmapped.fastq"

    output: "{source}_HZ_{seq_id}_S{sample}_L001_unmapped.kreport"

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=huge-long"

    shell: """kraken2 --db k2db --report {output}\
            --paired {input.in1} {input.in2}"""

# -------------------------- Step 5: validate classification --------------------------

rule divide_database:
    input: "k2db/taxonomy"

    conda: "envs/ev_analysis.yaml"

    output: "k2db/database100mers.kmer_distrib"

    threads: 16

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "bracken-build -d k2db -t {threads} -k 35 -l 100"

rule abundance_estimation:
    input:
        bracken_build="k2db/database100mers.kmer_distrib",
        report="{source}_HZ_{seq_id}_S{sample}_L001_unmapped.kreport"

    output: "{source}_HZ_{seq_id}_S{sample}_L001_unmapped.bracken"

    conda: "envs/ev_analysis.yaml"

    resources:
        slurm_extra = "--qos=huge-long"

    shell: "bracken -d k2db -i {input.report} -o {output}"
