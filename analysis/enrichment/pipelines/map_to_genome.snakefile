
rule map_reads_to_genome:
    input:
        # R1 = lambda wc: config["sample_sets"][wc.sample_set]["samples"][wc.sample]["R1"],
        # R2 = lambda wc: config["sample_sets"][wc.sample_set]["samples"][wc.sample]["R2"],
        R1 = lambda wc: determine_abundance_pipeline_input_files(wc, config)[0],
        R2 = lambda wc: determine_abundance_pipeline_input_files(wc, config)[1],
        genome_ref = lambda wc: config['genomes'][wc.genome]['reference']
    output:
        sorted_bam = config['scratch_dir'] + "/{sample_set}/aligned/{genome}/{sample}/{sample}_{fastq_type}_aligned_to_genome.sorted.bam",
        mapped_r1 = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/{genome}/{sample}/{sample}_{fastq_type}_aligned_to_{genome}_R1.fastq"),
        mapped_r2 = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/{genome}/{sample}/{sample}_{fastq_type}_aligned_to_{genome}_R2.fastq"),
        mapped_r1_gz = config['scratch_dir'] + "/{sample_set}/pre-processed/{genome}/{sample}/{sample}_{fastq_type}_aligned_to_{genome}_R1.fastq.gz",
        mapped_r2_gz = config['scratch_dir'] + "/{sample_set}/pre-processed/{genome}/{sample}/{sample}_{fastq_type}_aligned_to_{genome}_R2.fastq.gz",
    envmodules:
        "minimap2",
        "samtools",
    resources:
        project = config["project"],
        queue = "normal",
        mem = "24G",
        walltime = "04:00:00",
        cpus = 4
    shell:
        """
        minimap2 -t {resources.cpus} -ax sr --secondary=yes {input.genome_ref} {input.R1} {input.R2} | samtools view -bS -F 4 - | samtools sort -o {output.sorted_bam}

        samtools view -b {output.sorted_bam} | \
        samtools fastq -@ {resources.cpus} -1 {output.mapped_r1} -2 {output.mapped_r2} -0 /dev/null -s /dev/null -
        gzip -c {output.mapped_r1} > {output.mapped_r1_gz}
        gzip -c {output.mapped_r2} > {output.mapped_r2_gz}
        samtools index {output.sorted_bam}
        """

rule compute_coverage_across_genome:
    input:
        bam = rules.map_reads_to_genome.output.sorted_bam,
    output:
        bai = config['scratch_dir'] + "/{sample_set}/aligned/{genome}/{sample}/{sample}_{fastq_type}_aligned_to_genome.sorted.bam.bai",
        coverage = config['scratch_dir'] + "/{sample_set}/aligned/{genome}/{sample}/{sample}_{fastq_type}_coverage.tsv.gz",
    envmodules:
        "bedtools",
        "samtools",
    resources:
        project = config["project"],
        mem = "4G",
        walltime = "02:00:00",
        cpus = 1,
    shell:
        """
        samtools index {input.bam}
        bedtools genomecov -ibam {input.bam} -d | gzip > {output.coverage}
        """

rule find_kma_assigned_reads:
    input:
        R1 = lambda wc: config["sample_sets"][wc.sample_set]["samples"][wc.sample]["R1"],
        R2 = lambda wc: config["sample_sets"][wc.sample_set]["samples"][wc.sample]["R2"],
        kma_output = lambda wc: config["scratch_dir"] + f"/{wc.sample_set}/abundance/{wc.target_set}/kma/{wc.sample}.zip"
    output:
        read_list = temp(config["scratch_dir"] + "/{sample_set}/pre-processed/kma_assigned/{target_set}/{genome}/{sample}/{sample}_reads_mapped_to_{genome}_kma_targets.txt"),
        R1_filtered = config['scratch_dir'] + "/{sample_set}/pre-processed/kma_assigned/{target_set}/{genome}/{sample}/{sample}_assigned_to_kma_targets_R1.fastq.gz",
        R2_filtered = config['scratch_dir'] + "/{sample_set}/pre-processed/kma_assigned/{target_set}/{genome}/{sample}/{sample}_assigned_to_kma_targets_R2.fastq.gz",
    envmodules:
        "seqtk",
    resources:
        project = config["project"],
        queue = "normal",
        mem = "4G",
        walltime = "04:00:00",
        cpus = 1
    params:
        kma_location = lambda wc: config["scratch_dir"] + f"/{wc.sample_set}/abundance/{wc.target_set}/kma/{wc.sample}",
        kma_frag_file  = lambda wc: config["scratch_dir"] + f"/{wc.sample_set}/abundance/{wc.target_set}/kma/{wc.sample}/{wc.sample}_kma.frag.gz",
        targets = lambda wc: config["genomes"][wc.genome]['targets'],
        script_dir = config['snakemake_srcdir'],
        script_module = "pipelines.scripts.map_to_genome.extract_reads_mapped_to_targets_in_kma",
        venv = config['venv'],
    shell:
        """
        # unzip kma output and extract mapped read names
        unzip -d {params.kma_location} {input.kma_output}
        # move the files in sub_dirs in kma location up to kma location for easier access
        find {params.kma_location} -mindepth 2 -type f -exec mv {{}} {params.kma_location} \;
        # get read ids to those aligned to targets in frag file
        source {params.venv}
        export PYTHONPATH={params.script_dir}:$PYTHONPATH
        python -m {params.script_module} --frag_file {params.kma_frag_file} --targets {params.targets} --read_list {output.read_list}
        # extract reads from fastqs
        seqtk subseq {input.R1} {output.read_list} | gzip -c > {output.R1_filtered}
        seqtk subseq {input.R2} {output.read_list} | gzip -c > {output.R2_filtered}
        # clean up unzipped
        rm -r {params.kma_location}
        """


rule map_kma_assigned_reads_to_genome:
    input:
        R1 = rules.find_kma_assigned_reads.output.R1_filtered,
        R2 = rules.find_kma_assigned_reads.output.R2_filtered,
        genome_ref = lambda wc: config['genomes'][wc.genome]['reference'],
    output:
        sorted_bam = config['scratch_dir'] + "/{sample_set}/align_kma_assigned/{target_set}/{genome}/{sample}/{sample}_kma_mapped_to_genome.sorted.bam",
        sorted_bai = config['scratch_dir'] + "/{sample_set}/align_kma_assigned/{target_set}/{genome}/{sample}/{sample}_kma_mapped_to_genome.sorted.bam.bai",
    envmodules:
        "minimap2",
        "samtools",
    resources:
        project = config["project"],
        queue = "normal",
        mem = "4G",
        walltime = "04:00:00",  
        cpus = 1
    params:
        targets = lambda wc: config["genomes"][wc.genome]['targets'],
    shell:
        """
        # align
        if [ ! -s {input.R1} ] || [ ! -s {input.R2} ]; then
            echo "No reads after KMA filtering; writing empty BAM"
            samtools view -H {input.genome_ref} | samtools view -bS -o {output.sorted_bam} -
            exit 0
        else
            minimap2 -ax sr --secondary=yes {input.genome_ref} {input.R1} {input.R2} | samtools view -bS -F 4 - | samtools sort -o {output.sorted_bam}
        fi
        samtools index {output.sorted_bam}
        """
