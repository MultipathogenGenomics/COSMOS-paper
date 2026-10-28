
rule trim:
    input:
        R1 = lambda wc: config['sample_sets'][wc.sample_set]['samples'][wc.sample]['R1'],
        R2 = lambda wc: config['sample_sets'][wc.sample_set]['samples'][wc.sample]['R2'],
        adapters = lambda wc: config["sample_sets"][wc.sample_set]['adapters']
    output:
        r1 = temp(config["scratch_dir"] + "/{sample_set}/pre-processed/trimmed/{sample}/{sample}_R1_trimmed.fastq"),
        r2 = temp(config["scratch_dir"] + "/{sample_set}/pre-processed/trimmed/{sample}/{sample}_R2_trimmed.fastq"),
        r1_unpaired = temp(config["scratch_dir"] + "/{sample_set}/pre-processed/trimmed/{sample}/{sample}_R1_unpaired.fastq"),
        r2_unpaired = temp(config["scratch_dir"] + "/{sample_set}/pre-processed/trimmed/{sample}/{sample}_R2_unpaired.fastq"),
        r1_gz = config["scratch_dir"] + "/{sample_set}/pre-processed/trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
        r2_gz = config["scratch_dir"] + "/{sample_set}/pre-processed/trimmed/{sample}/{sample}_R2_trimmed.fastq.gz"    
    resources:
        project = config["project"],
        queue = "normal",
        mem = "12G",
        walltime = "04:00:00",
        cpus = 4
    envmodules:
        "Trimmomatic"
    params:
        min_length = 50,
        output_dir = lambda wc, output: os.path.dirname(output.r1)
    shell:
        """
        mkdir -p {params.output_dir}
        trimmomatic PE -threads {resources.cpus} \
        {input.R1} {input.R2} \
        {output.r1} {output.r1_unpaired} \
        {output.r2} {output.r2_unpaired} \
        ILLUMINACLIP:{input.adapters}:2:10:7:1:true \
        MINLEN:{params.min_length}
        gzip -c {output.r1} > {output.r1_gz} 
        gzip -c {output.r2} > {output.r2_gz}
        """


