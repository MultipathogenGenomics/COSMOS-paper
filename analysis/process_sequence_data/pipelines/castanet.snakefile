import os

rule run_castanet:
    input:
        castanet_targets = lambda wc: config["sample_sets"][wc.sample_set]["castanet_targets"],
        reads_r1 = lambda wc: ancient(get_castanet_fastq_paths(wc, config)[0]),
        reads_r2 = lambda wc: ancient(get_castanet_fastq_paths(wc, config)[1]),
    output:
        # consenus_dir = directory(config["scratch_dir"] + "/{sample_set}/castanet/{sample}/consensus_data"),
        # depth_dir = directory(config["scratch_dir"] + "/{sample_set}/castanet/{sample}/Depth_output"),
        # grouped_reads_dir = directory(config["scratch_dir"] + "/{sample_set}/castanet/{sample}/grouped_reads"),
        # hashes_dir = directory(config["scratch_dir"] + "/{sample_set}/castanet/{sample}/hashes"),
        castanet_done = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/.castanet_done",
        bam_file = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/{sample}.bam",
    resources:
        project = config["project"],
        queue = "normal",
        walltime = lambda wc, input: assign_wt_ram_cpu_based_on_num_reads(wc, [input.reads_r1, input.reads_r2])[0],
        mem = lambda wc, input: assign_wt_ram_cpu_based_on_num_reads(wc, [input.reads_r1, input.reads_r2], boost=100)[1],
        cpus = lambda wc, input: assign_wt_ram_cpu_based_on_num_reads(wc, [input.reads_r1, input.reads_r2])[2] + 8,
        castanet_slots = 1,
        cluster_extra = " -ljobfs=20GB"
    envmodules:
        "Trimmomatic/0.39",
        "bwa-mem2/2.2.1",
        "kraken2/2.1.4",
        "viralconsensus/0.0.6",
        "samtools/1.19",
        "mafft/7.505"
    params:
        python_module = lambda wc: castanet[wc.version].get("python_module", ""), # can't use lambda in envmodules
        venv = lambda wc: castanet[wc.version]["venv"],
        pythonpath = lambda wc: castanet[wc.version]["python_path"],
        pythondir = lambda wc: castanet[wc.version]["python_dir"],
        post_filter = lambda wc: str(config["sample_sets"][wc.sample_set].get("post_filter", True)),
        reads_folder = lambda wc, input: os.path.dirname(input.reads_r1),
        sample = "{sample}",
        tmp_dir = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/tmp",
        output_tmp_dir_for_castanet = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/tmp/{sample}",
        output_final_dir = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}",
    shell:
        """
        set -euo pipefail

        # load python
        if [ -n "{params.python_module}" ]; then
            module load {params.python_module}
        fi

        # if a fiinal results folder exists make sure it's empty
        if [ -d "{params.output_final_dir}/{params.sample}" ]; then
            rm -rf "{params.output_final_dir}/{params.sample}/"*
        fi

        mkdir -p {params.tmp_dir}

        # ensure tmp dir gets removed on ANY exit (success or failure)
        tmpdir="{params.output_tmp_dir_for_castanet}"
        cleanup() {{
            echo "Cleaning up temporary directory: $tmpdir"
            rm -rf "$tmpdir"
        }}
        trap cleanup EXIT

        # run castanet
        source {params.venv}
        export PYTHONPATH={params.pythondir}:${{PYTHONPATH:-}}
        castanet_log={params.tmp_dir}/{params.sample}_castanet.log
        {params.pythonpath} -m app.castanet_lite \
            -ExpName {params.sample} \
            -ExpDir {params.reads_folder} \
            -SaveDir {params.output_tmp_dir_for_castanet} \
            -RefStem {input.castanet_targets} \
            -DoKrakenPrefilter False \
            -DoTrimming False \
            -DoConsensus True \
            -PostFilt {params.post_filter} \
            -NThreads {resources.cpus} \
            2>&1 | tee "$castanet_log"

        # # castanet can occasionally log fatal errors but still exit 0
        # if grep -Eq "ERROR:root:Traceback|SystemError:|Errcode: [1-9][0-9]*|^  Error:" "$castanet_log"; then
        #     echo "Castanet reported a fatal error in log; failing rule" >&2
        #     exit 1
        # fi

        # move final outputs into place
        mv -f "$tmpdir/{params.sample}" "{params.output_final_dir}"

        # remove *_filt.fastq if they exist as very large, but don't fail if they don't
        rm -f "{params.output_final_dir}"/{params.sample}/*_filt.fastq || true

        # # ensure expected consensus output exists before marking success as this is final stage
        # if [ ! -s "{params.output_final_dir}/{params.sample}/consensus_seq_stats.csv" ]; then
        #     echo "Expected consensus output missing or empty; failing rule" >&2
        #     exit 1
        # fi

        # mark complete
        touch "{output.castanet_done}"
        """

rule tar_castanet_output_dirs:
    input:
        config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/.castanet_done",
        config["scratch_dir"] + "/{sample_set}/pre-processed/qc/{sample}/{sample}_qc.csv" # ensure QC is done before tarring
    output:
        consenus_tar = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/consensus_data.tar.gz",
        depth_tar = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/Depth_output.tar.gz",
        grouped_reads_tar = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/grouped_reads.tar.gz",
        hashes_tar = config["scratch_dir"] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/hashes.tar.gz"
    resources:
        project = config["project"],
        queue = "normal",
        mem = "1G",
        walltime = "01:00:00",
        cpus = 1
    run:
        for dir_path, tar_path in zip([output.replace(".tar.gz", "") for output in output], output):
            shell(f"tar -czf {tar_path} -C {os.path.dirname(dir_path)} {os.path.basename(dir_path)}")
            shell(f"rm -rf {dir_path}")

rule extract_mosquito_cox1_reads_from_castanet:
    input:
        bam = ancient(config['scratch_dir'] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/{sample}.bam")
    output:
        r1_extracted = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_extracted/{fastq_type}/{sample}/{cox1_state}/{sample}_mosquito_cox1_{cox1_state}_R1.fastq"),
        r2_extracted = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_extracted/{fastq_type}/{sample}/{cox1_state}/{sample}_mosquito_cox1_{cox1_state}_R2.fastq"),
        r1_extracted_gz = config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_extracted/{fastq_type}/{sample}/{cox1_state}/{sample}_mosquito_cox1_{cox1_state}_R1.fastq.gz",
        r2_extracted_gz = config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_extracted/{fastq_type}/{sample}/{cox1_state}/{sample}_mosquito_cox1_{cox1_state}_R2.fastq.gz"
    resources:
        project = config["project"],
        queue = "normal",
        walltime = "02:00:00",
        mem = "1G",
        cpus = 1
    params:
        # mosquito_contigs = lambda wc: return_castanet_faidx_contigs(wc, config)
        # mosquito_contigs = lambda wc, input: get_mosquito_contigs_from_castanetv92_bam(wc, config, input.bam)
        mosquito_contig_string = "mosquitocox1"
    envmodules:
        "samtools",
        "bedtools"
    shell:
        """
        contigs=$(samtools view -H {input.bam} | awk '/^@SQ/{{gsub(/SN:/,""); print $2}}')
        contigs_to_extract=$(echo "$contigs" | grep -F "{params.mosquito_contig_string}" || true)

        samtools index {input.bam}
        if [ -z "$contigs_to_extract" ]; then
            : > {output.r1_extracted}
            : > {output.r2_extracted}
        else
            samtools view -b -F 4 {input.bam} $contigs_to_extract | samtools sort -n - | bedtools bamtofastq -i - -fq {output.r1_extracted} -fq2 {output.r2_extracted}
        fi
        gzip -c {output.r1_extracted} > {output.r1_extracted_gz}
        gzip -c {output.r2_extracted} > {output.r2_extracted_gz}
        """

rule report_unreachable_fraction_if_timmed:
    input:
        bam = lambda wc: expand(config['scratch_dir'] + "/{{sample_set}}/castanet/{version}/{fastq_type}/{sample}/{sample}.bam", sample = config['sample_sets'][wc.sample_set]['samples'])
    output:
        unreachable_fraction = config['scratch_dir'] + "/{sample_set}/pre-processed/{fastq_type}/{sample_set}_unreachable_fraction_if_trimmed.castanet{version}.tsv"
    resources:
        project = config["project"],
        queue = "normal",
        walltime = "02:00:00",
        mem = "2G",
        cpus = 1
    envmodules:
        "samtools",
    params:
        full_contigs = lambda wc: config["target_sets"]["castanet"],
        trimmed_contigs = lambda wc: [
            f"{full}:{start}-{end + 1}"
            for full, (start, end) in zip(
                config["target_sets"]["castanet"],
                (return_region_from_trimmed_id(region) 
                for region in config["target_sets"]["castanet_trimmed"])
            )
        ]
    shell:
        """
        echo -e "sample\ttotal_COX1_reads\textracted_reads_post_trim\tunreachable_fraction" > {output.unreachable_fraction}
        for sample_bam in {input.bam}; do
            sample=$(basename $(dirname $sample_bam))
            all_reads=$(samtools view -c $sample_bam {params.full_contigs})
            reachable_reads=$(samtools view -c -F 4 $sample_bam {params.trimmed_contigs})
            if [ $all_reads -eq 0 ]; then
                unreachable_fraction=0
            else
                unreachable_fraction=$(echo "scale=4; ($all_reads - $reachable_reads) / $all_reads" | bc)
            fi
            echo -e "$sample\t$all_reads\t$reachable_reads\t$unreachable_fraction" >> {output.unreachable_fraction}
        done
        """

rule extract_deduplicated_reads: # only from mosquito cox1 for now
    input:
        bam = ancient(config['scratch_dir'] + "/{sample_set}/castanet/{version}/{fastq_type}/{sample}/{sample}.bam")
    output:
        tmp_bam = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_dedup/{fastq_type}/{sample}/{cox1_state}/{sample}_dedup.tmp.bam"),
        tmp_bai = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_dedup/{fastq_type}/{sample}/{cox1_state}/{sample}_dedup.tmp.bam.bai"),
        r1_dedup = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_dedup/{fastq_type}/{sample}/{cox1_state}/{sample}_dedup_mosquito_cox1_{cox1_state}_R1.fastq"),
        r2_dedup = temp(config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_dedup/{fastq_type}/{sample}/{cox1_state}/{sample}_dedup_mosquito_cox1_{cox1_state}_R2.fastq"),
        r1_dedup_gz = config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_dedup/{fastq_type}/{sample}/{cox1_state}/{sample}_dedup_mosquito_cox1_{cox1_state}_R1.fastq.gz",
        r2_dedup_gz = config['scratch_dir'] + "/{sample_set}/pre-processed/castanet{version}_dedup/{fastq_type}/{sample}/{cox1_state}/{sample}_dedup_mosquito_cox1_{cox1_state}_R2.fastq.gz"
    resources:
        project = config["project"],
        queue = "normal",
        walltime = "02:00:00",
        mem = "1G",
        cpus = 1
    envmodules:
        "samtools",
        "bedtools",
        "picard"
    params:
        mosquito_contigs = lambda wc: return_castanet_faidx_contigs(wc, config),
        picard_jar = config["picard_jar"],
    shell:
        """
        set -euo pipefail
        # mark duplicate reads using picard (write to temp BAM so we can index before region filtering)
        java -jar {params.picard_jar} MarkDuplicates \
            I={input.bam} \
            O={output.tmp_bam} \
            M=/dev/stdout \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            CREATE_INDEX=false
        samtools index {output.tmp_bam}
        samtools view -b -F 1028 {output.tmp_bam} {params.mosquito_contigs} | \
        samtools sort -n - | \
        bedtools bamtofastq -i - -fq {output.r1_dedup} -fq2 {output.r2_dedup}
        gzip -c {output.r1_dedup} > {output.r1_dedup_gz}
        gzip -c {output.r2_dedup} > {output.r2_dedup_gz}
        """
        