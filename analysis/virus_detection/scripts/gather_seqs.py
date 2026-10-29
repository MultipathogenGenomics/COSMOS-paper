import os
from Bio import SeqIO
import polars as pl
from analysis.common.io import parse_json

def find_sample_entry_in_multikey_dict(sample_dict: dict, sample: str) -> dict:
    for nest in sample_dict:
        samples = sample_dict[nest].get("samples", sample_dict[nest])
        for item in samples:
            if item == sample:
                return samples[item]
    raise KeyError(f"Sample {sample} not found in sample config")

def gather_seqs_from_castanet_into_fasta(virus_name: str, castanet_nc2_file: str, sample_config: str, output_fasta: str, virus_acronyms: dict) -> None:
    # find samples with the virus detected in the castanet nc2 file
    # read the castanet nc2 file
    nc2_df = pl.read_csv(castanet_nc2_file, separator="\t")
    # convert cols to floats
    for col in nc2_df.columns:
        if col != "sampleid":
            nc2_df = nc2_df.with_columns(pl.col(col).cast(pl.Float64))
    # get sampleid where virus is detected (i.e. nc2 > 0)
    detected_samples = set(nc2_df.filter(pl.col(virus_name) > 0.5)["sampleid"].to_list())
    # parse sample config to get file path
    sample_config = parse_json(sample_config, remove_nest=True)
    # gather the sequences from the castanet output files
    with open(output_fasta, 'w') as out_f:
        for sample in detected_samples:
            sample_entry = find_sample_entry_in_multikey_dict(sample_config, sample)
            fastq_file = sample_entry['R1']
            fastq_dir = os.path.dirname(fastq_file)
            head_dir = os.path.dirname(fastq_dir)
            castanet_consensus_folder = os.path.join(head_dir, "castanet", "v9.2", "trimmed", f"{sample}", "consensus_sequences")
            consensus_file = os.path.join(castanet_consensus_folder, f"{virus_name}_remapped_consensus_sequence.fasta")
            if os.path.exists(consensus_file):
                for record in SeqIO.parse(consensus_file, "fasta"):
                    record.id = f"{sample}_{record.id.lstrip('>')}"
                    record.description = ""
                    SeqIO.write(record, out_f, "fasta")
            else:
                print(f"Consensus file not found for sample {sample} and virus {virus_name}: {consensus_file}")
