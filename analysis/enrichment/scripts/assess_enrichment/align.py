# align reads to mt genomes

import os
import subprocess
from analysis.common.io import parse_json

def read_fasta(fasta_file: str) -> dict[str, str]:
    records = {}
    current_name = None
    current_sequence = []
    with open(fasta_file) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_sequence)
                current_name = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)
    if current_name is not None:
        records[current_name] = "".join(current_sequence)
    return records

def write_fasta_record(handle, name: str, sequence: str, line_width: int = 80) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), line_width):
        handle.write(f"{sequence[start:start + line_width]}\n")

def select_mt_genome_for_species(species: str, reference_config: dict, genus_defaults: dict) -> str:
    return select_mt_genome_record_for_species(species, reference_config, genus_defaults)["mt_genome"]

def select_mt_genome_record_for_species(species: str, reference_config: dict, genus_defaults: dict) -> dict:
    if species in reference_config['references']:
        return {
            "sample": species,
            "genus": species.split("_")[0],
            "reference_name": species,
            "mt_genome": reference_config['references'][species]['mt_genome'],
            "ref_is_species": True,
        }
    else:
        genus = species.split("_")[0]
        if genus in genus_defaults:
            reference_fasta = genus_defaults[genus]
            reference_name = os.path.basename(reference_fasta).split("_mt_genome.fasta")[0]
            return {
                "sample": species,
                "genus": genus,
                "reference_name": reference_name,
                "mt_genome": reference_fasta,
                "ref_is_species": False,
            }
        else:
            raise ValueError(f"No mt genome found for species {species} and no default for genus {genus}")

def align_reads_with_minimap(r1: str, r2: str, genome: str, output_prefix: str) -> None:
    command = f"minimap2 -ax sr {genome} {r1} {r2} | samtools view -bS - | samtools sort -o {output_prefix}.bam"
    subprocess.run(command, shell=True, check=True)

def align_reads_to_mt_genomes(reads_config: str, reference_config: str, genus_defaults: dict, output_dir: str) -> None:
    reads_config = parse_json(reads_config, remove_nest=True)
    reference_config = parse_json(reference_config, remove_nest=True)
    write_alignment_reference_map_from_configs(reads_config, reference_config, genus_defaults, output_dir)
    for species in reads_config['samples']:
        reads = reads_config['samples'][species]
        mt_genome = select_mt_genome_for_species(species, reference_config, genus_defaults)
        r1 = reads['R1']
        r2 = reads['R2']
        output_prefix = f"{output_dir}/{species}"
        align_reads_with_minimap(r1, r2, mt_genome, output_prefix)

def write_alignment_reference_map_from_configs(reads_config: dict, reference_config: dict, genus_defaults: dict, output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    output_file = f"{output_dir}/alignment_reference_map.tsv"
    rows = []
    for species in reads_config["samples"]:
        reference_record = select_mt_genome_record_for_species(species, reference_config, genus_defaults)
        rows.append({
            "sample": species,
            "bam": f"{output_dir}/{species}.bam",
            "reference_name": reference_record["reference_name"],
            "reference_fasta": reference_record["mt_genome"],
            "ref_is_species": str(reference_record["ref_is_species"]).lower(),
            "genus": reference_record["genus"],
        })
    with open(output_file, "w") as out_file:
        out_file.write("sample\tbam\treference_name\treference_fasta\tref_is_species\tgenus\n")
        for row in rows:
            out_file.write("\t".join([row["sample"], row["bam"], row["reference_name"], row["reference_fasta"], row["ref_is_species"], row["genus"]]) + "\n")
    return output_file

def write_alignment_reference_map(reads_config: str, reference_config: str, genus_defaults: dict, output_dir: str) -> str:
    reads_config = parse_json(reads_config, remove_nest=True)
    reference_config = parse_json(reference_config, remove_nest=True)
    return write_alignment_reference_map_from_configs(reads_config, reference_config, genus_defaults, output_dir)

def run_msa_of_mt_genomes(reference_config: str, output_dir: str) -> None:
    reference_config = parse_json(reference_config, remove_nest=True)
    os.makedirs(output_dir, exist_ok=True)
    joint_fasta = f"{output_dir}/references_joint.fasta"
    with open(joint_fasta, 'w') as out_f:
        for ref in reference_config["references"]:
            mt_genome = reference_config["references"][ref]["mt_genome"]
            records = read_fasta(mt_genome)
            if len(records) != 1:
                raise ValueError(f"Expected one sequence in {mt_genome}, found {len(records)}")
            sequence = next(iter(records.values()))
            write_fasta_record(out_f, ref, sequence)
    command = f"mafft --auto {joint_fasta} > {output_dir}/references_aligned.fasta"
    subprocess.run(command, shell=True, check=True)
