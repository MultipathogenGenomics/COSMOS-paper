
def count_reads_in_fastq(fastq_file: str) -> int:
    import gzip
    with gzip.open(fastq_file, 'rt') as f:
        return sum(1 for line in f) // 4
    
def filter_fasta_by_ids(fasta_file: str, ids_to_keep: set, output_file: str) -> None:
    from Bio import SeqIO
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in ids_to_keep:
                SeqIO.write(record, out_f, 'fasta')

def get_dict_key_from_value(d: dict, value):
    for key, val in d.items():
        if val == value:
            return key
    raise ValueError(f"Value {value} not found in dictionary")