import hashlib

from meta_deepFRI.utils.fasta_file_io import load_fasta_file, write_fasta_file


def hash_sequence_id(seq_id: str):
    return hashlib.sha256(bytes(seq_id, encoding='utf-8')).hexdigest()


def encode_faa_ids(path):
    seq_records = load_fasta_file(path)
    hash_lookup_dict = {}

    for i in range(len(seq_records)):
        seq_id_hash = hash_sequence_id(seq_records[i].id)
        while seq_id_hash in hash_lookup_dict.keys():
            seq_id_hash = hash_sequence_id(seq_records[i].id + "1")
        hash_lookup_dict[seq_id_hash] = seq_records[i].id
        seq_records[i].id = seq_id_hash

    faa_hashed_ids_path = str(path) + ".hashed_ids"
    write_fasta_file(seq_records, faa_hashed_ids_path)
    return faa_hashed_ids_path, hash_lookup_dict
