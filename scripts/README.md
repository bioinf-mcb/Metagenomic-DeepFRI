
# Build database:


## structure_files.parse_structure_file.process_structure_file:
1. reads the structure file extracting sequence and atom positions
2. skip short sequences and truncate ones that size is over max_protein_length inside target_db_config.json
3. save extracted data:
- sequence - protein_id.faa file containing single sequence
```
f.write(f">{protein_id}\n{sequence}\n")
```
SEQ_ATOMS_DATASET_PATH / SEQUENCES / (protein_id + ".faa")
- atom positions - binary file containing positions of all atoms and correlated amino acid chain index
SEQ_ATOMS_DATASET_PATH / ATOMS / (protein_id + ".bin")
For more information on how those binary are saved check out source code at CPP_lib/atoms_file_io.h
