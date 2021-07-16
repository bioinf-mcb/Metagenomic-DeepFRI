import time

from CONFIG import *
from utils import run_command
import pandas as pd

from Bio import pairwise2
from Bio import SeqIO

from CPP_lib.libAtomDistanceIO import load_contact_map
from matplotlib import pyplot as plt

def main_loop():
    # while True:
    faa_files = list(QUERY_FOLDER_PATH.glob("**/*.faa"))
    # if len(faa_files) == 0:
    #     time.sleep(1)
    #     continue

    query_file = faa_files[0]

    creation_time = str(time.time())
    work_path = (WORK_PATH / creation_time)

    while work_path.exists():
        creation_time = str(int(time.time()))
        work_path = (WORK_PATH / creation_time)
    work_path.mkdir()

    target_databases = sorted(list(TARGET_MMSEQS_DATABASE_PATH.iterdir()))[-1]
    # todo if there are more than 3 databases remove oldest? count how many loops are using specific database?

    run_command(f"mmseqs createdb {query_file} {work_path / 'queryDB'} --dbtype 1")
    run_command(f"mmseqs search {work_path / 'queryDB'} {target_databases / DEFAULT_MMSEQS_NAME} {work_path / 'resultDB'} {TMP_FOLDER_PATH}")
    run_command(f"mmseqs convertalis {work_path / 'queryDB'} {target_databases / DEFAULT_MMSEQS_NAME} {work_path / 'resultDB'} {work_path / 'resultDB.m8'}")

    column_names= ["query","target","identity","alignment_length","mismatches","gap_openings","query_start","query_end","target_start","target_end","e_value","bit_score"]
    output_data = pd.read_csv(work_path / 'resultDB.m8', sep="\t", names=column_names)
    # delete work_path as we no longer need mmseqs query and result database




    output_data = output_data[output_data.identity > 0.95]

    all_query_sequences = dict()
    with open(query_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            all_query_sequences[record.id] = record.seq

    for i in range(len(output_data)):
        query_id = output_data["query"].iloc[i]
        target_id = output_data["target"].iloc[i]
        query_sequence = all_query_sequences[query_id]

        with open(CONTACT_MAP_DATASET_PATH/"seq"/(target_id + ".faa")) as f:
            for record in SeqIO.parse(f, "fasta"):
                target_sequence = record.seq

        target_cmap = load_contact_map(str(CONTACT_MAP_DATASET_PATH/"cmap"/(target_id + ".bin")),ANGSTROM_CONTACT_THRESHOLD)
        query_cmap = load_contact_map(str(CONTACT_MAP_DATASET_PATH/"cmap"/(query_id + ".bin")), ANGSTROM_CONTACT_THRESHOLD)

        fig = plt.figure(figsize=(16, 8))
        fig.add_subplot(1, 2, 1)
        plt.imshow(query_cmap, vmin=0, vmax=1)
        fig.add_subplot(1, 2, 2)
        plt.imshow(target_cmap, vmin=0, vmax=1)
        plt.savefig(f'/home/soliareofastora/data/TMP/testtest/{str(i)}.png')
        plt.show()


        alignments = pairwise2.align.globalxx(query_sequence, target_sequence)
        for match in alignments:
            print(match)




if __name__ == '__main__':
    main_loop()