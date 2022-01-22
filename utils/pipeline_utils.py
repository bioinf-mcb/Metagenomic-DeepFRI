import json
import os

from utils.utils import chunks, merge_files_binary


def split_query_into_jobs(query_records, work_path, parallel_jobs):
    job_paths = []
    jobs_records = chunks(query_records, parallel_jobs)
    for i in range(len(jobs_records)):
        job_path = (work_path / str(i))
        job_path.mkdir()
        job_paths.append(job_path)
        with open(job_path / "job_sequences.faa", "w") as f:
            for record in jobs_records[i]:
                f.write(f">{record.id}\n{record.seq}\n")
    return job_paths


def merge_result_files(work_path, finished_path):
    finished_path.mkdir(parents=True)

    os.system(f"cp {work_path}/merged_query_sequences.faa {finished_path}")
    os.system(f"cp -r {work_path}/query_files {finished_path}")

    alignments = {}
    for alignment_file in list(work_path.glob("*/alignments.json")):
        alignments.update(json.load(open(alignment_file, "r")))
    json.dump(alignments, open(finished_path / "alignments.json", "w"), indent=4, sort_keys=True)

    merge_files_binary(list(work_path.glob("*/mmseqs2_search_results.m8")), finished_path / "mmseqs2_search_results.m8")

    for deepfri_file in list(work_path.glob("*/results*")):
        if not (finished_path / deepfri_file.name).exists():
            os.system(f"cp {deepfri_file} {finished_path}")
            continue
        with open(deepfri_file, "r") as source:
            with open(finished_path / deepfri_file.name, "a") as dst:
                dst.writelines(source.readlines()[2:])

    for jobs_metadata_file in list(work_path.glob("*/metadata*")):
        metadata_store_path = finished_path / "jobs_metadata" / jobs_metadata_file.parent.name
        metadata_store_path.mkdir(parents=True, exist_ok=True)
        os.system(f"cp {jobs_metadata_file} {metadata_store_path}")
