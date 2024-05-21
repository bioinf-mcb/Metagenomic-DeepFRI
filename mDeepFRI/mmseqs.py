import csv
import os
import tempfile
from dataclasses import dataclass
from functools import partial
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Annotated, Dict, Iterable, List, Literal

import numpy as np
import numpy.lib.recfunctions as rfn
from pysam import FastxFile, tabix_compress

import mDeepFRI
from mDeepFRI.utils import retrieve_fasta_entries_as_dict, run_command


@dataclass
class ValueRange:
    min: float
    max: float


def _createdb(sequences_file, db_path):
    """
    Converts FASTA file to a DB format needed for MMseqs2 search/
    This function should generate five files,
    e.g. queryDB, queryDB_h and its corresponding index file queryDB.index,
    queryDB_h.index and queryDB.lookup from the FASTA QUERY.fasta input sequences.

    sequence_file (str): path to FASTA file.
    db_path (str): path to output db file.

    Returns:
        None
    """

    run_command(f"mmseqs createdb {sequences_file} {db_path} --dbtype 1")


def _createindex(db_path: str, threads: int = 1):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(
            f"mmseqs createindex {db_path} {tmp_path} --threads {threads}")


def _search(query_db: str,
            target_db: str,
            result_db: str,
            mmseqs_max_eval: float = 10e-5,
            sensitivity: Annotated[float, ValueRange(min=1.0, max=7.5)] = 5.7,
            threads: int = 1):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_command(
            f"mmseqs search -e {mmseqs_max_eval} --threads {threads} "
            f"-s {sensitivity} {query_db} {target_db} {result_db} {tmp_path}")


def _convertalis(
    query_db: str,
    target_db: str,
    result_db: str,
    output_file: str,
    columns: Literal["query", "target", "fident", "alnlen", "mismatch",
                     "gapopen", "qstart", "qend", "tstart", "tend", "qcov",
                     "tcov", "evalue", "bits", "qseq", "tseq", "qheader",
                     "theader", "qaln", "taln", "qframe", "tframe", "mismatch",
                     "qcov", "tcov", "qset", "qsetid", "tset", "tsetid",
                     "taxid", "taxname", "taxlineage", "qorfstart", "qorfend",
                     "torfstart", "torfend", "ppos"] = [
                         "query", "target", "fident", "alnlen", "mismatch",
                         "gapopen", "qstart", "qend", "tstart", "tend", "qcov",
                         "tcov", "evalue", "bits"
                     ]):

    args = ",".join(columns)
    run_command(
        f"mmseqs convertalis {query_db} {target_db} {result_db} {output_file} --format-mode 4 "
        f"--format-output {args}")


class MMSeqsSearchResult(np.recarray):
    """
    Class for handling MMSeqs2 search results. The results are stored in a TSV file.
    Inherits from numpy.recarray.

    Args:
        data (np.recarray): MMSeqs2 search results.
        query_fasta (str): Path to query FASTA file.
        database (str): Path to MMSeqs2 database.

    Attributes:
        data (np.recarray): MMSeqs2 search results.
        query_fasta (str): Path to query FASTA file.
        database (str): Path to MMSeqs2 database.
        columns (np.array): Array with column names.

    Example:

            >>> from mDeepFRI.mmseqs import MMSeqsSearchResult
            >>> result = MMSeqsSearchResult.from_filepath("path/to/file.tsv")
            >>> # sort file by identity and select seq1 hits only
            >>> result[::-1].sort(order=["fident"])
            >>> seq1_hits = result[result["query"] == "seq1"]
            >>> # save results
            >>> result.save("path/to/file.tsv")
    """
    def __init__(self, result_arr, query_fasta=None, database=None):
        self.result_arr = result_arr
        self.query_fasta = Path(query_fasta).resolve() if query_fasta else None
        self.database = Path(database).resolve() if database else None

    def __new__(cls, result_arr, query_fasta=None, database=None):
        obj = np.asarray(result_arr).view(cls)
        obj.query_fasta = query_fasta
        obj.database = database
        return obj

    @property
    def columns(self) -> np.ndarray:
        return np.array(self.result_arr.dtype.names)

    def save(self, filepath, filetype: Literal["tsv", "npz"] = "tsv"):
        """
        Save search results to TSV or NumPy compressed file.

        Args:
            filepath (str): Path to output file.
            filetype (str): File type to save. Options: "tsv" or "npz".

        Returns:
            None

        Example:

            >>> from mDeepFRI.mmseqs import MMSeqsSearchResult
            >>> result = MMSeqsSearchResult.from_filepath("path/to/file.tsv")
            >>> result.save("path/to/file.tsv")
        """
        # append query file column to the result array
        if self.query_fasta:
            query_col = np.array([self.query_fasta] * len(self.result_arr),
                                 dtype="U")
            self.result_arr = rfn.append_fields(self, "query_fasta", query_col)

        # append database column to the result array
        if self.database:
            db_col = np.array([self.database] * len(self.result_arr),
                              dtype="U")
            self.result_arr = rfn.append_fields(self.result_arr, "database",
                                                db_col)

        if filetype == "tsv":
            with open(filepath, "w", newline="") as f:
                # write tsv
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(self.result_arr.dtype.names)
                for row in self.result_arr:
                    writer.writerow(row)

        elif filetype == "npz":
            np.savez_compressed(filepath, self.result_arr)

        else:
            raise ValueError("File type should be 'tsv' or 'npz'.")

    def apply_filters(self,
                      min_cov: float = 0.0,
                      min_ident: float = 0.0,
                      min_bits: float = 0) -> "MMSeqsSearchResult":
        """
        Filter search results optionally by coverage, identity and bit score.

        Args:
            min_cov (float): Minimum coverage of query and target sequences.
            min_ident (float): Minimum identity of query and target sequences.
            min_bits (float): Minimum bit score.

        Returns:
            MMSeqsSearchResult: Filtered MMSeqs2 search results.

        Example:
            >>> from mDeepFRI.mmseqs import MMSeqsSearchResult
            >>> result = MMSeqsSearchResult.from_filepath("path/to/file.tsv")
            >>> filtered = result.apply_filters(min_cov=30, min_ident=0.4, min_bits=50)
        """

        self.result_arr = self.result_arr[(self.result_arr["qcov"] >= min_cov)
                                          &
                                          (self.result_arr["tcov"] >= min_cov)]
        self.result_arr = self.result_arr[
            self.result_arr["fident"] >= min_ident]
        self.result_arr = self.result_arr[self.result_arr["bits"] >= min_bits]

        return MMSeqsSearchResult(self.result_arr, self.query_fasta,
                                  self.database)

    def select_best_matches(self,
                            k: int = 5,
                            threads: int = 1) -> "MMSeqsSearchResult":
        """
        Selects k best matches for each query sequence based on bit score and identity.

        Args:
            k (int): Number of best matches to select.
            threads (int): Number of threads to use.

        Returns:
            MMSeqsSearchResult: MMSeqs2 search results with best matches.
        """
        def select_top_k(query, db, k=30):
            return db[db["query"] == query][:k]

        # sort by bit score
        self.result_arr.sort(order=["query", "bits", "fident"],
                             kind="quicksort")
        select_top = partial(select_top_k, db=self.result_arr[::-1], k=k)
        # select k best hits
        top_k = []
        with ThreadPool(threads) as p:
            top_k = p.map(select_top, np.unique(self.result_arr["query"]))

        return MMSeqsSearchResult(np.concatenate(top_k), self.query_fasta,
                                  self.database)

    @classmethod
    def from_filepath(cls,
                      filepath,
                      query_fasta=None,
                      database=None) -> "MMSeqsSearchResult":
        """
        Load search results from TSV file.

        Args:
            filepath (str): Path to TSV file from convertalis.
            query_fasta (str): Path to query FASTA file (optional).
            database (str): Path to MMSeqs2 database (optional).

        Returns:
            MMSeqsSearchResult: MMSeqs2 search results.

        Example:

                >>> from mDeepFRI.mmseqs import MMSeqsSearchResult
                >>> result = MMSeqsSearchResult.from_filepath("path/to/file.tsv")
        """

        result_arr = np.recfromcsv(filepath,
                                   delimiter="\t",
                                   encoding="utf-8",
                                   names=True)
        return cls(result_arr, query_fasta, database)


class QueryFile:
    """
    Class for handling FASTA files with sequences to query against MMSeqs2 database.

    Args:
        filepath (str): Path to FASTA file.

    Attributes:
        filepath (str): Path to FASTA file.
        sequences (Dict[str, str]): Dictionary with sequence IDs as keys and sequences as values.
    """
    def __init__(self, filepath: str) -> None:
        self.filepath: str = filepath
        self.sequences: Dict[str, str] = {}
        self.filtered_out: List[str] = []

    def __repr__(self) -> str:
        return f"QueryFile(filepath={self.filepath})"

    def __str__(self) -> str:
        return f"QueryFile(filepath={self.filepath})"

    def __setitem__(self, key, value) -> None:
        self.sequences[key] = value

    def __getitem__(self, key) -> None:
        return self.sequences[key]

    def load_ids(self, ids: Iterable[str]) -> None:
        """
        Load sequences by ID from FASTA file. The file is indexed with `samtools faidx`
        to speed up the process. Sequences are stored in a dictionary with
        sequence IDs as keys and sequences as values.

        Note:
            This method allows to load only sequences with specified IDs,
            which can be useful when working with large FASTA files.
            `samtools faidx` works with uncompressed FASTA files, and files compressed with bgzip.
            If compression is wrong, the file will be automatically re-compressed.

        Args:
            ids (List[str]): List of sequence IDs to load.

        Returns:
            None

        Raises:
            ValueError: If sequence with specified ID is not found in FASTA file.

        Example:
                >>> from mDeepFRI.mmseqs import QueryFile
                >>> query_file = QueryFile("path/to/file.fasta")
                >>> query_file.load_ids(["seq1", "seq2"])
        """
        # check if exists
        filepath = Path(self.filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File {self.filepath} not found.")
        self.sequences = retrieve_fasta_entries_as_dict(filepath, ids)

    def load_sequences(self, ids: Iterable[str] = None) -> None:
        """
        Load sequences from FASTA file. Sequences are stored in a dictionary with sequence
        IDs as keys and sequences as values.



        Note:
            This method should be called only if maniuplating sequences directly is needed.

        Returns:
            None

        Example:

            >>> from mDeepFRI.mmseqs import QueryFile
            >>> query_file = QueryFile("path/to/file.fasta")
            >>> query_file.load_sequences()
        """

        if ids:
            self.load_ids(ids)
        else:
            with FastxFile(self.filepath) as f:
                for entry in f:
                    self.sequences[entry.name] = entry.sequence

    def remove_sequences(self, ids: List[str]):
        """
        Remove sequences by ID.

        Args:
            ids (List[str]): List of sequence IDs to remove.

        Returns:
            None

        Example:

            >>> from mDeepFRI.mmseqs import QueryFile
            >>> query_file = QueryFile("path/to/file.fasta")
            >>> query_file.load_sequences()
            >>> query_file.remove_sequences(["seq1", "seq2"])
        """
        for seq_id in ids:
            try:
                self.sequences.pop(seq_id, None)
            except KeyError:
                raise ValueError(
                    f"Sequence with ID {seq_id} not found in {self.filepath}")

    def filter_sequences(self, condition: callable = None):
        """
        Filter sequences by a custom condition.

        Args:
            condition (callable): A lambda function that takes a sequence as input and returns a boolean value.

        Returns:
            None

        Example:

            >>> from mDeepFRI.mmseqs import QueryFile
            >>> query_file = QueryFile("path/to/file.fasta")
            >>> query_file.load_sequences()
            >>> query_file.filter_sequences(lambda seq: 50 <= len(seq) <= 200)
            >>> query_file.filter_sequences(lambda seq: 'ATG' in seq)  # filter by presence of 'ATG' motif
            >>> query_file.filter_sequences(lambda seq: seq.count('N') < 5) # filter by number of 'N' characters
        """
        # check if sequences were loaded
        if not self.sequences:
            raise ValueError(
                "No sequences loaded. Use load_sequences() or load_ids() method to load sequences from FASTA file."
            )

        filtered_sequences = self.sequences.copy()
        if condition:
            filtered_sequences = {
                k: v
                for k, v in filtered_sequences.items() if condition(v)
            }
            self.filtered_out = list(
                set(self.sequences.keys()) - set(filtered_sequences.keys()))
        else:
            self.filtered_out = []

        self.sequences = filtered_sequences

    def search(self,
               database_path: str,
               eval: float = 10e-5,
               sensitivity: Annotated[float,
                                      ValueRange(min=1.0, max=7.5)] = 5.7,
               index_target: bool = False,
               threads: int = 1):
        """
        Queries sequences against MMSeqs2 database. The search results are stored in a tabular format.

        Args:
            database_path (str): Path to MMSeqs2 database or database FASTA.
            eval (float): Maximum e-value for MMSeqs2 search.
            sensitivity (float): Sensitivity value for MMSeqs2 search.
            index_target (bool): Create index for target database. Advised for repeated searches.
            threads (int): Number of threads to use.

        Returns:
            MMSeqsSearchResult: MMSeqs2 search results.

        Example:

                >>> from mDeepFRI.mmseqs import QueryFile
                >>> query_file = QueryFile("path/to/file.fasta")
                >>> query_file.load_sequences()
                >>> query_file.filter_sequences(min_length=50, max_length=200)
                >>> result = query_file.search("path/to/database")
        """
        # check sensitivity values
        if not 1.0 <= sensitivity <= 7.5:
            raise ValueError(
                "Sensitivity value should be between 1.0 and 7.5.")

        with tempfile.TemporaryDirectory() as tmp_path:
            if self.sequences:
                fasta_path = Path(tmp_path) / "filtered_query.fa"
                with open(fasta_path, "w") as f:
                    for seq_id, seq in self.sequences.items():
                        f.write(f">{seq_id}\n{seq}\n")
            else:
                fasta_path = self.filepath

            # create query db
            input_db_path = Path(tmp_path) / "query.mmseqsDB"
            _createdb(fasta_path, input_db_path)

            # create target db
            with open(database_path, "r") as f:
                first_line = f.readline()

            if first_line.startswith(">"):
                target_db_path = Path(database_path).with_suffix(".mmseqsDB")
                _createdb(database_path, target_db_path)
                if index_target:
                    _createindex(target_db_path, threads)

            else:
                target_db_path = database_path

            result_db = Path(tmp_path) / "search_resultDB"
            _search(input_db_path, target_db_path, result_db, eval,
                    sensitivity, threads)

            output_file = Path(tmp_path) / "search_results.tsv"
            _convertalis(input_db_path, target_db_path, result_db, output_file)

            result = MMSeqsSearchResult.from_filepath(
                output_file,
                query_fasta=fasta_path,
                database=target_db_path,
            )

        return result


def extract_fasta_foldcomp(foldcomp_db: str,
                           output_file: str,
                           threads: int = 1):

    ESM_DATABASES = ["highquality_clust30", "esmatlas", "esmatlas_v2023_02"]
    """
    Extracts FASTA from Foldcomp database and compresses it using block gzip.

    Args:
        foldcomp_db (str): Path to Foldcomp database.
        output_file (str): Path to output FASTA file.
        threads (int): Number of threads to use.

    Returns:
        str: Path to gzipped FASTA file.
    """

    foldcomp_bin = Path(mDeepFRI.__path__[0]).parent / "foldcomp_bin"
    database_name = Path(foldcomp_db).stem

    # run command
    run_command(
        f"{foldcomp_bin} extract --fasta -t {threads} {foldcomp_db} {output_file}"
    )

    if database_name in ESM_DATABASES:
        # use sed to correct the headers
        os.system(
            fr"sed -i 's/^>\(ESMFOLD V0 PREDICTION FOR \)\(.*\)/>\2/' {output_file}"
        )

    # gzip fasta file
    tabix_compress(output_file, str(output_file) + ".gz", force=True)
    # remove unzipped file
    os.remove(output_file)
    # remove possible previous index, might lead to errror
    try:
        os.remove(str(output_file) + ".gz.fai")
        os.remove(str(output_file) + ".gz.gzi")
    except FileNotFoundError:
        pass

    return Path(str(output_file) + ".gz")


def validate_mmseqs_database(database: str):
    """
    Check if MMSeqs2 database is intact.

    Args:
        database (str): Path to MMSeqs2 database.

    Returns:
        bool: True if database is intact.

    """

    # Verify all the files for MMSeqs2 database
    mmseqs2_ext = [
        ".index", ".dbtype", "_h", "_h.index", "_h.dbtype", ".idx",
        ".idx.index", ".idx.dbtype", ".lookup", ".source"
    ]

    target_db = Path(database)
    is_valid = True
    for ext in mmseqs2_ext:
        if not os.path.isfile(f"{target_db}{ext}"):
            is_valid = False
            break

    return is_valid
