import gzip
import logging
from pathlib import Path

from pysam import tabix_compress

import mDeepFRI
from mDeepFRI.database import Database
from mDeepFRI.mmseqs import createdb, createindex
from mDeepFRI.utils import download_file

logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger(__name__)


def create_pdb_mmseqs():

    PDB100 = "https://wwwuser.gwdg.de/~compbiol/colabfold/pdb100_230517.fasta.gz"
    # check if pdb exists in a build dir
    build_dir = Path(mDeepFRI.__path__[0]).parent
    pdb100_path = str(build_dir / "pdb100_230517.fasta.gz")
    # remove additional suffix
    pdb100_path = Path(pdb100_path.split(".")[0])
    uncompressed_path = pdb100_path.with_suffix(".fasta")
    compressed_path = pdb100_path.with_suffix(".fasta.gz")

    if not (compressed_path).exists():
        logger.info("Downloading PDB100 database.")
        download_file(PDB100, compressed_path)

        # re-compress with bgzip (tabix_compress)
        with gzip.open(compressed_path,
                       "rb") as f_in, open(uncompressed_path, "wb") as f_out:
            f_out.write(f_in.read())

        tabix_compress(uncompressed_path, compressed_path, force=True)

        # remove uncompressed
        uncompressed_path.unlink()

    # create an MMSeqs database from PDB100
    # in a build directory
    pdb100_mmseqs = build_dir / "pdb100_230517.mmseqsDB"
    # check if database exists
    if not pdb100_mmseqs.exists():
        logging.info("Creating MMSeqs2 database from PDB100.")
        createdb(compressed_path, pdb100_mmseqs)
        createindex(pdb100_mmseqs)

    pdb_db = Database(foldcomp_db=pdb100_path.stem,
                      sequence_db=compressed_path,
                      mmseqs_db=pdb100_mmseqs)

    return pdb_db
