import logging
from pathlib import Path

import click

from mDeepFRI import __version__
from mDeepFRI.database import build_database
from mDeepFRI.pipeline import metagenomic_deepfri
from mDeepFRI.utils.utils import download_model_weights

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


@click.group()
@click.option("--debug/--no-debug", default=False)
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx, debug):
    """mDeepFRI"""

    ctx.ensure_object(dict)
    ctx.obj["debug"] = debug

    if debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Debug mode is on.")
        logging.getLogger("requests").setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logger.info("Debug mode is off.")
        logging.getLogger("requests").setLevel(logging.INFO)


@cli.command
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(exists=False),
    help="Path to folder where the database will be created.",
)
@click.pass_context
def get_models(ctx, output):
    """Download model weights for mDeepFRI."""

    if ctx.obj["debug"] is True:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("requests").setLevel(logging.DEBUG)
        logging.getLogger("urllib3").setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logging.getLogger("requests").setLevel(logging.INFO)
        logging.getLogger("urllib3").setLevel(logging.INFO)

    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    download_model_weights(output_path)


# add subcommand for building database
@cli.command()
@click.option(
    "-i",
    "--input",
    required=True,
    multiple=True,
    type=click.Path(exists=True),
    help="Path to structure files or folders containing structure files.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(exists=False),
    help="Path to folder where the database will be created.",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    type=int,
    help="Number of threads to use. Default is 1.",
)
@click.option(
    "-m",
    "--max-protein-length",
    default=1_000,
    type=int,
    help=
    "If protein is longer than this value, it will be truncated. Default is 1000 aa.",
)
@click.pass_context
def build_db(ctx, input, output, threads, max_protein_length):
    """Build a database for meta_deepFRI."""

    if ctx.obj["debug"] is True:
        logger.setLevel(logging.DEBUG)

    input_seqs = [Path(seqs) for seqs in input]
    output_path = Path(output)
    build_database(input_paths=input_seqs,
                   output_path=output_path,
                   overwrite=True,
                   threads=threads,
                   max_protein_length=max_protein_length)


@cli.command()
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="Path to an input protein sequences (FASTA file, may be gzipped).",
)
@click.option(
    "-d",
    "--db-path",
    required=True,
    type=click.Path(exists=True),
    help="Path to a structures database folder.",
)
@click.option(
    "-w",
    "--weights",
    required=True,
    type=click.Path(exists=True),
    help="Path to a folder containing model weights.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(exists=False),
    help="Path to output file.",
)
@click.option(
    "-f",
    "--output-format",
    default="tsv",
    type=click.Choice(["tsv", "csv", "json"]),
    help="Output format. Default is TSV.",
)
@click.option(
    "-p",
    "--processing-modes",
    default=["bp", "cc", "ec", "mf"],
    type=click.Choice(["bp", "cc", "ec", "mf"]),
    multiple=True,
    help="Processing modes. Default is all "
    "(biological process, cellular component, enzyme comission, molecular function).",
)
@click.option(
    "-a",
    "--angstrom-contact-thresh",
    default=6,
    type=float,
    help="Angstrom contact threshold. Default is 6.",
)
@click.option(
    "--generate-contacts",
    default=2,
    type=int,
    help="Gap fill threshold during contact map alignment.",
)
@click.option(
    "--mmseqs-min-bit-score",
    default=None,
    type=float,
    help="Minimum bit score for MMseqs2 alignment.",
)
@click.option(
    "--mmseqs-max-evalue",
    default=0.001,
    type=float,
    help="Maximum e-value for MMseqs2 alignment.",
)
@click.option(
    "--mmseqs-min-identity",
    default=0.5,
    type=float,
    help="Minimum identity for MMseqs2 alignment.",
)
@click.option(
    "--alignment-matrix",
    default="blosum62",
    type=click.Choice(
        ["blosum62", "pam250", "blosum80", "blosum45", "pam120", "pam160"]),
    help="Alignment matrix for contact map alignment.",
)
@click.option(
    "--alignment-gap-open",
    default=10,
    type=int,
    help="Gap open penalty for contact map alignment.",
)
@click.option(
    "--alignment-gap-extend",
    default=1,
    type=int,
    help="Gap extend penalty for contact map alignment.",
)
@click.option(
    "--alignment-min-identity",
    default=0.5,
    type=float,
    help="Minimum identity for contact map alignment.",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    type=int,
    help="Number of threads to use. Default is 1.",
)
@click.pass_context
def predict_function(ctx, input, db_path, weights, output, output_format,
                     processing_modes, angstrom_contact_thresh,
                     generate_contacts, mmseqs_min_bit_score,
                     mmseqs_max_evalue, mmseqs_min_identity, alignment_matrix,
                     alignment_gap_open, alignment_gap_extend,
                     alignment_min_identity, threads):
    """Predict protein function from sequence."""
    if ctx.obj["debug"] is True:
        logger.setLevel(logging.DEBUG)

    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    metagenomic_deepfri(Path(input), Path(db_path), Path(weights), output_path,
                        output_format, processing_modes,
                        angstrom_contact_thresh, generate_contacts,
                        mmseqs_min_bit_score, mmseqs_max_evalue,
                        mmseqs_min_identity, alignment_matrix,
                        alignment_gap_open, alignment_gap_extend,
                        alignment_min_identity, threads)


if __name__ == "__main__":
    cli(prog_name="mDeepFRI")
