import logging
from pathlib import Path

import click

from mDeepFRI import __version__
from mDeepFRI.pipeline import predict_protein_function
from mDeepFRI.utils import download_model_weights, generate_config_json

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


@click.group()
@click.option("--debug/--no-debug", default=False)
@click.version_option(version=__version__)
@click.pass_context
def main(ctx, debug):
    """mDeepFRI"""

    loggers = [
        logging.getLogger(name) for name in logging.root.manager.loggerDict
    ]
    for log in loggers:
        if debug:
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.INFO)


@main.command
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(exists=False),
    help="Path to folder where the database will be created.",
)
@click.option("-v",
              "--version",
              required=True,
              type=click.Choice(["1.0", "1.1"]),
              help="Version of the model.")
@click.pass_context
def get_models(ctx, output, version):
    """Download model weights for mDeepFRI."""

    logger.info("Downloading DeepFRI models.")
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    download_model_weights(output_path, version)
    generate_config_json(output_path, version)
    logger.info(f"DeepFRI models v{version} downloaded to {output_path}.")


@main.command()
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
    multiple=True,
    help="Path to a structures database compessed with FoldComp.",
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
    "-p",
    "--processing-modes",
    default=["bp", "cc", "ec", "mf"],
    type=click.Choice(["bp", "cc", "ec", "mf"]),
    multiple=True,
    help="Processing modes. Default is all"
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
    "--mmseqs-min-bitscore",
    default=None,
    type=float,
    help="Minimum bitscore for MMseqs2 alignment.",
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
    "--top-k",
    default=5,
    type=int,
    help="Number of top MMSeqs2 alignment for"
    "precise pairwise alignment check. Default is 30.",
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
    "--remove-intermediate",
    default=False,
    type=bool,
    is_flag=True,
    help="Remove intermediate files.",
)
@click.option(
    "--overwrite",
    default=False,
    type=bool,
    is_flag=True,
    help="Overwrite existing files.",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    type=int,
    help="Number of threads to use. Default is 1.",
)
@click.pass_context
def predict_function(ctx, input, db_path, weights, output, processing_modes,
                     angstrom_contact_thresh, generate_contacts,
                     mmseqs_min_bitscore, mmseqs_max_evalue,
                     mmseqs_min_identity, top_k, alignment_gap_open,
                     alignment_gap_extend, alignment_min_identity,
                     remove_intermediate, overwrite, threads):
    """Predict protein function from sequence."""
    logger.info("Starting Metagenomic-DeepFRI.")

    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)

    predict_protein_function(
        query_file=input,
        databases=db_path,
        weights=weights,
        output_path=output_path,
        deepfri_processing_modes=processing_modes,
        angstrom_contact_threshold=angstrom_contact_thresh,
        generate_contacts=generate_contacts,
        mmseqs_min_bitscore=mmseqs_min_bitscore,
        mmseqs_max_eval=mmseqs_max_evalue,
        mmseqs_min_identity=mmseqs_min_identity,
        top_k=top_k,
        alignment_gap_open=alignment_gap_open,
        alignment_gap_continuation=alignment_gap_extend,
        identity_threshold=alignment_min_identity,
        remove_intermediate=remove_intermediate,
        overwrite=overwrite,
        threads=threads)


if __name__ == "__main__":
    main()
