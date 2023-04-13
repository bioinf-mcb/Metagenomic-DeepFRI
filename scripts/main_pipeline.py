import argparse
from pathlib import Path

from meta_deepFRI.metagenomic_deepfri import metagenomic_deepfri


def parse_args():

    parser = argparse.ArgumentParser(
        description="This script contains all the logic to run DeepFri's GCN or CNN experiments.")

    # logic described here is implemented in parse_input_paths
    parser.add_argument("-i",
                        "--input",
                        nargs='+',
                        required=True,
                        default=None,
                        help="List of filepaths containing query protein FASTA files (.faa or .faa.gz)")
    parser.add_argument("-db", "--db_path", required=True, default=None, help="Database path.")
    parser.add_argument("-c", "--config", required=True, default=None, help="Path to a config file from models folder.")
    parser.add_argument("-o", "--output_path", required=True, default=None, help="Output directory path.")
    parser.add_argument("--output_format",
                        nargs="*",
                        choices=['tsv', 'json', 'csv'],
                        default=['tsv'],
                        help="Output format.")

    # DeepFRI parameters
    parser.add_argument("--processing_modes", nargs="*", choices=['mf', 'bp', 'cc', 'ec'],
                        default=['mf', 'bp', 'cc', 'ec'],
                        help="DeepFRI processing modes. mf - molecular function, " \
                            "bp - biological process, cc - cellular component, " \
                            "ec - enzyme comission.")

    # Contact map alignment parameters
    parser.add_argument("--angstrom_contact_thresh",
                        required=False,
                        default=6,
                        type=int,
                        help="Angstrom contact threshold.")
    parser.add_argument("--generate_contacts",
                        required=False,
                        default=2,
                        type=int,
                        help="Gap fill during contact map alignment.")

    # MMSeqs2 parameters
    parser.add_argument("--mmseqs_min_bit_score",
                        required=False,
                        default=None,
                        type=float,
                        help="Minimum bit score for MMSeqs2 search.")
    parser.add_argument("--mmseqs_max_evalue",
                        required=False,
                        default=None,
                        type=float,
                        help="Maximum e-value for MMSeqs2 search.")
    parser.add_argument("--mmseqs_min_identity",
                        required=False,
                        default=0.5,
                        type=float,
                        help="Minimum sequence identity for MMSeqs2 search.")

    # Pairwise alignment parameters
    parser.add_argument("--alignment_match",
                        required=False,
                        default=2,
                        type=float,
                        help="Match score for pairwise sequence alignment.")
    parser.add_argument("--alignment_missmatch",
                        required=False,
                        default=-1,
                        type=float,
                        help="Missmatch score for pairwise sequence alignment.")
    parser.add_argument("--alignment_gap_open",
                        required=False,
                        default=-1,
                        type=float,
                        help="Gap open score for pairwise sequence alignment.")
    parser.add_argument("--alignment_gap_continuation",
                        required=False,
                        default=-0.1,
                        type=float,
                        help="Gap continuation score for pairwise sequence alignment.")
    parser.add_argument("--alignment_min_identity",
                        required=False,
                        default=0.5,
                        type=float,
                        help="Minimum sequence identity for pairwise sequence alignment.")
    parser.add_argument("-t",
                        "--threads",
                        required=False,
                        default=1,
                        type=int,
                        help="Number of threads to use for parallel processing.")

    return parser.parse_args()


def main():
    args = parse_args()

    metagenomic_deepfri(Path(args.input), Path(args.db_path), Path(args.config), Path(args.output_path),
                        args.processing_modes, args.angstrom_contact_thresh, args.generate_contacts, args.output_format,
                        args.mmseqs_min_bit_score, args.mmseqs_max_evalue, args.mmseqs_min_identity,
                        args.alignment_match, args.alignment_missmatch, args.alignment_gap_open,
                        args.alignment_gap_continuation, args.alignment_min_identity, args.threads)


if __name__ == '__main__':
    main()
