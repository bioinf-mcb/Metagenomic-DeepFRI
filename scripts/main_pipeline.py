import argparse
from pathlib import Path

from meta_deepFRI.metagenomic_deepfri import metagenomic_deepfri


def parse_args():
    parser = argparse.ArgumentParser(
    description="This script for prediction of protein function from sequence using" \
                "structure database, contact map alignment and DeepFRI model.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # logic described here is implemented in parse_input_paths
    parser.add_argument("-i",
                        "--input",
                        required=True,
                        default=None,
                        help="Filepath to a query protein FASTA files (.faa or .faa.gz)")
    parser.add_argument("-db", "--db_path", required=True, default=None, help="Database path.")
    parser.add_argument("-o", "--output_path", required=True, default=None, help="Output directory path.")
    parser.add_argument("--output_format",
                        nargs="*",
                        choices=['tsv', 'json', 'csv'],
                        default=['tsv'],
                        help="Output format.")

    # DeepFRI parameters
    parser.add_argument("-w",
                        "--weights",
                        required=False,
                        default="trained_models",
                        help="Path to a folder containing downloaded deepfri weights and model_config.json.")
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
    parser.add_argument("--alignment_matrix",
                        required=False,
                        default="blosum62",
                        type=str,
                        help="Substitution matrix for pairwise sequence alignment.")
    parser.add_argument("--alignment_gap_open",
                        required=False,
                        default=11,
                        type=int,
                        help="Gap open score penalty for pairwise sequence alignment.")
    parser.add_argument("--alignment_gap_extend",
                        required=False,
                        default=1,
                        type=int,
                        help="Gap continuation score penalty for pairwise sequence alignment.")
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
    print("""
    If you use this software please cite:
    - Gligorijević et al. "Structure-based protein function prediction using graph convolutional networks" Nat. Comms. (2021). https://doi.org/10.1038/s41467-021-23303-9
    - Steinegger & Söding "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets" Nat. Biotechnol. https://doi.org/10.1038/nbt.3988
    - Maranga et al. "Comprehensive Functional Annotation of Metagenomes and Microbial Genomes Using a Deep Learning-Based Method" mSystems (2023) https://doi.org/10.1128/msystems.01178-22
    - Daily "Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments" BMC Bioinform. (2016) https://doi.org/10.1186/s12859-016-0930-z
    """)
    args = parse_args()
    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    metagenomic_deepfri(Path(args.input), Path(args.db_path), Path(args.weights), output_path, args.output_format,
                        args.processing_modes, args.angstrom_contact_thresh, args.generate_contacts,
                        args.mmseqs_min_bit_score, args.mmseqs_max_evalue, args.mmseqs_min_identity,
                        args.alignment_matrix, args.alignment_gap_open, args.alignment_gap_extend,
                        args.alignment_min_identity, args.threads)


if __name__ == '__main__':
    main()
