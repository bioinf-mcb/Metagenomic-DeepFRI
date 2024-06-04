# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.5...HEAD

## [Release]

## [1.1.6] - 2024-06-04
[1.1.6]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.5...v1.1.6

### Added
- MMSeqs2 API + tests
- rudimentary documentation website
- fixed bugs of v1.1.5


## [1.1.5] - 2024-04-01
[1.1.5]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.4...v1.1.5

### Added
- improved logging
- support for ESM databases (closes #80)
- corrected structures retrieval from PDB (closes #81, #83)
    - correct non-standard amino acids
    - remove base pairs

### Changed
- retrieval of structures from PDB & FoldComp databases
- optimized contact map alignment via Cython - improved search in large databases (ESM)


## [v1.1.4] - 2024-03-23
[v1.1.4]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.3...v1.1.4

### Added
- `--skip-pdb` flag in case of unexpected PDB bugs
- user parameter for allowed protein lengths `--min-length` and `--max-length`

### Fixes
- installation of PDB database
- alignment of faulty PDB structures
    - cases are reported via `logging`, require case-by-case investigation
    - faulty structures are then aligned to a predicted database, which are consistent
    - around 130/1300 alignments to PDB fail (10%)

## [v1.1.3] - 2024-03-13
[v1.1.3]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.2...v1.1.3

### Fixes
- bug where MMseqs2 returned empty results search [#73](https://github.com/bioinf-mcb/Metagenomic-DeepFRI/issues/73)

## [v1.1.2] - 2024-02-22
[v1.1.2]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.1...v1.1.2

### Added
- chunked loading of MMSeqs2 preselected targets - improves runtime

### Fixed
- C++17 lowered to C++14 for compatibility
- pyopal API v0.5 integrated
- removed `scipy` distance calculation with function
- replaced `biopython` with `biotite` for structure manipulation - better runtime and easier deploy
- tests fixed & added to CI/CD

## [v1.1.1] - 2024-03-01
[v1.1.1]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.0.0...v1.1.1

### Added

- Automatic setup of the PDB database during the first run
- Automatic parsing of both mmCIF and PDB files
- Removed `torch` dependency - significant improvement in startup times
- added models of v1.1 - available for the download via `get-models` subcommand
- CI/CD
