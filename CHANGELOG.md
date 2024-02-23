# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
[Unreleased]: https://github.com/bioinf-MCB/Metagenomic-DeepFRI/compare/v1.1.2...HEAD

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
