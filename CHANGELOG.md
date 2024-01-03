# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [v1.1.1] - 2024-03-01

### Added
- Automatic setup of the PDB database during the first run
- Automatic parsing of both mmCIF and PDB files
- Removed `torch` dependency - significant improvement in startup times
- added models of v1.1 - available for the download via `get-models` subcommand
