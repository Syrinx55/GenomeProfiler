# GenomeProfiler

**W.I.P.**

Welcome to GenomeProfiler! GenomeProfiler is a pipeline for the profiling of prokaryotes and plasmids, including gene annotation across a wide range of features. GenomeProfiler provides features and nucleotide coordinates of metadata features for easy use in visual profiling and other analyses. Current features include antibiotic resistance genes, virulence genes, mobile genetic elements, finding related plasmids (for plasmid input), and more!

## Install

**NOTE:** Instructions are currently for `linux-64` and `osx-64`.

### Install Git

Git is used to download GenomeProfiler.

Follow official instructions to install Git: <https://git-scm.com/downloads>

### Install Miniforge

Follow official instructions to install Miniforge: <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>

Verify that Conda is available by printing its version:
```bash
conda --version
```

### Install GenomeProfiler

Clone GenomeProfiler from GitHub:
```bash
git clone https://github.com/Syrinx55/GenomeProfiler.git
```

Create and update `genome-profiler` Conda environment:
```bash
conda env update --prune -f GenomeProfiler/environment.yml
```

Environment variables used by GenomeProfiler may be written in `.env`.
Copy template `.env` into current directory:
```bash
cp GenomeProfiler/template/.env .
```

Open `.env` in a text editor (GNU nano depicted here):
```bash
nano .env
```

**TODO:** Provide instructions for getting values to use.

Delete the sample values (after `=`) and replace them with your own values. Save and exit.

Activate `genome-profiler` Conda environment:
```bash
conda activate genome-profiler
```

Download resources (databases etc.) used by GenomeProfiler:
```bash
GenomeProfiler/genome_profiler --setup
```

## Usage

Activate `genome-profiler` Conda environment:
```bash
conda activate genome-profiler
```

Print help message:
```bash
GenomeProfiler/genome_profiler --help
```

## Semantic Versioning

This repository adheres to Semantic Versioning 2.0.0: <https://semver.org/>.

The current version of this repository is stated in `./VERSION.txt`.
