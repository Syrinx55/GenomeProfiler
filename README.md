# GenomeProfiler

**W.I.P.**

Welcome to GenomeProfiler! GenomeProfiler is a pipeline for the profiling of prokaryotes and plasmids, including gene annotation across a wide range of features. GenomeProfiler provides features and nucleotide coordinates of metadata features for easy use in visual profiling and other analyses. Current features include antibiotic resistance genes, virulence genes, mobile genetic elements, finding related plasmids (for plasmid input), and more!

## Install

**NOTE:** Instructions are currently for `linux-64` and `osx-64` architectures.

### Install Git

Git is used to download GenomeProfiler.

Follow [official instructions](https://git-scm.com/downloads) to install Git.

Verify that Git is available by printing its version:
```bash
git --version
```

### Install Miniforge

Follow [official instructions](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) to install Miniforge.

Verify that Conda is available by printing its version:
```bash
conda --version
```

### Register Email with Entrez

Registering an email with Entrez is not necessary for use, and the `GENPROF_ENTREZ_EMAIL` environment variable (see [Install GenomeProfiler](#install-genomeprofiler)) may be set to any valid email address.
However, registration will allow for the processing of more accessions simultaneously in a given run.
To register, visit the [NCBI Entrez website](https://www.ncbi.nlm.nih.gov/search/) and click "Log in".
Create an account using your preferred email. If you do not wish to create an NCBI account, simply use a valid email.
A later step will instruct you to include this email address in your `.env` file.

### Get IslandViewer Token

*(For the* `islandviewer` *tool.)*

Visit the [IslandViewer website](https://www.pathogenomics.sfu.ca/islandviewer/) to create an account and obtain an API token.
Navigate to the "Login" tab and create an account or sign into an existing one.
Once you are signed in, navigate to the "Jobs" tab and click "HTTP API Token".
A later step will instruct you to include this token in your `.env` file.
*Note that this API key expires every 30 days and needs to be renewed.*

### Install GenomeProfiler

Navigate to the directory in which you want to install GenomeProfiler.

Clone GenomeProfiler from GitHub:
```bash
git clone https://github.com/Syrinx55/GenomeProfiler.git
```

Switch to release version 0.4.1:
```bash
git checkout release-0.4.1
```

Create and update the `genome-profiler` Conda environment:
```bash
conda env update --prune -f GenomeProfiler/environment.yml
```

Environment variables used by GenomeProfiler may be written in `.env`.
Copy the template `.env` into the current directory:
```bash
cp GenomeProfiler/template/.env .
  ```

Open `.env` in a text editor (GNU nano depicted here):
```bash
nano .env
```

The entries of the `.env` file are as follows: 
- `GENPROF_ENTREZ_EMAIL`: your valid email address.
- `GENPROF_ISLANDVIEWER_AUTH_TOKEN`: your API key from IslandViewer.

Delete the sample values (after `=`) and replace them with your own values. Save and exit.

Activate the `genome-profiler` Conda environment:
```bash
conda activate genome-profiler
```

Install resources (databases etc.) used by GenomeProfiler into the current directory:
```bash
GenomeProfiler/genome_profiler --setup
```

## Usage

### Graphical Mode

**W.I.P.**

Open the `GenomeProfiler` directory within a graphical file manager (e.g. Finder). Double-click `genome_profiler_gui`.

By default, GenomeProfiler will save its output within the parent directory of the `GenomeProfiler` directory.

### Command Line Interface

Activate `genome-profiler` Conda environment:
```bash
conda activate genome-profiler
```

Print the help message:
```bash
GenomeProfiler/genome_profiler --help
```

By default, GenomeProfiler will save its output within the current directory.

## Semantic Versioning

This repository adheres to Semantic Versioning 2.0.0: <https://semver.org/>.

The current version of this repository is stated in `./VERSION.txt`.

## Acknowledgements

- Benjamin Alford ([Ssolar5](https://github.com/ssolar5)): Packaging for distribution, command-line interface, documentation, cleanup.
