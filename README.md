# GenomeProfiler [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/genome_profiler/README.html) [![Anaconda-Server Badge](https://anaconda.org/bioconda/genome_profiler/badges/license.svg)](https://anaconda.org/bioconda/genome_profiler) [![Anaconda-Server Badge](https://anaconda.org/bioconda/genome_profiler/badges/platforms.svg)](https://anaconda.org/bioconda/genome_profiler) [![Anaconda-Server Badge](https://anaconda.org/bioconda/genome_profiler/badges/downloads.svg)](https://anaconda.org/bioconda/genome_profiler)

Welcome to GenomeProfiler! GenomeProfiler is an automated pipeline for the profiling of prokaryotes and plasmids, including gene annotation across a wide range of features. GenomeProfiler provides features and nucleotide coordinates of metadata features for easy use in visual profiling and other analyses. Current features include antibiotic resistance genes, virulence genes, mobile genetic elements, finding related plasmids (for plasmid input), and more.

## Install

**NOTE:** Instructions are currently for `linux-64` and `osx-64` architectures.

### Install Miniforge

Follow [official instructions](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install) to install Miniforge.

Verify that Conda is available by printing its version:
```bash
conda --version
```

### Install GenomeProfiler

If you do not already have a Conda environment, create one:
```bash
conda create --name myenv
```

Install the `genome_profiler` package from the `bioconda` channel:
```bash
conda install --name myenv -c bioconda genome_profiler
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
*Note that this API token expires every 30 days and needs to be renewed.*

### Install Resources Used by GenomeProfiler

Navigate to the directory in which you want to install GenomeProfiler.

Environment variables used by GenomeProfiler may be written in `.env`.
Create a file named `.env` and open it in a text editor (GNU nano depicted here):
```bash
nano .env
```

Insert the following text into `.env`:

```bash
# genome_profiler
GENPROF_ENTREZ_EMAIL=
GENPROF_ISLANDVIEWER_AUTH_TOKEN=
```

The entries of `.env` are as follows: 
- `GENPROF_ENTREZ_EMAIL`: your valid email address.
- `GENPROF_ISLANDVIEWER_AUTH_TOKEN`: your API key from IslandViewer.

Insert your own value after each corresponding `=`. Save and exit.

Activate the Conda environment which `genome_profiler` is installed to:
```bash
conda activate myenv
```

Install resources (databases etc.) used by GenomeProfiler to `./data_GenomeProfiler`:
```bash
genome_profiler --setup
```

## Usage

### Graphical Mode

Open GenomeProfiler in graphical mode:
```bash
genome_profiler --gui
```

By default, GenomeProfiler will save its output to `./output_GenomeProfiler`.

### Command Line Interface

Print the GenomeProfiler help message:
```bash
genome_profiler --help
```

By default, GenomeProfiler will save its output to `./output_GenomeProfiler`.

## Semantic Versioning

This repository adheres to Semantic Versioning 2.0.0: <https://semver.org/>.

## Acknowledgements

- Benjamin Alford ([Ssolar5](https://github.com/ssolar5)): Packaging for distribution, command-line interface, documentation, cleanup.
