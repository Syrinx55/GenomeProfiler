# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import os
import sys
import argparse
from datetime import datetime

def generate_output_dir(base_dir, accession):
    # Returns a timestamped subdirectory path for a given accession.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(base_dir, f"{accession}_{timestamp}")
    os.makedirs(path, exist_ok=True)
    return path


def parse_cli_args():
    parser = argparse.ArgumentParser(description="GenomeProfiler CLI Pipeline")
    parser.add_argument("accessions", nargs="+", help="One or more NCBI accession IDs")
    parser.add_argument(
        "--include-tools",
        nargs="*",
        default=[],
        help="Only run these tools (e.g., abricate ectyper mobileog)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="Number of parallel workers (default: 2)",
    )
    parser.add_argument(
        "--output-dir",
        help="Base directory to store results (default from config)",
    )
    parser.add_argument(
        "--fasta",
        help="Optional local FASTA file to use instead of fetching from NCBI",
    )
    parser.add_argument(
        "--genbank",
        help="Optional local GenBank file to use instead of fetching from NCBI",
    )
    parser.add_argument(
        "--accession", help="(internal) Single accession for subprocess execution"
    )
    return parser.parse_args()
