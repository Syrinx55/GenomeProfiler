# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import os
import argparse
from datetime import datetime


def generate_output_dir(base_dir, accession):
    # Returns a timestamped subdirectory path for a given accession.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(base_dir, f"{accession}_{timestamp}")
    os.makedirs(path, exist_ok=True)
    return path


def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description="GenomeProfiler CLI - A tool for automated genome feature profiling and visualization.",
        epilog=""" 
Available Tools (use with --include-tools):
  abricate       - Run abricate for resistance gene profiling
  isescan        - Run ISEScan to identify insertion sequences
  mobileog       - Run MobileOG-db analysis
  ectyper        - Run ECTyper for serotyping
  tncentral      - Run BLASTn search against TnCentral
  integron       - Run Integron Finder
  islandviewer   - Submit genome to IslandViewer for island prediction
  plsdb          - Use PLSDB for plasmid similarity screening
  phastest       - Run local PHASTEST search for phage regions
  parser         - Parse tool outputs for visualization (enabled by default in GUI)

Examples:
  python cli_main.py NZ_CP000000
  python cli_main.py NZ_CP000000 --include-tools abricate parser --timestamped-output
  python cli_main.py --help
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "accessions",
        nargs="*",
        help="NCBI accession(s) to process (optional if using --fasta or --genbank)",
    )
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
    args = parser.parse_args()
    if not args.accessions and not args.fasta and not args.genbank:
        parser.error(
            "You must provide at least one accession or use --fasta/--genbank."
        )
    return args


parse_args_from_cli = parse_command_line_args
