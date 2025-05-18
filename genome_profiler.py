import argparse
from multiprocessing import cpu_count
# from configparser import ConfigParser


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="GenomeProfiler",
        description="A tool to automate genome feature profiling and visualization.",
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
        "-g",
        "--gui",
        action="store_true",
        help="Run GenomeProfiler in graphical mode. NOTE: This will ignore all other arguments.",
    )
    parser.add_argument(
        "--setup",
        action="store_true",
        help="Download resources required by tools, such as databases.",
    )
    parser.add_argument(
        "-t",
        "--include-tools",
        nargs="*",
        default=[],
        help="Run these tools on the accessions.",
    )
    parser.add_argument(
        "-a",
        "--all-tools",
        action="store_true",
        help="Run all available tools on the accessions."
    )
    parser.add_argument(
        "--exclude-tools",
        nargs="*",
        default=[],
        help="Do not run these tools on the accessions. NOTE: If a tool is included and then excluded, it will not be re-included."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=cpu_count(),
        help=f"Number of parallel workers (default: {cpu_count()})",
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
        "accessions",
        nargs="*",
        help="NCBI accession(s) to process (optional if using --fasta or --genbank)",
    )

    return parser


def parse_cli_args():
    pass


def main():
    # FIXME debug
    arg_parser = create_argument_parser()
    arg_parser.print_help()
    parsed_args = arg_parser.parse_args()
    print(parsed_args)


if __name__ == "__main__":
    main()
