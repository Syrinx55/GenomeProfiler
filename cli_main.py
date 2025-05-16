# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import os
import sys
from datetime import datetime
from dotenv import load_dotenv
from Bio import Entrez, SeqIO
from configparser import ConfigParser
from collection_pipeline import process_accession, validate_environment
from concurrent.futures import ThreadPoolExecutor
import subprocess
from pathlib import Path
from cli_interface import (
    parse_args_from_cli,
    generate_output_dir,
)

CONFIG_FILE = "config_genomeprofiler.ini"
SECTION = "brig_settings"


def extract_protein_to_nuc_coords(genbank_path):
    mapping = {}
    for record in SeqIO.parse(genbank_path, "genbank"):
        for feat in record.features:
            if feat.type == "CDS" and "protein_id" in feat.qualifiers:
                protein_id = feat.qualifiers["protein_id"][0]
                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                strand = feat.location.strand
                mapping[protein_id] = (start, end, strand)
    return mapping


def convert_phastest_hits_to_nucleotide(hit_file, genbank_path, out_file):
    protein_map = extract_protein_to_nuc_coords(genbank_path)
    import pandas as pd

    df = pd.read_csv(hit_file, sep="\t", header=None)
    df.columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]

    def get_coords(qid):
        base_id = qid.split("|")[0] if "|" in qid else qid
        return protein_map.get(base_id, (None, None, None))

    df[["nuc_start", "nuc_end", "strand"]] = df["qseqid"].apply(
        lambda x: pd.Series(get_coords(x))
    )
    df.to_csv(out_file, sep="\t", index=False)


def load_config():
    config = ConfigParser()
    if not config.read(CONFIG_FILE):
        raise FileNotFoundError(f"Missing config file: {CONFIG_FILE}")
    if not config.has_section(SECTION):
        raise KeyError(f"Missing section: {SECTION}")

    required_keys = [
        "output_base",
        "max_workers",
        "sleep_interval",
        "abricate_path",
        "integron_path",
        "isescan_path",
        "ectyper_path",
        "islandviewer_api_submit",
        "prodigal_path",
        "diamond_path",
        "mobileog_db_faa",
        "mobileog_db_csv",
        "nuccore_csv",
        "assembly_csv",
        "plsdb_sketch_path",
    ]
    optional_defaults = {
        "plsdb_timeout": "30",
        "plsdb_max_results": "100",
        "ectyper_cores": "2",
    }

    missing = [key for key in required_keys if not config.has_option(SECTION, key)]
    if missing:
        raise ValueError(f"Missing required config keys: {', '.join(missing)}")

    for key, value in optional_defaults.items():
        if not config.has_option(SECTION, key):
            config.set(SECTION, key, value)

    return config[SECTION]


# FIXME declaration shadows import (this declaration uses a name that already
# exists in an import statement)
def generate_output_dir(base, acc):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(base, f"{acc}_{timestamp}")
    os.makedirs(path, exist_ok=True)
    return path


def launch_with_logging(acc, config, args, output_dir, timestamped):
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, f"{acc}_pipeline.log")

    script_path = os.path.abspath(__file__)
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    command = [
        sys.executable,
        script_path,
        "--accession",
        acc,
        "--output-dir",
        output_dir,
    ]
    if args.include_tools:
        command += ["--include-tools"] + args.include_tools
    if args.fasta:
        command += ["--fasta", args.fasta]
    if args.genbank:
        command += ["--genbank", args.genbank]
    if timestamped:
        command += ["--timestamped-output"]

    with open(log_path, "w") as log_file:
        process = subprocess.Popen(
            command,
            stdout=log_file,
            stderr=log_file,
            env=env,
        )
        process.wait()
        return process.returncode == 0


def main():
    if "--accession" in sys.argv:
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--accession", required=False)
        parser.add_argument("--output-dir")
        parser.add_argument(
            "--include-tools",
            nargs="*",
            help="Tools to include: abricate, isescan, mobileog, ectyper, tncentral, integron, islandviewer, plsdb, phastest, parser",
        )
        parser.add_argument("--fasta")
        parser.add_argument("--genbank")
        parser.add_argument(
            "--timestamped-output",
            action="store_true",
            help="Enable timestamped output directories",
        )
        args = parser.parse_args()
        if not (args.accession or args.fasta or args.genbank):
            print(
                "[ERROR] You must provide either an accession or a local FASTA/GenBank file."
            )
            sys.exit(1)

        config = load_config()
        load_dotenv()
        validate_environment(config)
        Entrez.email = os.environ["GENPROF_ENTREZ_EMAIL"]
        acc_out = (
            generate_output_dir(config["output_base"], args.accession)
            if args.timestamped_output
            else os.path.join(config["output_base"], args.accession)
        )
        accession_to_use = args.accession if args.accession else "manual"
        process_accession(
            accession=accession_to_use,
            config=config,
            included_tools=args.include_tools,
            debug_textbox=None,
            fasta_override=args.fasta,
            genbank_override=args.genbank,
            output_override=acc_out,
        )
        if "parser" in (args.include_tools or []):
            from data_parser import run_parser

            run_parser(Path(acc_out))
        return
        
    config = load_config()
    validate_environment(config)
    Entrez.email = os.environ["GENPROF_ENTREZ_EMAIL"]

    args = parse_args_from_cli()
    output_base = args.output_dir or config["output_base"]

    print(
        f"\n[GenomeProfiler CLI] Starting analysis with tools: {args.include_tools}"
    )
    print(f"[INFO] Using {args.workers} worker(s)")
    print(f"[INFO] Output base directory: {output_base}\n")
    if args.include_tools and "phastest" in args.include_tools:
        print("[INFO] Phastest tool is enabled via --include-tools")

    accessions = args.accessions
    futures = []
    future_map = {}

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for acc in accessions:
            acc_out = (
                generate_output_dir(output_base, acc)
                if hasattr(args, "timestamped_output")
                else os.path.join(output_base, acc)
            )
            
    print("\n[GenomeProfiler CLI] All processing complete.\n")


if __name__ == "__main__":
    main()
