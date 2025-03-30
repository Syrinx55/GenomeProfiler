# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import os
import re
import csv
import sys
import argparse
import time
import requests
import subprocess
import json
import gzip
import zipfile
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from configparser import ConfigParser
from ratelimit import limits, sleep_and_retry
from time import sleep
from plsdbapi import query
from requests_toolbelt.multipart.encoder import MultipartEncoder
from island_viewer import submit_islandviewer_job


# Point to the config file and settings
# This file contains your paths as well as your entrez email and API key
CONFIG_FILE = "config_plasmidviz.ini"
SECTION = "brig_settings"


def load_config():
    config = ConfigParser()
    if not config.read(CONFIG_FILE):
        raise FileNotFoundError(f"Missing config file: {CONFIG_FILE}")

    required_keys = [
        "entrez_email",
        "output_base",
        "max_workers",
        "sleep_interval",
        "abricate_path",
        "integron_path",
        "isescan_path",
        "ectyper_path",
        "islandviewer_api_submit",
        "islandviewer_auth_token",
        "prodigal_path",
        "diamond_path",
        "mobileog_db_faa",
        "mobileog_db_csv",
        # New lines: needed for the two-step merging
        "nuccore_csv",  # Must point to 'nuccore.csv' with columns [NUCCORE_ACC, ASSEMBLY_UID]
        "assembly_csv",  # Must point to 'assembly.csv' with columns [ASSEMBLY_UID, ASSEMBLY_ACC]
        "plsdb_sketch_path",  # The local .msh file for MASH
    ]

    optional_defaults = {
        "plsdb_timeout": "30",
        "plsdb_max_results": "100",
        "ectyper_cores": "2",
    }

    if not config.has_section(SECTION):
        raise KeyError(f"Missing section: {SECTION}")

    missing = [key for key in required_keys if not config.has_option(SECTION, key)]
    if missing:
        raise ValueError(f"Missing required config keys: {', '.join(missing)}")

    # Set defaults for optional parameters
    for key, value in optional_defaults.items():
        if not config.has_option(SECTION, key):
            config.set(SECTION, key, value)

    return config[SECTION]


def validate_environment(config):
    checks = [
        ("abricate", [config["abricate_path"], "--version"]),
        ("integron_finder", [config["integron_path"], "--version"]),
        ("isescan", [config["isescan_path"], "--version"]),
        ("ectyper", [config["ectyper_path"], "--version"]),
        ("diamond", [config["diamond_path"], "version"]),
        ("prodigal", [config["prodigal_path"], "-v"]),
    ]

    missing = []
    for name, cmd in checks:
        try:
            subprocess.run(
                cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing.append(name)

    if missing:
        raise EnvironmentError(f"Missing tools: {', '.join(missing)}")


@sleep_and_retry
@limits(calls=1, period=1)
def api_request(
    method, url, config, **kwargs
):  # Submit API request (ensure email is correct)
    headers = kwargs.pop("headers", {})
    headers.update(
        {"User-Agent": f"BRIG-Automator/1.0 (contact: {config['entrez_email']})"}
    )

    try:
        response = requests.request(method, url, headers=headers, timeout=30, **kwargs)
        response.raise_for_status()
        return response
    except requests.exceptions.RequestException as e:
        print(f"API Error: {str(e)}")
        return None


def fetch_ncbi_data(accession, config):
    try:
        # 1. FASTA
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="fasta", retmode="text"
        )
        fasta = handle.read()
        # 2. GenBank (full w/ parts)
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gbwithparts", retmode="text"
        )
        genbank = handle.read()
        return fasta, genbank
    except Exception as e:
        print(f"NCBI fetch error for {accession}: {str(e)}")
        return None, None


def ensure_blastdb_exists(db_fasta, db_out):
    # The BLAST db should be a nucleotide db.
    db_index = f"{db_out}.nhr"
    if not os.path.exists(db_index):
        print(f"[INFO] BLAST database not found at {db_index}. Creating database...")
        cmd = [
            "makeblastdb",
            "-in",
            db_fasta,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            db_out,
        ]
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] BLAST database created at {db_out}")
    else:
        print(f"[INFO] BLAST database already exists at {db_index}")


def run_blastn_tncentral(
    query_fasta, tncentral_db, output_file, evalue=1e-5, num_threads=4
):
    blastn_cmd = [
        "blastn",
        "-query",
        query_fasta,
        "-db",
        tncentral_db,
        "-out",
        output_file,
        "-outfmt",
        "6",  # Tabular format; adjust if you need a different format
        "-evalue",
        str(evalue),
        "-num_threads",
        str(num_threads),
    ]
    print(f"[INFO] Running BLASTN command: {' '.join(blastn_cmd)}")
    subprocess.run(blastn_cmd, check=True)
    print(f"[SUCCESS] BLASTN search complete. Results saved to {output_file}")


def run_abricate(fasta_path, output_dir, config):
    os.environ["BLAST_N_THREADS"] = "4"
    dbs = [
        "vfdb",
        "plasmidfinder",
        "card",
        "argannot",
        "resfinder",
    ]  # Can be adjusted to add/omit certain dbs from search

    for db in dbs:
        output_file = os.path.join(output_dir, f"{db}.tsv")
        cmd = [
            config["abricate_path"],
            "--db",
            db,
            "--nopath",
            "--minid",
            "80",  # Set minimum identity
            "--mincov",
            "80",  # Set minimum coverage
            fasta_path,
        ]
        subprocess.run(cmd, stdout=open(output_file, "w"), check=True)


def run_integron_finder(fasta_path, output_dir, config):
    try:
        cmd = [config["integron_path"], fasta_path, "--local", "--outdir", output_dir]
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Integron Finder failed: {str(e)}")


def run_isescan(fasta_path, output_dir, config):
    cmd = [
        config["isescan_path"],
        "--seqfile",
        fasta_path,
        "--output",
        output_dir,
        "--nthread",
        "4",
    ]
    subprocess.run(cmd, check=True)


def poll_islandviewer_job(token, config, sleep_interval=3, timeout=10000):
    status_url_template = config.get("islandviewer_api_status")
    if not status_url_template:
        raise ValueError("Missing 'islandviewer_api_status' in config.")

    status_url = status_url_template.format(token=token)
    print(f"[INFO] Polling job status at {status_url}...")
    start_time = time.time()
    while time.time() - start_time < timeout:
        response = requests.get(
            status_url, headers={"x-authtoken": config.get("islandviewer_auth_token")}
        )
        response.raise_for_status()
        status_data = response.json()
        current_status = status_data.get("status")
        print(f"[INFO] Current job status: {current_status}")
        if current_status.lower() in ["complete", "completed"]:
            download_url = status_data.get("download_url")
            if not download_url:
                # If the API did not supply a download URL, build it from a template.
                download_url_template = config.get("islandviewer_api_download")
                if download_url_template:
                    download_url = download_url_template.format(token=token)
                else:
                    raise ValueError(
                        "Job completed but no download URL was provided and no download template is set in config."
                    )
            return download_url
        elif current_status.lower() == "failed":
            raise Exception(f"Job failed: {status_data.get('message')}")
        time.sleep(sleep_interval)
    raise TimeoutError("Timeout waiting for IslandViewer job to complete.")


def download_islandviewer_results(download_url, output_dir, accession, config):
    output_file = os.path.join(output_dir, f"{accession}_islandviewer.tsv")
    headers = {"x-authtoken": config.get("islandviewer_auth_token")}
    print(f"[INFO] Downloading IslandViewer results from {download_url}...")
    response = requests.get(download_url, headers=headers)
    response.raise_for_status()
    with open(output_file, "w") as f:
        f.write(response.text)
    print(f"[SUCCESS] IslandViewer results downloaded to {output_file}")
    return output_file


def run_mash_screen_query(
    query_fasta,
    mash_sketch_path,
    output_tsv,
    max_pvalue=0.1,
    min_ident=0.99,
    threads=2,
    winner_takes_all=False,
):  # Parameters emulate default MASH settings on the PLSDB web tool.
    w_flag = ["-w"] if winner_takes_all else []
    cmd = [
        "mash",
        "screen",
        mash_sketch_path,
        query_fasta,
        "-v",
        str(max_pvalue),
        "-i",
        str(min_ident),
        "-p",
        str(threads),
    ] + w_flag

    print(f"Running MASH screen: {' '.join(cmd)}")
    with open(output_tsv, "w") as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
    print(f"MASH screen output saved to {output_tsv}")


def is_mash_table_empty(mash_output_tsv):
    try:
        with open(mash_output_tsv, "r") as f:
            # Collect all non-empty lines (after stripping whitespace)
            lines = [line.strip() for line in f if line.strip()]
        print(f"[DEBUG] Mash table line count: {len(lines)}")
        # If the file is empty, further processing with ECTyper is skipped.
        return len(lines) < 1
    except Exception as e:
        print(f"[ERROR] Unable to read mash output: {e}")
        return True


def download_plasmid_from_entrez(
    nuccore_acc, output_path_fasta, output_path_genbank, email="entrez_email"
):
    print(f"[INFO] Downloading plasmid {nuccore_acc} from Entrez...")
    Entrez.email = email
    try:
        # Download FASTA data
        handle = Entrez.efetch(
            db="nucleotide", id=nuccore_acc, rettype="fasta", retmode="text"
        )
        fasta_data = handle.read()
        with open(output_path_fasta, "w") as f:
            f.write(fasta_data)
        print(f"[SUCCESS] Plasmid FASTA {nuccore_acc} saved to {output_path_fasta}")

        # Download GenBank data
        handle = Entrez.efetch(
            db="nucleotide", id=nuccore_acc, rettype="gbwithparts", retmode="text"
        )
        genbank_data = handle.read()
        with open(output_path_genbank, "w") as f:
            f.write(genbank_data)
        print(f"[SUCCESS] Plasmid GenBank {nuccore_acc} saved to {output_path_genbank}")

    except Exception as e:
        print(f"[ERROR] Failed to download plasmid {nuccore_acc} from Entrez: {e}")


def download_plasmid_accessions(mash_output_tsv, output_dir, email):
    # Create a folder for plasmid sequences inside the given output directory
    plasmid_folder = os.path.join(output_dir, "plasmid_sequences")
    os.makedirs(plasmid_folder, exist_ok=True)
    print(f"[INFO] Plasmid download folder created at {plasmid_folder}")

    # Open the MASH output TSV file and read all lines.
    with open(mash_output_tsv, "r") as f:
        lines = f.readlines()

    # Assume the first line is a header and process the remaining lines.
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) < 5:
            continue  # Skip lines that don't have enough columns
        nuccore_acc = parts[4]  # Column 5: the nuccore accession

        # Debug print to check extracted accession
        print(f"[DEBUG] Extracted accession: {nuccore_acc}")

        # Validate accession format; adjust regex as needed to match expected formats.
        # This regex will allow for formats like NZ_CP149268.1, CP117215.1, NC_025138.1, etc.
        if not re.match(r"^(?:[A-Z]{1,2}_)?[A-Z]{2}\d+\.\d+$", nuccore_acc):
            print(f"[ERROR] Invalid NCBI accession format: {nuccore_acc}")
            continue  # Skip invalid accessions

        # Define output file paths for FASTA and GenBank data.
        output_fasta = os.path.join(plasmid_folder, f"{nuccore_acc}.fasta")
        output_genbank = os.path.join(plasmid_folder, f"{nuccore_acc}.gbk")

        # Call the combined function that downloads both FASTA and GenBank records.
        download_plasmid_from_entrez(nuccore_acc, output_fasta, output_genbank, email)


def run_ectyper(accession, output_dir, config):  # ECTyper is run on original accession
    handle = Entrez.efetch(
        db="nucleotide", id=accession, rettype="fasta", retmode="text"
    )
    fasta_data = handle.read()
    asm_path = os.path.join(output_dir, f"{accession}.fa")
    with open(asm_path, "w") as fa:
        fa.write(fasta_data)

    cmd = [
        config["ectyper_path"],
        "-i",
        asm_path,
        "-o",
        os.path.join(output_dir, accession),
        "--verify",
        "--cores",
        config.get("ectyper_cores", "2"),
    ]
    subprocess.run(cmd, check=True)


#  Two-Step Merging Function
def add_assembly_acc_two_step(
    mash_output_tsv, nuccore_csv, assembly_csv, out_tsv="mash_with_assembly.tsv"
):
    mash_data = []
    with open(mash_output_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("pvalue"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            mash_data.append(
                {
                    "pvalue": parts[0],
                    "proportion": parts[1],
                    "distance": parts[2],
                    "query_id": parts[3],
                    "ref_id": parts[4],  # Will be matched to NUCCORE_ACC
                }
            )
    mash_df = pd.DataFrame(mash_data)
    # rename 'ref_id' -> 'NUCCORE_ACC'
    mash_df.rename(columns={"ref_id": "NUCCORE_ACC"}, inplace=True)

    # load nuccore + assembly
    nuccore_df = pd.read_csv(nuccore_csv)
    assembly_df = pd.read_csv(assembly_csv)

    # Suppose nuccore.csv has columns: [NUCCORE_ACC, ASSEMBLY_UID]
    needed_nuccore = ["NUCCORE_ACC", "ASSEMBLY_UID"]
    sub_nuccore = nuccore_df[[c for c in needed_nuccore if c in nuccore_df.columns]]

    merged_step1 = pd.merge(mash_df, sub_nuccore, on="NUCCORE_ACC", how="left")

    # Suppose assembly.csv has columns: [ASSEMBLY_UID, ASSEMBLY_ACC]
    needed_assembly = ["ASSEMBLY_UID", "ASSEMBLY_ACC"]
    sub_assembly = assembly_df[[c for c in needed_assembly if c in assembly_df.columns]]

    merged_step2 = pd.merge(merged_step1, sub_assembly, on="ASSEMBLY_UID", how="left")

    merged_step2.to_csv(out_tsv, sep="\t", index=False)
    print(f"Final MASH table with assembly accession saved to {out_tsv}")


def is_cached(accession, dirs):
    # Check if processing can be skipped based on cached results.
    expected_outputs = [
        os.path.join(dirs["ectyper"], "ectyper_serotypes.tsv"),
        os.path.join(dirs["plsdb"], f"{accession}_mash_with_assembly.tsv"),
        os.path.join(dirs["assembly_fastas"]),
    ]
    return all(os.path.exists(f) for f in expected_outputs)


def download_assemblies(mash_asm_output, output_dir):
    # Download genome assemblies using the NCBI Datasets CLI tool.
    os.makedirs(output_dir, exist_ok=True)
    mash_df = pd.read_csv(mash_asm_output, sep="\t")

    for accession in mash_df["ASSEMBLY_ACC"].dropna().unique():
        output_path = os.path.join(output_dir, f"{accession}.zip")
        cmd = [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome",
            "--filename",
            output_path,
        ]
        print(f"[INFO] Downloading {accession} using datasets CLI...")
        try:
            subprocess.run(cmd, check=True)
            print(f"[SUCCESS] {accession} downloaded successfully!")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed to download {accession}: {e}")


def extract_assemblies(output_dir):
    # Extract downloaded genome FASTA files.
    fasta_files = []
    for file in os.listdir(output_dir):
        if file.endswith(".zip"):
            zip_path = os.path.join(output_dir, file)
            extract_path = os.path.join(output_dir, file.replace(".zip", ""))
            try:
                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    zip_ref.extractall(extract_path)
                print(f"[SUCCESS] Extracted to {extract_path}")
            except zipfile.BadZipFile:
                print(f"[ERROR] Invalid ZIP file: {zip_path}")
    for root, _, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".fna") or file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))
    return fasta_files


def run_ectyper_on_assemblies(fasta_files, output_dir, config, combined_results_file):
    # Run ECTyper on all extracted FASTA files and store results in a single file.
    os.makedirs(output_dir, exist_ok=True)

    with open(combined_results_file, "w") as combined_file:
        first_file = True  # To track the first file and include the header

        for fasta in fasta_files:
            base_name = os.path.basename(fasta).replace(".fna", "").replace(".fa", "")
            output_dir_path = os.path.join(output_dir, f"{base_name}_ectyper")
            os.makedirs(output_dir_path, exist_ok=True)  # Ensure the directory exists
            output_tsv = os.path.join(output_dir_path, "output.tsv")

            print(f"[INFO] Running eCTyper on {fasta}")

            cmd = [
                config["ectyper_path"],
                "-i",
                fasta,
                "-o",
                output_dir_path,
                "--verify",
                "--cores",
                str(config.get("ectyper_cores", "2")),
            ]

            try:
                subprocess.run(cmd, check=True)

                if os.path.exists(output_tsv):
                    with open(output_tsv, "r") as infile:
                        if first_file:
                            combined_file.write(infile.read())  # Write with header
                            first_file = False
                        else:
                            next(infile)  # Skip header for subsequent files
                            combined_file.write(infile.read())

                    print(f"[SUCCESS] ECTyper completed for {fasta}")

                else:
                    print(
                        f"[ERROR] ECTyper did not generate expected output: {output_tsv}"
                    )

            except subprocess.CalledProcessError as e:
                print(f"[ERROR] ECTyper failed on {fasta}: {e}")


def extract_serotype_from_subdir(subdir_path):
    # Extract the Serotype from the output.tsv in the subdirectory.
    output_tsv = os.path.join(subdir_path, "output.tsv")

    # Ensure the output file exists
    if not os.path.exists(output_tsv):
        print(f"[ERROR] {output_tsv} not found.")
        return None

    # Read the output.tsv and extract the 'Serotype' column
    df = pd.read_csv(output_tsv, sep="\t")

    # Check if 'Serotype' column exists
    if "Serotype" in df.columns:
        # Return the first non-null serotype value (assuming it's consistent for each assembly)
        serotype = df["Serotype"].dropna().iloc[0]
        return serotype
    else:
        print(f"[WARNING] 'Serotype' column not found in {output_tsv}.")
        return None


def create_serotype_mapping(mash_asm_output, ectyper_results_file, output_file):
    try:
        # 1. Load Data
        df_mash = pd.read_csv(mash_asm_output, sep="\t")
        ectyper_df = pd.read_csv(ectyper_results_file, sep="\t")

        # 2. Normalize Columns (case-insensitive)
        df_mash.columns = df_mash.columns.str.strip().str.lower()
        ectyper_df.columns = ectyper_df.columns.str.strip().str.lower()

        # 3. Handle Missing "file" Column
        # Use "name" if "file" is missing
        filename_col = "file" if "file" in ectyper_df.columns else "name"
        print(f"[DEBUG] Using column for filenames: {filename_col}")

        # 4. Extract Assembly Accession from "Name" (e.g., "GCF_027925845.1_ASM2792584v1_genomic" → "GCF_027925845.1")
        ectyper_df["assembly_acc"] = ectyper_df[filename_col].apply(
            lambda x: "_".join(
                x.split("_")[:2]
            )  # Split on underscores and take first two parts
        )

        # 5. Remove Version Suffix (e.g., "GCF_027925845.1" → "GCF_027925845")
        ectyper_df["assembly_acc_base"] = (
            ectyper_df["assembly_acc"].str.split(".").str[0]
        )
        df_mash["assembly_acc_base"] = df_mash["assembly_acc"].str.split(".").str[0]

        # 6. Map Serotypes
        serotype_map = ectyper_df.set_index("assembly_acc_base")["serotype"].to_dict()
        df_mash["serotype"] = df_mash["assembly_acc_base"].map(serotype_map)

        # 7. Save Results
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        df_mash[["nuccore_acc", "assembly_acc", "serotype"]].to_csv(
            output_file, sep="\t", index=False
        )
        print(f"[SUCCESS] Saved serotype mapping to {output_file}")

    except Exception as e:
        print(f"[ERROR] Failed to create mapping: {str(e)}")
        raise


def combine_ectyper_results(ectyper_output_dir, combined_output_file):
    with open(combined_output_file, "w") as outfile:
        first_file = True
        for subdir in os.listdir(ectyper_output_dir):
            subdir_path = os.path.join(ectyper_output_dir, subdir)
            if os.path.isdir(subdir_path) and subdir.endswith("_genomic_ectyper.json"):
                tsv_file_path = os.path.join(
                    subdir_path, "output.tsv"
                )  # Check for output.tsv file in each subdir
                if os.path.exists(tsv_file_path):
                    print(f"[INFO] Processing {tsv_file_path}...")
                    with open(tsv_file_path, "r") as infile:
                        if first_file:
                            outfile.write(
                                infile.read()
                            )  # Write header from the first file
                            first_file = False
                        else:
                            next(infile)  # Skip the header in subsequent files
                            outfile.write(infile.read())
                else:
                    print(f"[WARNING] No output.tsv found in {subdir_path}, skipping.")


def extract_sorted_cds(genbank_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    pos_path = os.path.join(output_dir, "positive_strand.fna")
    neg_path = os.path.join(output_dir, "negative_strand.fna")

    with open(pos_path, "w") as pos_f, open(neg_path, "w") as neg_f:
        for record in SeqIO.parse(genbank_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cds_seq = feature.extract(record.seq)
                    if feature.location.strand == 1:
                        pos_f.write(
                            f">{record.id}_{feature.qualifiers.get('locus_tag', ['unknown'])[0]}\n"
                        )
                        pos_f.write(str(cds_seq) + "\n")
                    elif feature.location.strand == -1:
                        cds_seq = cds_seq.reverse_complement()
                        neg_f.write(
                            f">{record.id}_{feature.qualifiers.get('locus_tag', ['unknown'])[0]}\n"
                        )
                        neg_f.write(str(cds_seq) + "\n")


def build_mobileog_db(output_dir, config):
    # Builds the MobileOG Diamond database if it doesn't already exist.
    # Define the path for the database; we use a fixed name (e.g., mobileOG-db.dmnd)
    db_path = os.path.join(output_dir, "mobileOG-db")
    # BLAST databases for nucleotides produce .nhr, .nin, .nsq files.
    # For a protein database (which Diamond uses), expect .dmnd file.
    if not os.path.exists(db_path + ".dmnd"):
        print(
            f"[INFO] MobileOG Diamond database not found at {db_path}. Building database..."
        )
        cmd = [
            config["diamond_path"],
            "makedb",
            "--in",
            config["mobileog_db_faa"],
            "-d",
            db_path,
        ]
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] MobileOG Diamond database created at {db_path}")
    else:
        print(f"[INFO] MobileOG Diamond database already exists at {db_path}")
    return db_path


def run_mobileog_search(protein_file, db_path, output_file, config):
    # Runs a Diamond blastp search of the protein_file against the MobileOG database.
    cmd = [
        config["diamond_path"],
        "blastp",
        "--query",
        protein_file,
        "--db",
        db_path,
        "--outfmt",
        "6",
        "--evalue",
        "1e-20",  # adjust parameters as needed
        "--id",
        "90",
        "--query-cover",
        "90",
        "--max-target-seqs",
        "1",
        "--out",
        output_file,
        "--threads",
        str(config.get("max_workers", "2")),
    ]
    print(f"[INFO] Running Diamond search: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[SUCCESS] Diamond search complete. Results saved to {output_file}")


def run_mobileogdb_genbank(genbank_path, output_dir, config):
    # Step 1: Extract proteins from the GenBank file
    protein_file = os.path.join(output_dir, "genbank_proteins.faa")
    with open(protein_file, "w") as out_f:
        for record in SeqIO.parse(genbank_path, "genbank"):
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                if "translation" not in feat.qualifiers:
                    continue
                protein = feat.qualifiers["translation"][0]
                locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
                product = feat.qualifiers.get("product", ["unknown"])[0]
                header = f">{locus_tag}|{product}"
                out_f.write(f"{header}\n{protein}\n")
    print(f"[INFO] Extracted proteins written to {protein_file}")

    # Step 2: Build (or verify) the MobileOG Diamond database
    db_path = build_mobileog_db(output_dir, config)

    # Step 3: Run Diamond search against the built MobileOG database
    diamond_out = os.path.join(output_dir, "mobileOG_results.tsv")
    run_mobileog_search(protein_file, db_path, diamond_out, config)


def analyze_plasmid_features(genbank_path, mobileog_tsv, output_dir):
    plasmid_features = []
    for record in SeqIO.parse(genbank_path, "genbank"):
        for feat in record.features:
            if feat.type == "plasmid":
                plasmid_features.append(
                    {"location": str(feat.location), "qualifiers": str(feat.qualifiers)}
                )
    with open(os.path.join(output_dir, "plasmid_features.json"), "w") as f:
        json.dump(plasmid_features, f)

    mobileog_results = pd.read_csv(mobileog_tsv)
    plasmid_genes = mobileog_results[mobileog_results["plasmid"].notna()]
    plasmid_genes.to_csv(os.path.join(output_dir, "plasmid_associated_genes.tsv"))


def setup_directories(accession, config):
    base_dir = os.path.join(config["output_base"], accession)
    dirs = {
        "base": base_dir,
        "cds": os.path.join(base_dir, "CDS"),
        "ncbi": os.path.join(base_dir, "ncbi"),
        "tncentral": os.path.join(base_dir, "tncentral"),
        "abricate": os.path.join(base_dir, "abricate"),
        "islandviewer": os.path.join(base_dir, "islandviewer"),
        "integron": os.path.join(base_dir, "integron_finder"),
        "isescan": os.path.join(base_dir, "isescan"),
        "plsdb": os.path.join(base_dir, "plsdb"),
        "ectyper": os.path.join(base_dir, "ectyper"),
        "ectyper_assemblies": os.path.join(base_dir, "ectyper_assemblies"),
        "mobileogdb": os.path.join(base_dir, "mobileogdb"),
        "assembly_fastas": os.path.join(base_dir, "assembly_fastas"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


# Main pipeline function
def process_accession(accession, config):
    start_time = datetime.now()
    dirs = setup_directories(accession, config)

    # 1. Fetch data from NCBI
    fasta, genbank = fetch_ncbi_data(accession, config)
    if not fasta or not genbank:
        return False

    fasta_path = os.path.join(dirs["ncbi"], f"{accession}.fasta")
    genbank_path = os.path.join(dirs["ncbi"], f"{accession}.gbk")
    with open(fasta_path, "w") as f:
        f.write(fasta)
    with open(genbank_path, "w") as f:
        f.write(genbank)

    try:
        print(f"[{accession}] Running BLASTn against tncentral")
        tncentral_fasta = config.get("tncentral_fasta")
        tncentral_db = config.get("tncentral_db")

        # Ensure the BLAST database exists:
        ensure_blastdb_exists(tncentral_fasta, tncentral_db)

        blast_output = os.path.join(
            dirs["tncentral"], f"{accession}_tncentral_blast.tsv"
        )
        run_blastn_tncentral(
            fasta_path,
            tncentral_db,
            blast_output,
            evalue=1e-5,
            num_threads=int(config["max_workers"]),
        )

        print(f"[{accession}] Running Abricate...")
        run_abricate(fasta_path, dirs["abricate"], config)

        print(f"[{accession}] Running ISEScan...")
        try:
            run_isescan(fasta_path, dirs["isescan"], config)
        except Exception as e:
            print(f"[{accession}] ISEScan failed, continuing: {e}")

        print(f"[{accession}] Extracting and sorting CDS features...")
        extract_sorted_cds(genbank_path, dirs["cds"])

        print(f"[{accession}] Running mobileOG-db analysis...")
        run_mobileogdb_genbank(genbank_path, dirs["mobileogdb"], config)

        print(f"[{accession}] Extracting and sorting CDS features again...")

        print(f"[{accession}] Running Integron Finder...")
        run_integron_finder(fasta_path, dirs["integron"], config)

        print(f"[{accession}] Running local Mash screen for PLSDB hits...")
        mash_sketch = config["plsdb_sketch_path"]
        mash_output = os.path.join(dirs["plsdb"], f"{accession}_mash.tsv")
        run_mash_screen_query(
            query_fasta=fasta_path,
            mash_sketch_path=mash_sketch,
            output_tsv=mash_output,
            max_pvalue=0.1,
            min_ident=0.99,
            threads=int(config["max_workers"]),
            winner_takes_all=False,
        )

        # Check if the Mash table is empty
        if is_mash_table_empty(mash_output):
            print(
                f"[{accession}] Mash table is empty. Skipping PLSDB-dependent steps (assemblies download, plasmid accessions, and ectyper on extracted assemblies)."
            )
        else:
            # Continue with PLSDB-dependent steps:
            mash_asm_output = os.path.join(
                dirs["plsdb"], f"{accession}_mash_with_assembly.tsv"
            )
            add_assembly_acc_two_step(
                mash_output_tsv=mash_output,
                nuccore_csv=config["nuccore_csv"],
                assembly_csv=config["assembly_csv"],
                out_tsv=mash_asm_output,
            )
            print(f"[{accession}] Using MASH output file: {mash_asm_output}")

            print(f"[{accession}] Downloading related plasmids: {mash_asm_output}")
            email = config["entrez_email"]
            download_plasmid_accessions(mash_output, dirs["plsdb"], email)

            print(f"[{accession}] Downloading assemblies using NCBI Datasets CLI...")
            download_assemblies(mash_asm_output, dirs["assembly_fastas"])

            print(f"[{accession}] Extracting FASTA files...")
            fasta_files = extract_assemblies(dirs["assembly_fastas"])

            print(f"[{accession}] Running eCTyper on extracted assemblies...")
            combined_ectyper_results_file = os.path.join(
                dirs["ectyper"], "ectyper_results_combined.tsv"
            )
            run_ectyper_on_assemblies(
                fasta_files, dirs["ectyper"], config, combined_ectyper_results_file
            )

            if not os.path.exists(combined_ectyper_results_file):
                print(
                    f"[ERROR] Combined eCTyper results file was not created: {combined_ectyper_results_file}"
                )

            print(f"[{accession}] Mapping serotypes to MASH output...")
            ectyper_results_file = os.path.join(
                dirs["ectyper"], "ectyper_results_combined.tsv"
            )
            serotype_mapping_file = os.path.join(
                dirs["ectyper"], "ectyper_serotypes.tsv"
            )
            create_serotype_mapping(
                mash_asm_output, ectyper_results_file, serotype_mapping_file
            )

            print(f"[{accession}] Running ectyper on the original accession...")
            run_ectyper(accession, dirs["ectyper"], config)

            print(f"[{accession}] Submitting genome to IslandViewer HTTP API...")
        response_data = submit_islandviewer_job(genbank_path, config)
        print(f"[{accession}] IslandViewer submission response: {response_data}")
        token = response_data.get("token")
        if not token:
            raise ValueError("No token returned from IslandViewer submission.")

        # Poll job for the download URL
        download_url = poll_islandviewer_job(
            token, config, sleep_interval=int(config.get("sleep_interval", "10"))
        )
        print(f"[{accession}] Download URL obtained: {download_url}")

        # Download the results and assign to islandviewer_tsv
        islandviewer_tsv = download_islandviewer_results(
            download_url, dirs["islandviewer"], accession, config
        )
        print(f"[{accession}] IslandViewer TSV available at {islandviewer_tsv}")

    except Exception as e:
        print(f"[{accession}] Pipeline failed: {str(e)}")
        return False

    duration = datetime.now() - start_time
    print(f"[{accession}] Completed in {duration}")
    return True


def main():
    try:
        config = load_config()
        validate_environment(config)
        Entrez.email = config["entrez_email"]

        parser = argparse.ArgumentParser(description="BRIG Automation Pipeline")
        parser.add_argument("accessions", nargs="+", help="NCBI accession numbers")
        parser.add_argument(
            "--workers",
            type=int,
            default=int(config["max_workers"]),
            help="Number of parallel workers",
        )
        args = parser.parse_args()

        print(
            f"Processing {len(args.accessions)} accessions with {args.workers} workers"
        )

        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(process_accession, acc, config): acc
                for acc in args.accessions
            }
            success_count = 0
            for future in as_completed(futures):
                acc = futures[future]
                try:
                    if future.result():
                        success_count += 1
                        print(f"[{acc}] Successfully processed")
                    else:
                        print(f"[{acc}] Failed processing")
                except Exception as e:
                    print(f"[{acc}] Critical error: {str(e)}")

        print(
            f"\nProcessing complete: {success_count}/{len(args.accessions)} succeeded"
        )

    except Exception as e:
        print(f"Fatal error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
