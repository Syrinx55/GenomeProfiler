# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import os
import re
import time
import requests
import subprocess
import json
import zipfile
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from datetime import datetime
from ratelimit import limits, sleep_and_retry
from island_viewer import submit_islandviewer_job


@sleep_and_retry
@limits(calls=1, period=1)
def api_request(
    method, url, config, **kwargs
):  # Submit API request (ensure email is correct)
    headers = kwargs.pop("headers", {})
    headers.update(
        {
            "User-Agent": f"GenomeProfiler/1.0 (contact: {os.environ['GENPROF_ENTREZ_EMAIL']})"
        }
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
        Entrez.email = os.environ["GENPROF_ENTREZ_EMAIL"]

        # First, try to resolve accession via esearch if it's not already a valid ID
        print(f"[INFO] Resolving accession via Entrez search for: {accession}")
        search_handle = Entrez.esearch(db="nucleotide", term=accession, retmode="xml")
        search_results = Entrez.read(search_handle)
        ids = search_results.get("IdList", [])
        if not ids:
            print(f"[ERROR] No matching accession found for: {accession}")
            return None, None
        resolved_id = ids[0]
        print(f"[SUCCESS] Using resolved NCBI ID: {resolved_id}")

        # Fetch FASTA
        handle = Entrez.efetch(
            db="nucleotide", id=resolved_id, rettype="fasta", retmode="text"
        )
        fasta = handle.read()
        # Fetch GenBank
        handle = Entrez.efetch(
            db="nucleotide", id=resolved_id, rettype="gbwithparts", retmode="text"
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
        "ecoli_vf",
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
            status_url,
            headers={"x-authtoken": os.environ["GENPROF_ISLANDVIEWER_AUTH_TOKEN"]},
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
    headers = {"x-authtoken": os.environ["GENPROF_ISLANDVIEWER_AUTH_TOKEN"]}
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
                    # Extract gene name
                    gene_name = feature.qualifiers.get("gene", ["unknown_gene"])[0]

                    # Extract start and stop positions
                    start = (
                        int(feature.location.start) + 1
                    )  # Convert to 1-based indexing
                    stop = int(feature.location.end)  # End is exclusive in Biopython

                    # Extract the CDS sequence
                    cds_seq = feature.extract(record.seq)

                    # Determine the strand and write to the appropriate file
                    if feature.location.strand == 1:
                        pos_f.write(
                            f">{record.id}_{gene_name}_start{start}_stop{stop}\n"
                        )
                        pos_f.write(str(cds_seq) + "\n")
                    elif feature.location.strand == -1:
                        cds_seq = cds_seq.reverse_complement()
                        neg_f.write(
                            f">{record.id}_{gene_name}_start{start}_stop{stop}\n"
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
    """
    Extract proteins from the GenBank file to include
    nucleotide coordinates. Additionally, store the strand (1 for + and -1 for -)
    in the CDS mapping for later integration with MobileOG results.
    """
    protein_file = os.path.join(output_dir, "genbank_proteins.faa")
    cds_mapping = {}  # Mapping: composite_id -> {start: int, stop: int, strand: int}

    with open(protein_file, "w") as out_f:
        for record in SeqIO.parse(genbank_path, "genbank"):
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                if "translation" not in feat.qualifiers:
                    continue
                protein = feat.qualifiers["translation"][0]
                locus_tag = feat.qualifiers.get("locus_tag", ["unknown"])[0]
                product = feat.qualifiers.get("product", ["unknown"])[0]
                start = int(feat.location.start) + 1  # 1-based indexing
                stop = int(feat.location.end)  # end position as given
                strand = feat.location.strand  # typically 1 or -1
                # Create a composite id without strand (since it can be stored separately)
                composite_id = f"{record.id}_{locus_tag}_start{start}_stop{stop}"
                header = f">{composite_id}|{product}"
                out_f.write(f"{header}\n{protein}\n")
                # Store the coordinates and strand in the mapping dictionary
                cds_mapping[composite_id] = {
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                }

    print(
        f"[INFO] Extracted proteins with CDS coordinates and strand info to {protein_file}"
    )

    # Build (or verify) the MobileOG Diamond database and run the search.
    db_path = build_mobileog_db(output_dir, config)
    diamond_out = os.path.join(output_dir, "mobileOG_results.tsv")
    run_mobileog_search(protein_file, db_path, diamond_out, config)

    # Return both the Diamond results path and the CDS mapping (with strand info).
    return diamond_out, cds_mapping


def integrate_mobileog_annotations(mobileog_results_tsv, cds_mapping, output_file):
    """
    Integrates the MobileOG Diamond output with CDS nucleotide coordinates and strand info.
    Extract the composite identifier and join with the CDS mapping.
    """
    # Read the MobileOG Diamond results.
    # Diamond outfmt 6 output columns:
    # qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
    df = pd.read_csv(mobileog_results_tsv, sep="\t", header=None)
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

    # The query header (qseqid) is in the form: composite_id|product
    # Extract the composite_id by splitting on the pipe.
    df["composite_id"] = df["qseqid"].apply(lambda x: x.split("|")[0])

    # Map nucleotide coordinates and strand from the cds_mapping dictionary.
    df["nuc_start"] = df["composite_id"].apply(
        lambda x: cds_mapping.get(x, {}).get("start")
    )
    df["nuc_stop"] = df["composite_id"].apply(
        lambda x: cds_mapping.get(x, {}).get("stop")
    )
    # Convert strand 1 and -1 to human-readable '+' and '-' if desired.
    df["strand"] = df["composite_id"].apply(
        lambda x: (
            "+"
            if cds_mapping.get(x, {}).get("strand") == 1
            else "-"
            if cds_mapping.get(x, {}).get("strand") == -1
            else None
        )
    )

    # Save the integrated results to a TSV file.
    df.to_csv(output_file, sep="\t", index=False)
    print(
        f"[INFO] Integrated MobileOG annotations (with CDS positions and strand) saved to {output_file}"
    )


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


def setup_directories(accession, config, output_override=None):
    base_dir = output_override or os.path.join(config["output_base"], accession)
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
        "mobileogdb": os.path.join(base_dir, "mobileogdb"),
        "assembly_fastas": os.path.join(base_dir, "assembly_fastas"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


def extract_accession_from_file(file_path):
    # Parse the accession from a local FASTA or GenBank file.
    try:
        with open(file_path) as handle:
            record = next(
                SeqIO.parse(
                    handle, "fasta" if file_path.endswith(".fasta") else "genbank"
                )
            )
            return record.id.split()[0]
    except Exception as e:
        print(f"[ERROR] Could not extract accession from {file_path}: {e}")
        return None


# Main pipeline function
def process_accession(
    accession,
    config,
    included_tools=None,
    debug_textbox=None,
    fasta_override=None,
    genbank_override=None,
    output_override=None,
    timestamp_output=False,
):
    def log(msg):
        if debug_textbox:
            debug_textbox.insert("end", msg + "\n")
            debug_textbox.see("end")
        else:
            print(msg)

    start_time = datetime.now()
    fasta = None
    genbank = None

    if timestamp_output and not output_override:
        resolved_name = accession or "unknown"
        timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        output_override = os.path.join(
            config["output_base"], f"{resolved_name}_{timestamp}"
        )

    # Use output_override in setup
    dirs = setup_directories(accession, config, output_override=output_override)

    fasta_path = os.path.join(dirs["ncbi"], f"{accession}.fasta")
    genbank_path = os.path.join(dirs["ncbi"], f"{accession}.gbk")

    # FASTA and GenBank override handling with match check
    if (
        fasta_override
        and os.path.exists(fasta_override)
        and genbank_override
        and os.path.exists(genbank_override)
    ):
        log(f"[INFO] Using uploaded FASTA file: {fasta_override}")
        log(f"[INFO] Using uploaded GenBank file: {genbank_override}")

        with open(fasta_override) as f_fa, open(genbank_override) as f_gb:
            fasta = f_fa.read()
            genbank = f_gb.read()

        # Extract ID from FASTA header
        fasta_id = None
        try:
            first_line = fasta.splitlines()[0]
            if first_line.startswith(">"):
                fasta_id = first_line[1:].split()[0].strip()
        except Exception as e:
            log(f"[WARNING] Failed to extract ID from FASTA: {e}")

        # Extract ACCESSION from GenBank
        genbank_id = None
        try:
            for line in genbank.splitlines():
                if line.startswith("ACCESSION"):
                    parts = line.split()
                    if len(parts) > 1:
                        genbank_id = parts[1].strip()
                        break
        except Exception as e:
            log(f"[WARNING] Failed to extract accession from GenBank: {e}")

        # If both match, skip NCBI query; else try fallback via filename
        if fasta_id and genbank_id and fasta_id == genbank_id:
            accession = fasta_id
            log(
                f"[SUCCESS] FASTA and GenBank IDs match: {accession}. Skipping NCBI query."
            )
        else:
            # Fallback: try matching base filenames
            fasta_base = os.path.splitext(os.path.basename(fasta_override))[0]
            genbank_base = os.path.splitext(os.path.basename(genbank_override))[0]
            if fasta_base == genbank_base:
                accession = fasta_base
                log(
                    f"[SUCCESS] Using base filename match: {accession}. Skipping NCBI query."
                )
            else:
                log("[INFO] IDs do not match and fallback match also failed.")
                fasta, accession = None, None

    if not fasta and fasta_override and os.path.exists(fasta_override):
        log(f"[INFO] Using uploaded FASTA file: {fasta_override}")
        with open(fasta_override) as f:
            fasta = f.read()
        if not accession:
            try:
                first_line = fasta.splitlines()[0]
                if first_line.startswith(">"):
                    query_id = first_line[1:].split()[0].strip()
                    log(
                        f"[INFO] Attempting to resolve accession from header: {query_id}"
                    )
                    Entrez.email = os.environ["GENPROF_ENTREZ_EMAIL"]
                    search_handle = Entrez.esearch(
                        db="nucleotide", term=query_id, retmode="xml"
                    )
                    search_results = Entrez.read(search_handle)
                    ids = search_results.get("IdList", [])
                    if ids:
                        accession = ids[0].strip()
                        log(f"[SUCCESS] Found matching accession: {accession}")
                        fasta, _ = fetch_ncbi_data(accession, config)
                    else:
                        log(
                            f"[WARNING] No matching accession found for header: {query_id}"
                        )
                else:
                    log(f"[WARNING] FASTA header does not start with '>': {first_line}")
            except Exception as e:
                log(f"[WARNING] Could not resolve accession from FASTA header: {e}")
    elif not fasta and accession:
        fasta, _ = fetch_ncbi_data(accession, config)
        if not fasta:
            log(f"[ERROR] Could not fetch FASTA for {accession}")
            return False
    elif not fasta:
        log("[ERROR] No FASTA provided and no accession available.")
        return False

    # GenBank override handling
    if genbank_override and os.path.exists(genbank_override) and not genbank:
        log(f"[INFO] Using uploaded GenBank file: {genbank_override}")
        with open(genbank_override) as f:
            genbank = f.read()
        try:
            for line in genbank.splitlines():
                if line.startswith("ACCESSION"):
                    parts = line.split()
                    if len(parts) > 1:
                        accession = parts[1].strip()
                        accession = accession.strip()
                        log(f"[INFO] Extracted accession from GenBank: {accession}")
                        fetched_fasta, _ = fetch_ncbi_data(accession, config)
                        if not fetched_fasta:
                            log(f"[ERROR] Could not fetch FASTA for {accession}")
                            return False
                        fasta = fetched_fasta
                    break
        except Exception as e:
            log(f"[WARNING] Could not extract accession from GenBank: {e}")

        # Rebuild output directory based on updated accession
        if timestamp_output:
            timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
            output_override = os.path.join(
                config["output_base"], f"{accession}_{timestamp}"
            )
        dirs = setup_directories(accession, config, output_override=output_override)

        fasta_path = os.path.join(dirs["ncbi"], f"{accession}.fasta")
        genbank_path = os.path.join(dirs["ncbi"], f"{accession}.gbk")

        # Fetch the FASTA using the updated accession
        fetched_fasta, _ = fetch_ncbi_data(accession, config)
        if not fetched_fasta:
            log(f"[ERROR] Could not fetch FASTA for {accession}")
            return False
        fasta = fetched_fasta
    # If only a FASTA was uploaded, try to fetch the corresponding GenBank using the accession
    elif (
        fasta_override and os.path.exists(fasta_override) and not genbank and accession
    ):
        log(
            f"[INFO] Attempting to fetch GenBank using accession from FASTA: {accession}"
        )
        _, fetched_genbank = fetch_ncbi_data(accession, config)
        if not fetched_genbank:
            log(f"[ERROR] Could not fetch GenBank for {accession}")
            return False
        genbank = fetched_genbank
    elif accession and not genbank:
        accession = accession.strip()
        _, genbank = fetch_ncbi_data(accession, config)
        if not genbank:
            log(f"[ERROR] Could not fetch GenBank for {accession}")
            return False
    elif not genbank:
        log("[ERROR] No GenBank provided and no accession available.")
        return False

    with open(genbank_path, "w") as f:
        f.write(genbank)

    # Write the final, validated FASTA to file
    with open(fasta_path, "w") as f:
        f.write(fasta)

    try:
        log(f"[{accession}] Running BLASTn against tncentral")
        tncentral_fasta = config.get("tncentral_fasta")
        tncentral_db = config.get("tncentral_db")

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

        if not included_tools or "abricate" in included_tools:
            log(f"[{accession}] Running Abricate...")
            run_abricate(fasta_path, dirs["abricate"], config)

        if not included_tools or "isescan" in included_tools:
            log(f"[{accession}] Running ISEScan...")
            try:
                run_isescan(fasta_path, dirs["isescan"], config)
            except Exception as e:
                log(f"[{accession}] ISEScan failed, continuing: {e}")

        log(f"[{accession}] Extracting and sorting CDS features...")
        extract_sorted_cds(genbank_path, dirs["cds"])

        if not included_tools or "mobileog" in included_tools:
            log(
                f"[{accession}] Running mobileOG-db analysis with CDS coordinate enrichment..."
            )
            mobileog_results_tsv, cds_mapping = run_mobileogdb_genbank(
                genbank_path, dirs["mobileogdb"], config
            )
            integrated_output = os.path.join(
                dirs["mobileogdb"], "mobileOG_integrated.tsv"
            )
            integrate_mobileog_annotations(
                mobileog_results_tsv, cds_mapping, integrated_output
            )

        if not included_tools or "integron finder" in included_tools:
            log(f"[{accession}] Running Integron Finder...")
            run_integron_finder(fasta_path, dirs["integron"], config)

        if not included_tools or "plsdb" in included_tools:
            log(f"[{accession}] Running local Mash screen for PLSDB hits...")
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

            if not is_mash_table_empty(mash_output):
                mash_asm_output = os.path.join(
                    dirs["plsdb"], f"{accession}_mash_with_assembly.tsv"
                )
                add_assembly_acc_two_step(
                    mash_output_tsv=mash_output,
                    nuccore_csv=config["nuccore_csv"],
                    assembly_csv=config["assembly_csv"],
                    out_tsv=mash_asm_output,
                )
                log(f"[{accession}] Using MASH output file: {mash_asm_output}")

                email = os.environ["GENPROF_ENTREZ_EMAIL"]
                download_plasmid_accessions(mash_output, dirs["plsdb"], email)

                log(f"[{accession}] Downloading assemblies using NCBI Datasets CLI...")
                download_assemblies(mash_asm_output, dirs["assembly_fastas"])

                fasta_files = extract_assemblies(dirs["assembly_fastas"])
                log(f"[{accession}] Running eCTyper on extracted assemblies...")
                combined_ectyper_results_file = os.path.join(
                    dirs["ectyper"], "ectyper_results_combined.tsv"
                )
                run_ectyper_on_assemblies(
                    fasta_files, dirs["ectyper"], config, combined_ectyper_results_file
                )

                if os.path.exists(combined_ectyper_results_file):
                    ectyper_results_file = os.path.join(
                        dirs["ectyper"], "ectyper_results_combined.tsv"
                    )
                    serotype_mapping_file = os.path.join(
                        dirs["ectyper"], "ectyper_serotypes.tsv"
                    )
                    create_serotype_mapping(
                        mash_asm_output, ectyper_results_file, serotype_mapping_file
                    )
                else:
                    log(
                        f"[ERROR] Combined eCTyper results file was not created: {combined_ectyper_results_file}"
                    )

        if not included_tools or "ectyper" in included_tools:
            log(f"[{accession}] Running ectyper on the original accession...")
            run_ectyper(accession, dirs["ectyper"], config)

        if not included_tools or "islandviewer" in included_tools:
            log(f"[{accession}] Submitting genome to IslandViewer HTTP API...")
            response_data = submit_islandviewer_job(genbank_path, config)
            log(f"[{accession}] IslandViewer submission response: {response_data}")
            token = response_data.get("token")
            if not token:
                raise ValueError("No token returned from IslandViewer submission.")

            download_url = poll_islandviewer_job(
                token, config, sleep_interval=int(config.get("sleep_interval", "10"))
            )
            log(f"[{accession}] Download URL obtained: {download_url}")

            islandviewer_tsv = download_islandviewer_results(
                download_url, dirs["islandviewer"], accession, config
            )
            log(f"[{accession}] IslandViewer TSV available at {islandviewer_tsv}")

    except Exception as e:
        log(f"[{accession}] Pipeline failed: {str(e)}")
        return False

    duration = datetime.now() - start_time
    log(f"[{accession}] Completed in {duration}")
    return True
