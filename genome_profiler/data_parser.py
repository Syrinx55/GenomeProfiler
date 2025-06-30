# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

from Bio import SeqIO
import re
import os
import pandas as pd
import argparse
from pathlib import Path


def parse_cds_fasta(fasta_file, output_tsv):
    records = SeqIO.parse(fasta_file, "fasta")
    output_lines = []

    for record in records:
        # Expected header format: >accession_gene_start123_stop456
        match = re.search(r"start(\d+)_stop(\d+)", record.id)
        if match:
            start, stop = match.groups()
            gene = ""
            header_parts = record.id.split("_start")[0].split("_")
            if len(header_parts) > 1:
                gene_candidate = header_parts[-1]
                if not gene_candidate.lower().startswith("gene"):
                    gene = gene_candidate
            output_lines.append(f"{start}\t{stop}\t{gene}")
        else:
            print(f"[WARNING] Could not parse coordinates from: {record.id}")

    with open(output_tsv, "w") as out_f:
        for line in output_lines:
            out_f.write(line + "\n")


def parse_abricate_tsv(input_file, output_prefix):
    pos_lines = []
    neg_lines = []

    with open(input_file, "r") as in_f:
        header = in_f.readline().strip().split("\t")
        header_indices = {col: i for i, col in enumerate(header)}
        for line in in_f:
            fields = line.strip().split("\t")
            start = fields[header_indices.get("START", 0)]
            end = fields[header_indices.get("END", 1)]
            gene = fields[header_indices.get("GENE", 2)]
            strand = (
                fields[header_indices.get("STRAND", -1)]
                if "STRAND" in header_indices
                else "+"
            )

            line_out = f"{start}\t{end}\t{gene}"
            if strand == "-":
                neg_lines.append(line_out)
            else:
                pos_lines.append(line_out)

    for strand_label, lines in [("positive", pos_lines), ("negative", neg_lines)]:
        output_file = f"{output_prefix}_{strand_label}.tsv"
        with open(output_file, "w") as out_f:
            for line in lines:
                out_f.write(line + "\n")


def parse_isescan_csv(csv_file, output_file):
    df = pd.read_csv(csv_file)
    if {"isBegin", "isEnd", "cluster"}.issubset(df.columns):
        with open(output_file, "w") as out_f:
            for _, row in df.iterrows():
                out_f.write(f"{row['isBegin']}\t{row['isEnd']}\t{row['cluster']}\n")
    else:
        print(f"[ERROR] Missing required columns in {csv_file}")


def parse_mobileog_tsv(input_file, output_dir):
    df = pd.read_csv(input_file, sep="\t")
    required_cols = {"strand", "sseqid", "nuc_start", "nuc_stop"}
    if not required_cols.issubset(df.columns):
        print(f"[ERROR] Missing required columns in {input_file}")
        return

    # Extract gene name and category from sseqid
    def extract_info(sseqid):
        parts = sseqid.split("|")
        gene = (
            parts[1] if len(parts) > 1 and parts[1] not in {"NA", "NA:Keyword"} else ""
        )
        category = parts[3].lower() if len(parts) > 3 else "unknown"
        return gene, category

    # Normalize categories to fixed 5-category names
    category_mapping = {
        "stability/transfer/defense": "Stability_Transfer_Defense",
        "transfer": "Transfer",
        "replication/recombination/repair": "Replication_Recombination_Repair",
        "phage": "Phage",
        "integration/excision": "Integration_Excision",
    }

    for _, row in df.iterrows():
        gene, raw_category = extract_info(row["sseqid"])
        strand = row["strand"]
        start = row["nuc_start"]
        end = row["nuc_stop"]
        category = category_mapping.get(raw_category, "Other")
        output_path = os.path.join(output_dir, f"mobileog_{strand}_{category}.tsv")
        with open(output_path, "a") as f:
            f.write(f"{start}\t{end}\t{gene}\n")


def parse_islandviewer_tsv(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t")
    if {"Gene start", "Gene end"}.issubset(df.columns):
        with open(output_file, "w") as out_f:
            for _, row in df.iterrows():
                out_f.write(f"{row['Gene start']}\t{row['Gene end']}\n")
    else:
        print(f"[ERROR] Missing required columns in {input_file}")


def parse_tncentral_tsv(input_file, output_dir):
    df = pd.read_csv(input_file, sep="\t", header=None)
    if df.shape[1] < 5:
        print(f"[ERROR] Unexpected column count in {input_file}")
        return

    # Save all feature coordinates
    all_coords_file = os.path.join(output_dir, "tncentral_all.tsv")
    with open(all_coords_file, "w") as out_f:
        for _, row in df.iterrows():
            out_f.write(f"{row[6]}\t{row[7]}\t\n")

    # Filter and save ISFinder-specific results
    isfinder_df = df[df[1].str.startswith("ISFinder_")]
    isfinder_file = os.path.join(output_dir, "tncentral_isfinder.tsv")
    with open(isfinder_file, "w") as out_f:
        for _, row in isfinder_df.iterrows():
            name = row[1].replace("ISFinder_", "")
            out_f.write(f"{row[6]}\t{row[7]}\t{name}\n")


def parse_integron_finder_tsv(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", skiprows=2)
    if {"pos_beg", "pos_end", "annotation"}.issubset(df.columns):
        with open(output_file, "w") as out_f:
            for _, row in df.iterrows():
                out_f.write(
                    f"{row['pos_beg']}\t{row['pos_end']}\t{row['annotation']}\n"
                )
    else:
        print(f"[ERROR] Missing required columns in {input_file}")


def parse_isescan_inverted_repeats_gff(
    gff_file, output_positive, output_negative, output_unstranded
):
    with open(gff_file, "r") as in_f:
        lines = in_f.readlines()

    with open(output_positive, "w") as out_pos, open(
        output_negative, "w"
    ) as out_neg, open(output_unstranded, "w") as out_unstr:
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "terminal_inverted_repeat":
                continue
            start, end = parts[3], parts[4]
            strand1 = parts[6].strip()
            strand2 = parts[7].strip()
            strand = None
            for s in [strand1, strand2]:
                if s in {"+", "-"}:
                    strand = s
                    break
            if strand == "+":
                out_pos.write(f"{start}\t{end}\n")
            elif strand == "-":
                out_neg.write(f"{start}\t{end}\n")
            else:
                out_unstr.write(f"{start}\t{end}\n")


def run_parser(output_path):
    output_path = Path(output_path).expanduser().resolve()

    # Diagnostic output
    print(f"[DEBUG] output_path: {output_path}")

    base_path = output_path
    parsed_dir = base_path / "parsed_output"
    # FIXME changed
    parsed_dir.mkdir(parents=True, exist_ok=True)

    accession = base_path.name.split("_")[0]

    cds_pos_path = base_path / "CDS/positive_strand.fna"
    if cds_pos_path.exists():
        parse_cds_fasta(cds_pos_path, parsed_dir / "positive_cds.tsv")

    cds_neg_path = base_path / "CDS/negative_strand.fna"
    if cds_neg_path.exists():
        parse_cds_fasta(cds_neg_path, parsed_dir / "negative_cds.tsv")

    for db in ["argannot", "card", "plasmidfinder", "resfinder", "vfdb"]:
        abricate_path = base_path / f"abricate/{db}.tsv"
        if abricate_path.exists():
            parse_abricate_tsv(abricate_path, parsed_dir / f"{db}_parsed")

    isescan_csv_path = base_path / "isescan/ncbi" / f"{accession}.fasta.csv"
    if isescan_csv_path.exists():
        parse_isescan_csv(isescan_csv_path, parsed_dir / "isescan_parsed.tsv")

    mobileog_path = base_path / "mobileogdb/mobileOG_integrated.tsv"
    if mobileog_path.exists():
        parse_mobileog_tsv(mobileog_path, parsed_dir)

    islandviewer_path = base_path / "islandviewer" / f"{accession}_islandviewer.tsv"
    if islandviewer_path.exists():
        parse_islandviewer_tsv(
            islandviewer_path, parsed_dir / "islandviewer_parsed.tsv"
        )

    tncentral_path = base_path / "tncentral" / f"{accession}_tncentral_blast.tsv"
    if tncentral_path.exists():
        parse_tncentral_tsv(tncentral_path, parsed_dir)

    integron_finder_path = (
        base_path
        / "integron_finder"
        / f"Results_Integron_Finder_{accession}"
        / f"{accession}.integrons"
    )
    if integron_finder_path.exists():
        parse_integron_finder_tsv(
            integron_finder_path, parsed_dir / "integron_finder_parsed.tsv"
        )

    isescan_gff_path = base_path / "isescan/ncbi" / f"{accession}.fasta.gff"
    if isescan_gff_path.exists():
        parse_isescan_inverted_repeats_gff(
            isescan_gff_path,
            parsed_dir / "isescan_inverted_repeats_positive.tsv",
            parsed_dir / "isescan_inverted_repeats_negative.tsv",
            parsed_dir / "isescan_inverted_repeats_unstranded.tsv",
        )


def main():
    parser = argparse.ArgumentParser(
        description="Parse feature outputs from genome pipeline."
    )
    parser.add_argument(
        "--accession", type=str, help="NCBI accession (e.g. NZ_CP169015)"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config_genomeprofiler_REAL.ini",
        help="Path to config file",
    )
    parser.add_argument(
        "--skip-parsing", action="store_true", help="Skip parsing outputs."
    )
    args = parser.parse_args()

    if args.skip_parsing:
        print("[INFO] Skipping data parsing.")
        return

    if not args.accession:
        print("[ERROR] No accession provided. Parsing cannot proceed.")
        return

    run_parser(args.accession, args.config)


if __name__ == "__main__":
    main()
