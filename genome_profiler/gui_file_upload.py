# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

from tkinter import filedialog


def select_fasta_file():
    path = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*")],
    )
    if path:
        print(f"[INFO] FASTA file selected: {path}")
    return path


def select_genbank_file():
    path = filedialog.askopenfilename(
        title="Select GenBank File",
        filetypes=[("GenBank files", "*.gbk *.gb"), ("All files", "*")],
    )
    if path:
        print(f"[INFO] GenBank file selected: {path}")
    return path


def get_fasta_path(current_path):
    return current_path


def get_genbank_path(current_path):
    return current_path


def get_identifier_from_fasta_header(path):
    try:
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    return line[1:].strip().split()[0]
    except Exception:
        pass
    return "local"
