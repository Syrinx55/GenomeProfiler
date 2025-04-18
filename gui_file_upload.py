# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import tkinter as tk
from tkinter import filedialog

_fasta_path = None
_genbank_path = None


def select_fasta_file():
    global _fasta_path
    _fasta_path = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*")],
    )
    if _fasta_path:
        print(f"[INFO] FASTA file selected: {_fasta_path}")


def select_genbank_file():
    global _genbank_path
    _genbank_path = filedialog.askopenfilename(
        title="Select GenBank File",
        filetypes=[("GenBank files", "*.gbk *.gb"), ("All files", "*")],
    )
    if _genbank_path:
        print(f"[INFO] GenBank file selected: {_genbank_path}")


def get_fasta_path():
    return _fasta_path


def get_genbank_path():
    return _genbank_path
