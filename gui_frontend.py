# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent

import tkinter as tk
from tkinter import messagebox, scrolledtext, filedialog, ttk
from threading import Thread
import os
import json
from pathlib import Path
from collection_pipeline import process_accession
from Bio import Entrez
from data_parser import run_parser
from dotenv import load_dotenv

CONFIG_FILE = "config_genomeprofiler.ini"

SETTINGS_FILE = "user_settings.json"

TOOL_OPTIONS = [
    "abricate",
    "isescan",
    "mobileog",
    "ectyper",
    "tncentral",
    "integron finder",
    "islandviewer",
    "plsdb",
    "phastest",
]


def load_user_settings():
    if os.path.exists(SETTINGS_FILE):
        with open(SETTINGS_FILE, "r") as f:
            return json.load(f)
    return {}


def save_user_settings(settings):
    with open(SETTINGS_FILE, "w") as f:
        json.dump(settings, f)


def launch_pipeline(
    accessions,
    included_tools,
    workers,
    debug_textbox,
    fasta_path,
    genbank_path,
    progress,
    timestamp_output,
    parser_var,
    config,
):
    Entrez.email = os.environ["GENPROF_ENTREZ_EMAIL"]

    ascii_art = r"""

   ______                                ____             _____ __         
  / ____/__  ____  ____  ____ ___  ___  / __ \_________  / __(_) /__  _____
 / / __/ _ \/ __ \/ __ \/ __ `__ \/ _ \/ /_/ / ___/ __ \/ /_/ / / _ \/ ___/
/ /_/ /  __/ / / / /_/ / / / / / /  __/ ____/ /  / /_/ / __/ / /  __/ /    
\____/\___/_/ /_/\____/_/ /_/ /_/\___/_/   /_/   \____/_/ /_/_/\___/_/     
                                                                    
"""
    debug_textbox.insert("end", ascii_art + "\n")
    debug_textbox.see("end")
    debug_textbox.insert(tk.END, f"Launching pipeline on: {accessions}\n")
    debug_textbox.insert(tk.END, f"Included tools: {included_tools}\n")
    debug_textbox.insert(tk.END, f"Using {workers} worker(s)\n")
    debug_textbox.insert(tk.END, "\n[INFO] Running pipeline...\n\n")
    debug_textbox.see(tk.END)

    try:
        for acc in accessions:
            process_accession(
                acc,
                config,
                included_tools=included_tools,
                debug_textbox=debug_textbox,
                fasta_override=fasta_path,
                genbank_override=genbank_path,
                timestamp_output=timestamp_output,
            )
            try:
                if parser_var.get():
                    run_parser(Path(config["output_base"]) / acc)
            except Exception as parser_err:
                debug_textbox.insert(
                    tk.END, f"[ERROR] Parsing failed: {str(parser_err)}\n"
                )
                debug_textbox.see(tk.END)
    except Exception as e:
        debug_textbox.insert(tk.END, f"\n[ERROR] {str(e)}\n")
    debug_textbox.insert(tk.END, "\n[INFO] Done.\n")
    debug_textbox.see(tk.END)
    progress.stop()
    progress.pack_forget()


def get_local_identifier_from_fasta(path):
    try:
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    return line[1:].strip().split()[0]
    except Exception:
        pass
    return "local"


def start_execution(
    entry, checkboxes, worker_spin, fasta_path, genbank_path, timestamp_var, parser_var
):
    accessions = entry.get().strip()
    selected_tools = [tool for tool, var in checkboxes.items() if var.get() == 1]
    workers = int(worker_spin.get())

    if not accessions and not (fasta_path.get() or genbank_path.get()):
        messagebox.showwarning(
            "Input Required", "Please enter accessions or upload at least one file."
        )
        return

    settings = {"included_tools": selected_tools, "workers": workers}
    save_user_settings(settings)

    debug_win = tk.Toplevel()
    debug_win.title("GenomeProfiler Debug Console")
    debug_textbox = scrolledtext.ScrolledText(
        debug_win, width=100, height=30, bg="black", fg="lime"
    )
    debug_textbox.pack(fill="both", expand=True)

    progress = ttk.Progressbar(
        debug_win, orient="horizontal", length=300, mode="indeterminate"
    )
    progress.pack(pady=5)
    progress.start()

    accession_list = (
        accessions.split(",")
        if accessions
        else [get_local_identifier_from_fasta(fasta_path.get())]
    )

    thread = Thread(
        target=launch_pipeline,
        args=(
            accession_list,
            selected_tools,
            workers,
            debug_textbox,
            fasta_path.get() or None,
            genbank_path.get() or None,
            progress,
            timestamp_var.get(),
            parser_var,
        ),
    )
    thread.daemon = True
    thread.start()


def browse_file(var, label):
    path = filedialog.askopenfilename()
    if path:
        var.set(path)
        label.config(text=os.path.basename(path))


def main():
    root = tk.Tk()
    root.title("GenomeProfiler Pipeline GUI")

    user_settings = load_user_settings()

    tk.Label(root, text="Enter comma-separated NCBI Accessions:").pack(pady=5)
    accession_entry = tk.Entry(root, width=60)
    accession_entry.pack(pady=5)

    tk.Label(root, text="Upload Local FASTA File (Optional):").pack()
    fasta_path = tk.StringVar()
    fasta_label = tk.Label(root, text="None Selected")
    fasta_label.pack()
    tk.Button(
        root, text="Browse FASTA", command=lambda: browse_file(fasta_path, fasta_label)
    ).pack()

    tk.Label(root, text="Upload Local GenBank File (Optional):").pack()
    genbank_path = tk.StringVar()
    genbank_label = tk.Label(root, text="None Selected")
    genbank_label.pack()
    tk.Button(
        root,
        text="Browse GenBank",
        command=lambda: browse_file(genbank_path, genbank_label),
    ).pack()

    tk.Label(root, text="Select Tools to INCLUDE in the Pipeline:").pack(pady=5)
    checkboxes = {}
    tool_frame = tk.Frame(root)
    tool_frame.pack(pady=5)
    for i, tool in enumerate(TOOL_OPTIONS):
        var = tk.IntVar(
            value=1 if tool in user_settings.get("included_tools", TOOL_OPTIONS) else 0
        )
        chk = tk.Checkbutton(tool_frame, text=tool, variable=var)
        chk.grid(row=i // 4, column=i % 4, sticky="w", padx=5, pady=2)
        checkboxes[tool] = var

    tk.Label(root, text="Number of Workers:").pack(pady=5)
    worker_spin = tk.Spinbox(root, from_=1, to=8, width=5)
    worker_spin.delete(0, "end")
    worker_spin.insert(0, user_settings.get("workers", 2))
    worker_spin.pack(pady=5)

    parser_var = tk.IntVar()
    tk.Checkbutton(root, text="Parse data for visualization", variable=parser_var).pack(
        pady=5
    )

    timestamp_var = tk.IntVar()
    tk.Checkbutton(
        root, text="Use Timestamped Output Directory", variable=timestamp_var
    ).pack(pady=5)

    tk.Button(
        root,
        text="Run Pipeline",
        command=lambda: start_execution(
            accession_entry,
            checkboxes,
            worker_spin,
            fasta_path,
            genbank_path,
            timestamp_var,
            parser_var,
        ),
    ).pack(pady=10)

    root.mainloop()


if __name__ == "__main__":
    main()
