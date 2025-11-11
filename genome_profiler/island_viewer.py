# Copyright (c) 2025 Jacob Alford <jalford0000@gmail.com>
# SPDX-License-Identifier: BSD-2-Clause-Patent
# pylint: disable=undefined-variable

import os
import requests
import sys
from requests_toolbelt.multipart.encoder import MultipartEncoder


def submit_islandviewer_job(genbank_path, config):
    """
    Submits a genome file (in GenBank format) to IslandViewer using the HTTP API.

    The function sources the API token, submission URL, and email from the config.

    Expected config keys:
      - islandviewer_api_submit: The submission URL (e.g., "https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/")
      - islandviewer_auth_token: Your IslandViewer API authentication token.  This token can be created by visiting the site (https://www.pathogenomics.sfu.ca/islandviewer/) and creating an account, then generating a token.
      - entrez_email: Your email address.

    Parameters:
      genbank_path (str): Path to the GenBank file.
      config (SectionProxy): Configuration object (from ConfigParser).

    Returns:
      dict: Parsed JSON response from IslandViewer.
    """
    # Get values from the config file
    email = os.environ["GENPROF_ENTREZ_EMAIL"]
    token = os.environ["GENPROF_ISLANDVIEWER_AUTH_TOKEN"]

    if not email or not token:
        raise ValueError(
            "Missing required config parameters for IslandViewer submission."
        )

    submit_url = server = "https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/"
    ext = "/rest/submit/"

    multipart_data = MultipartEncoder(
        fields={
            "format_type": "GENBANK",  # Use "FASTA" if preferred
            "email_addr": email,
            "genome_file": ("filename", open(genbank_path, "rb"), "text/plain"),
        }
    )

    headers = {"Content-Type": multipart_data.content_type, "x-authtoken": token}

    print(f"[INFO] Submitting IslandViewer job for {genbank_path}...")
    r = requests.post(submit_url, headers=headers, data=multipart_data)

    # Debug: print raw response for troubleshooting
    print("[DEBUG] Raw response text:")
    print(r.text)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    try:
        decoded = r.json()
        print("[DEBUG] Parsed response:", repr(decoded))
    except ValueError as e:
        print(f"[ERROR] JSON decoding failed: {e}")
        sys.exit(1)

    return decoded


if __name__ == "__main__":
    config = {
        "islandviewer_api_submit": "https://www.pathogenomics.sfu.ca/islandviewer/rest/submit/",
        "islandviewer_auth_token": "your_authentication_token",
    }
    genbank_file = genbank_path = os.path.join(dirs["ncbi"], f"{accession}.gbk")
    # Replace with the actual path to your genome file
    response = submit_islandviewer_job(genbank_file, config)
    print(response)
