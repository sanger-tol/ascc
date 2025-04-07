#!/usr/bin/env python3

import os
import sys
import argparse
import textwrap

VERSION = "1.0.0"

DESCRIPTION = """
Script to check if the nt BLAST database has taxonomy included by verifying the presence of
taxonomy files (taxdb.btd, taxonomy4blast.sqlite3, taxdb.bti) in the database directory.

This check is specifically for the nt BLAST database used in the EXTRACT_NT_BLAST subworkflow,
not for other BLAST databases used elsewhere in the pipeline (VecScreen, PacBio barcodes check, etc.).

Version: {VERSION}
Written by Eerik Aunin (eeaunin)
"""


def check_nt_blast_taxonomy(db_path):
    """
    Check if the nt BLAST database at the given path has taxonomy included.
    Returns True if all required taxonomy files are present, False otherwise.
    """
    taxonomy_files = ["taxdb.btd", "taxdb.bti"]
    optional_files = ["taxonomy4blast.sqlite3"]  # This file might not be present in all installations

    missing_files = []
    for file in taxonomy_files:
        if not os.path.exists(os.path.join(db_path, file)):
            missing_files.append(file)

    if missing_files:
        sys.stderr.write(f"ERROR: The nt BLAST database at '{db_path}' does not have taxonomy included.\n")
        sys.stderr.write(f"The following required taxonomy files are missing: {', '.join(missing_files)}\n")
        sys.stderr.write(
            "Please rebuild the nt BLAST database with taxonomy included using the -parse_seqids and -taxid_map options with makeblastdb,\n"
        )
        sys.stderr.write(
            "or download the pre-built nt BLAST database from NCBI which already includes taxonomy information.\n"
        )
        sys.stderr.write(
            "This is required for the nt_blast step since the pipeline no longer uses accession2taxid files.\n"
        )
        return False

    # Check for optional files and warn if missing
    missing_optional = []
    for file in optional_files:
        if not os.path.exists(os.path.join(db_path, file)):
            missing_optional.append(file)

    if missing_optional:
        sys.stderr.write(
            f"WARNING: The following optional taxonomy files are missing from the nt BLAST database: {', '.join(missing_optional)}\n"
        )
        sys.stderr.write(
            "The database may still work, but it's recommended to have all taxonomy files for optimal performance.\n"
        )

    sys.stderr.write("nt BLAST database taxonomy files found. The database has taxonomy included.\n")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Check NT BLAST Taxonomy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(DESCRIPTION),
    )
    parser.add_argument("--db_path", required=True, help="Path to the nt BLAST database directory")
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    args = parser.parse_args()

    if not check_nt_blast_taxonomy(args.db_path):
        # Print status to stdout for Nextflow to capture
        print("nt_database_taxonomy_files_not_found")
        sys.exit(0)  # Exit normally
    else:
        # Print success status to stdout for Nextflow to capture
        print("nt_database_taxonomy_files_found")
