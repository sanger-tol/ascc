#!/usr/bin/env python3
"""
Script conversion of .sam file with mapped reads to sorted and indexed .bam file
"""
# MIT License
#
# Copyright (c) 2020-2021 Genome Research Ltd.
#
# Author: Eerik Aunin (ea10@sanger.ac.uk)
#
# This file is a part of the Genome Decomposition Analysis (GDA) pipeline.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import general_purpose_functions as gpf
import argparse


def main(sam_file, threads, fasta_path):
    sam_file_extensionless_path = sam_file.split(".sam")[0]
    bam_file = sam_file_extensionless_path + ".bam"
    sorted_bam_file = sam_file_extensionless_path + "_sorted.bam"
    sam_to_bam_command = None
    if fasta_path == "":
        sam_to_bam_command = "samtools view -Sb -@ {} {} > {}".format(threads, sam_file, bam_file)
    else:
        sam_to_bam_command = "samtools view -T {} -Sb -@ {} {} > {}".format(fasta_path, threads, sam_file, bam_file)
    gpf.run_system_command(sam_to_bam_command)
    sort_bam_command = "samtools sort -@ {} {} -o {}".format(threads, bam_file, sorted_bam_file)
    gpf.run_system_command(sort_bam_command)
    index_bam_command = "samtools index " + sorted_bam_file
    gpf.run_system_command(index_bam_command)
    remove_sam_command = "rm " + sam_file
    gpf.run_system_command(remove_sam_command)
    remove_unsorted_bam_command = "rm " + bam_file
    gpf.run_system_command(remove_unsorted_bam_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sam_file", type=str, help="Path to input SAM file")
    parser.add_argument("--threads", type=int, help="Number of threads (default: 1)", default=1)
    parser.add_argument(
        "--fasta_path",
        type=str,
        help="Optional (needed when the reads have been mapped to a genome larger than 4Gb with minimap2): path the FASTA file to which the reads were mapped",
        default="",
    )
    args = parser.parse_args()
    main(args.sam_file, args.threads, args.fasta_path)
