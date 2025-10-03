#!/usr/bin/env python3
"""
General purpose functions
File for functions that can be reused in many Python scripts
"""
# MIT License
#
# Copyright (c) 2020-2021 Genome Research Ltd.
#
# Author: Eerik Aunin (eeaunin@gmail.com)
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

import os
from os.path import isfile
import sys
import subprocess
import signal
from datetime import datetime
import argparse


def l(path):
    """
    Function for loading text file as a list and removing line breaks from line ends
    """
    lines = []
    if isfile(path):
        with open(path, "r") as in_file:
            lines = in_file.readlines()
            lines = [x.rstrip() for x in lines]
    else:
        sys.stderr.write("Error: file not found (" + path + ")\n")
        sys.exit(1)
    return lines


def ll(in_path):
    """
    Function for reading a text file line by line
    """
    if isfile(in_path):
        with open(in_path, "r") as in_file:
            for line in in_file:
                line = line.rstrip()
                yield line
    else:
        sys.stderr.write("Error: file not found (" + in_path + ")\n")
        sys.exit(1)


def export_list_as_line_break_separated_file(out_list_f, out_path):
    """
    Exports a list to a file, each item on a separate row
    """
    with open(out_path, "w") as out_file:
        for item in out_list_f:
            out_file.write(str(item))
            out_file.write("\n")


def spl(line, left_splitter, right_splitter, direction=0):
    """
    Function for cropping a string from left and right side
    Direction:  if 0: the string will be cropped first from left and then right
                if 1: the string will be cropped first from right and then left
    Returns None if the splitting cannot be done because a splitter or both splitters are not in the input string
    """
    out_line = None
    if left_splitter in line and right_splitter in line:
        if direction == 0:
            out_line = line.split(left_splitter)[1]
            out_line = out_line.split(right_splitter)[0]
        elif direction == 1:
            out_line = line.split(right_splitter)[0]
            out_line = out_line.split(left_splitter)[1]
    return out_line


def print_with_fixed_row_length(seq, max_length):
    """
    Input: 1) a string 2) maximum line length in output
    Output: the input string printed to STDOUT in chunks with fixed maximum line length
    """
    split_seq = [seq[i : i + max_length] for i in range(0, len(seq), max_length)]
    for line in split_seq:
        print(line)


def split_with_fixed_row_length(seq, max_length):
    """
    Input: 1) a string 2) maximum line length in output
    Output: the input string split in chunks with fixed maximum line length
    """
    split_seq = [seq[i : i + max_length] for i in range(0, len(seq), max_length)]
    return split_seq


def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence
    """
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n",
    }
    reverse_comp = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_comp


def read_fasta_in_chunks(in_path):
    """
    Input: path to FASTA file
    Output (iterator): string tuples where the first element is a FASTA header and the second element is the corresponding FASTA sequence
    """
    in_data = ll(in_path)
    current_seq_header = None
    seq = ""
    for line in in_data:
        if line != "":
            if line[0] == ">":
                if seq != "":
                    yield (current_seq_header, seq)
                seq = ""
                current_seq_header = line[1 : len(line)]
            else:
                seq += line
    if seq != "":
        yield (current_seq_header, seq)


def list_to_chunks(lst, n):
    """
    Yield successive n-sized chunks from lst
    https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def string_to_chunks(line, n):
    """
    Function for splitting a string every nth character
    https://stackoverflow.com/questions/9475241/split-string-every-nth-character
    """
    return [line[i : i + n] for i in range(0, len(line), n)]


def run_system_command(
    system_command, verbose=True, dry_run=False, tries=1, expected_exit_code=0
):
    """
    Executes a system command and checks its exit code
    """
    triggering_script_name = sys.argv[0].split("/")[-1]
    try_counter_string = ""
    if dry_run == False:
        for i in range(0, tries):
            if verbose == True:
                time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                if i > 0:
                    try_counter_string = ", try {}".format(i + 1)
                out_message = "<{}, {}{}> executing command: {}\n".format(
                    time_now, triggering_script_name, try_counter_string, system_command
                )
                sys.stderr.write(out_message)
            try:
                output = subprocess.check_output(
                    system_command,
                    stderr=subprocess.STDOUT,
                    shell=True,
                    timeout=None,
                    universal_newlines=True,
                )
                break
            except subprocess.CalledProcessError as exc:
                out_errormessage = (
                    "<"
                    + triggering_script_name
                    + "> "
                    + " exited with error code "
                    + str(exc.returncode)
                )
                if exc.output.isspace() == False:
                    out_errormessage += ". Error message: " + exc.output
                if i == tries - 1:
                    if exc.returncode != expected_exit_code:
                        sys.stderr.write(out_errormessage + "\n")
                        os.kill(os.getpid(), signal.SIGINT)


def check_if_file_exists(in_path):
    """
    Function for checking if a file exists
    """
    if os.path.isfile(in_path) == False:
        sys.stderr.write("The input file " + in_path + " was not found\n")
        sys.exit(1)


def get_file_paths(in_folder_path, extension):
    """
    Function for getting the paths to all files with a specific extension in a user-specified folder
    in_folder_path: path to the folder with input files
    extension: file extension of input files
    Output: paths to individual files with the specific extension (list)
    """
    onlyfiles = list()
    selected_file_paths = list()
    if os.path.isdir(in_folder_path):
        onlyfiles = [
            f
            for f in os.listdir(in_folder_path)
            if os.path.isfile(os.path.join(in_folder_path, f))
        ]
        for file_item in onlyfiles:
            if "." + extension in file_item:
                file_item_split = file_item.split(".")
                if file_item_split[len(file_item_split) - 1] == extension:
                    selected_file_paths.append(in_folder_path + "/" + file_item)
    else:
        sys.stderr.write("Error: folder not found (" + in_folder_path + ")\n")
        sys.exit(1)
    return selected_file_paths


def main():
    # Placeholder for accessing version.
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
