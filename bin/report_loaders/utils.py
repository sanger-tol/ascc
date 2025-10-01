#!/usr/bin/env python3
"""
Utility functions for ASCC HTML report generation.

This module contains utility functions used by other report loader modules.
Developed by Eerik Aunin (eeaunin@gmail.com)
"""
import os


def wrap_table_html(table_html):
    """Wrap HTML table with responsive container layers.

    Args:
        table_html (str): The HTML table string

    Returns:
        str: The wrapped HTML table
    """
    return f"""
    <div class="outer-container">
        <div class="table-responsive">
            <div class="table-wrapper">
                {table_html}
            </div>
        </div>
    </div>
    """


def find_files_in_dir(directory, pattern=None, extension=None):
    """Find files in a directory matching a pattern or extension.

    Args:
        directory (str): The directory to search in
        pattern (str, optional): A substring to match in filenames
        extension (str, optional): A file extension to match (e.g., '.txt')

    Returns:
        list: A list of full paths to matching files
    """
    if not directory or not os.path.exists(directory) or not os.path.isdir(directory):
        return []

    files = []
    for f in os.listdir(directory):
        file_path = os.path.join(directory, f)
        if os.path.isfile(file_path):
            if extension and f.endswith(extension):
                files.append(file_path)
            elif pattern and pattern in f:
                files.append(file_path)
            elif not extension and not pattern:
                files.append(file_path)

    return files
