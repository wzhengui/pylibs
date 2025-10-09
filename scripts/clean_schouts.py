#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Script Name: clean_schouts.py
#
# Description:
#   clear all the output files in the "output" directory
#
# Usage:
#   ./clean_schouts.py
#
# Author: Wenfan Wu, 10/09/2025
# -----------------------------------------------------------------------------
import os
import glob

## define necessary directories
workdir = os.getenv('PWD')  # main work directory
schout = os.path.join(workdir, "outputs")  # directory for schism outputs

## delete old schouts
ok = input("\nDo you wish to continue? (Y/n): ")

if ok.lower() == 'y':
    os.chdir(schout)
    
    files_to_delete = [
        "schout*", "*.nc", "*.out", "*.scribe", "*fatal*", 
        "max*", "*global*", "*.out.nml"
    ]

    for pattern in files_to_delete:
        for file in glob.glob(pattern):
            try:
                os.remove(file)
                print(f"Deleted: {file}")
            except FileNotFoundError:
                print(f"No {pattern} files found")
            except Exception as e:
                print(f"Error deleting {file}: {e}")

elif ok.lower() == 'n':
    print("Stop cleaning SCHOUTS")
else:
    print("Please input y/n")

