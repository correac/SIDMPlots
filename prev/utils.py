"""
Some helper functions .

Contains the argument parser and default parsing results.

See the README for available argument parsers.
"""

import argparse as ap
import glob
from typing import Optional

parser = ap.ArgumentParser(
    description="""General argument parser for SIDMPlots scripts."""
)

parser.add_argument(
    "-d",
    "--directory",
    help="Directory containing snapshots. Required.",
    type=str,
    required=True,
    nargs="*",
)

parser.add_argument(
    "-s",
    "--snapshot",
    help="Snapshot number to visualise. Required.",
    type=str,
    required=True,
    nargs="*",
)

parser.add_argument(
    "-n",
    "--run-names",
    help="Names of the runs for placement in legends.",
    type=str,
    required=False,
    nargs="*",
)

parser.add_argument(
    "-o",
    "--output",
    help="Output directory for the figures. If not present, the same directory as the snapshots are in is used.",
    required=False,
    type=str,
    default=None,
)

parser.add_argument(
    "-t",
    "--type",
    help="Type of simulation, Hydro or DMONLY. Default DMONLY.",
    required=False,
    type=str,
    default=None,
)

args = parser.parse_args()

if args.type is None:
    args.type = "DMONLY"

if args.output is None:
    args.output = args.directory[0]

