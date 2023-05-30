#!/usr/bin/env python3
import gzip
import json
from pathlib import Path
import pysam
from typing import Dict


def run():

    lengths = parse_fasta_dict()

    for vcf in Path("input").rglob("*.vcf.gz"):
        filter_variants(vcf, lengths)


def parse_fasta_dict() -> Dict[str, int]:
    # Parse the fasta dict
    print("Reading in fasta.dict")
    with pysam.AlignmentFile("fasta.dict", "r") as handle:
        lengths = {
            contig["SN"]: contig["LN"]
            for contig in handle.header.to_dict()["SQ"]
        }

    print("Read in contig lengths:")
    print(json.dumps(lengths, indent=4))
    return lengths


def filter_variants(vcf: Path, lengths: Dict[str, int]):

    chrs = set()
    nvars = 0

    print(f"Reading in {vcf}")

    # Direct the output to the current working directory
    print(f"Writing out to {vcf.name}")
    with gzip.open(vcf, 'rt') as input:
        with gzip.open(vcf.name, 'wt') as output:

            for line in input:

                if line.startswith("##contig"):
                    if get_chrname_header(line) in lengths:
                        output.write(line)
                        print(line)
                elif line.startswith("#"):
                    output.write(line)
                else:
                    chrname = get_chrname_var(line)
                    if chrname in lengths:
                        output.write(line)
                        nvars += 1
                        chrs.add(chrname)

    print(f"Found {nvars:,} variants from {len(chrs)} contigs")


def get_chrname_var(line):
    return line.split("\t", 1)[0]


def get_chrname_header(line: str):
    chrname = line.split("ID=", 1)[1]
    chrname = chrname.split(",", 1)[0]
    return chrname


if __name__ == "__main__":
    run()
