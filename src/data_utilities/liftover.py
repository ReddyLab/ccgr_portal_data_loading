#!/usr/bin/env python

import csv
import argparse

from pyliftre import Mapper


def merge_mapping(mapping):
    assert len(mapping) > 0
    new_mapping = list(mapping[0])
    for region in mapping[1:]:
        if new_mapping[2] - region[1] > 100:
            raise ValueError("big gap")
        new_mapping[2] = region[2]

    return new_mapping


class Experiment:
    def __init__(self, experiment_file, output, unmapped_file, mapper):
        self.experiment_file = csv.DictReader(experiment_file, delimiter="\t")
        self.writer = csv.DictWriter(
            output,
            fieldnames=(
                "chrom",
                "start",
                "end",
                "strand",
                "parent_chrom",
                "parent_start",
                "parent_end",
                "parent_strand",
                "facets",
                "misc",
            ),
            delimiter="\t",
        )
        self.writer.writeheader()
        self.unmapped_file = unmapped_file
        self.mapper = mapper

    def liftover(self):
        for row in self.experiment_file:
            strand = row["strand"]
            if strand == ".":
                strand = "+"
            new_coords = self.mapper.map_coordinates(row["chrom"], int(row["start"]), int(row["end"]), strand, "l")[
                1::2
            ]
            new_row = {
                "facets": row["facets"],
                "misc": row["misc"],
            }

            match len(new_coords):
                case 0:
                    if self.unmapped_file:
                        self.unmapped_file.write(f"base\t{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}")
                case _:
                    new_coords = merge_mapping(new_coords)
                    new_row.update(
                        {
                            "chrom": new_coords[0],
                            "start": new_coords[1],
                            "end": new_coords[2],
                            "strand": new_coords[3] if row["strand"] != "." else ".",
                        }
                    )

            p_strand = row["parent_strand"]
            if p_strand == ".":
                p_strand = "+"
            new_parent_coords = self.mapper.map_coordinates(
                row["parent_chrom"], int(row["parent_start"]), int(row["parent_end"]), p_strand, "l"
            )[1::2]
            match len(new_parent_coords):
                case 0:
                    if self.unmapped_file:
                        self.unmapped_file.write(
                            f"parent\t{row['parent_chrom']}\t{row['parent_start']}\t{row['parent_end']}\t{p_strand}"
                        )
                case _:
                    new_parent_coords = merge_mapping(new_parent_coords)
                    new_row.update(
                        {
                            "parent_chrom": new_parent_coords[0],
                            "parent_start": new_parent_coords[1],
                            "parent_end": new_parent_coords[2],
                            "parent_strand": new_parent_coords[3] if row["parent_strand"] != "." else ".",
                        }
                    )

            if "chrom" in new_row and "parent_chrom" in new_row:
                self.writer.writerow(new_row)


class Analysis:
    def __init__(self, analysis_file, output, unmapped_file, mapper):
        self.analysis_file = csv.DictReader(analysis_file, delimiter="\t")
        self.writer = csv.DictWriter(
            output,
            fieldnames=(
                "chrom",
                "start",
                "end",
                "strand",
                "gene_name",
                "gene_ensembl_id",
                "raw_p_val",
                "adj_p_val",
                "effect_size",
                "facets",
            ),
            delimiter="\t",
        )
        self.writer.writeheader()
        self.unmapped_file = unmapped_file
        self.mapper = mapper

    def liftover(self):
        for row in self.analysis_file:
            strand = row["strand"]
            if strand == ".":
                strand = "+"
            new_coords = self.mapper.map_coordinates(row["chrom"], int(row["start"]), int(row["end"]), strand, "l")[
                1::2
            ]

            match len(new_coords):
                case 0:
                    if self.unmapped_file:
                        self.unmapped_file.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}")
                case _:
                    try:
                        new_coords = merge_mapping(new_coords)
                    except:
                        n = [(f"{c[0]}:{c[1]}-{c[2]}:{c[3]} ({c[2] - c[1]})") for c in new_coords]
                        print(
                            f"{row['chrom']}:{row['start']}-{row['end']}:{strand} ({int(row['end']) - int(row['start'])}) -> {n}"
                        )
                        continue
                    self.writer.writerow(
                        {
                            "chrom": new_coords[0],
                            "start": new_coords[1],
                            "end": new_coords[2],
                            "strand": new_coords[3] if row["strand"] != "." else ".",
                            "gene_name": row["gene_name"],
                            "gene_ensembl_id": row["gene_ensembl_id"],
                            "raw_p_val": row["raw_p_val"],
                            "adj_p_val": row["adj_p_val"],
                            "effect_size": row["effect_size"],
                            "facets": row["facets"],
                        }
                    )


def get_args():
    parser = argparse.ArgumentParser("A tool for lifting over experiment and analysis data")
    parser.add_argument("-c", "--chain", required=True, help="The chain file containing the liftover data")
    parser.add_argument("-u", "--unmapped", type=argparse.FileType("w"), help="The file for logging unmapped regions")
    parser.add_argument(
        "-o", "--output", type=argparse.FileType("w"), help="The output file with the lifted over regions"
    )

    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("-e", "--experiment", type=argparse.FileType(), help="The experiment data file to lift over")
    inputs.add_argument("-a", "--analysis", type=argparse.FileType(), help="The analysis data file to lift over")

    return parser.parse_args()


def run(args):
    mapper = Mapper(args.chain)
    if args.experiment:
        regions = Experiment(args.experiment, args.output, args.unmapped, mapper)
    elif args.analysis:
        regions = Analysis(args.analysis, args.output, args.unmapped, mapper)

    regions.liftover()


if __name__ == "__main__":
    run(get_args())
