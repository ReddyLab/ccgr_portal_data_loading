#!/usr/bin/env python

import csv
import argparse

from collections import defaultdict
from importlib.resources import files

from pyliftre import Mapper


def merge_mapping(mapping):
    assert len(mapping) > 0
    new_mapping = list(mapping[0])
    for region in mapping[1:]:
        if new_mapping[2] - region[1] > 100 or new_mapping[3] != region[3]:
            raise ValueError("big gap")
        new_mapping[2] = region[2]

    return new_mapping


class Liftover:
    fieldnames = ()

    def __init__(self, liftover_file, output, unmapped_file, summary_file, mapper, **kwargs):
        self.liftover_file = csv.DictReader(liftover_file, delimiter="\t")
        self.writer = csv.DictWriter(
            output,
            fieldnames=self.fieldnames,
            delimiter="\t",
        )
        self.writer.writeheader()
        self.unmapped_file = unmapped_file
        self.summary_file = summary_file
        self.mapper = mapper
        self.mapped_count = 0
        self.unmapped_count = 0
        self.unmapped_genes = None

    def _handle_unmapped(self, line=None):
        self.unmapped_count += 1
        if line is not None:
            if self.unmapped_file:
                self.unmapped_file.write(line)
            else:
                print(line)

    def _write_mapped(self, data):
        self.writer.writerow(data)
        self.mapped_count += 1

    def write_summary(self):
        if self.summary_file:
            summary = f"Mapped Coordinates: {self.mapped_count}\nUnmapped Coordinates: {self.unmapped_count}\n"
            if self.unmapped_genes is not None:
                summary += f"Unmapped Genes: {self.unmapped_genes}\n"

            self.summary_file.write(summary)


class Experiment(Liftover):
    fieldnames = (
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
    )

    def liftover(self):
        for row in self.liftover_file:
            strand = row["strand"]
            if strand == ".":
                strand = "+"

            new_coords = self.mapper.map_coordinates(row["chrom"], int(row["start"]), int(row["end"]), strand, "l")
            if new_coords is None:
                self._handle_unmapped(f"base\t{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                continue
            new_coords = new_coords[1::2]

            new_row = {
                "facets": row["facets"],
                "misc": row["misc"],
            }

            match len(new_coords):
                case 0:
                    self._handle_unmapped(f"base\t{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                case _:
                    try:
                        new_coords = merge_mapping(new_coords)
                    except:
                        self._handle_unmapped(f"base\t{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                        continue

                    # Chromosome names longer than 10 chars aren't allowed in the database.
                    # It's not technically unmapped, but may as well be
                    if len(new_coords[0]) > 10:
                        self._handle_unmapped(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                        continue

                    new_row.update(
                        {
                            "chrom": new_coords[0],
                            "start": new_coords[1],
                            "end": new_coords[2],
                            "strand": new_coords[3] if row["strand"] != "." else ".",
                        }
                    )

            no_parent = False
            if row["parent_chrom"] != "":
                p_strand = row["parent_strand"]
                if p_strand == ".":
                    p_strand = "+"
                new_parent_coords = self.mapper.map_coordinates(
                    row["parent_chrom"], int(row["parent_start"]), int(row["parent_end"]), p_strand, "l"
                )

                if new_parent_coords is None:
                    self._handle_unmapped(
                        f"parent\t{row['parent_chrom']}\t{row['parent_start']}\t{row['parent_end']}\t{p_strand}\n"
                    )
                    continue

                new_parent_coords = new_parent_coords[1::2]
                match len(new_parent_coords):
                    case 0:
                        self._handle_unmapped(
                            f"parent\t{row['parent_chrom']}\t{row['parent_start']}\t{row['parent_end']}\t{p_strand}\n"
                        )
                    case _:
                        try:
                            new_parent_coords = merge_mapping(new_parent_coords)
                        except:
                            self._handle_unmapped(
                                f"parent\t{row['parent_chrom']}\t{row['parent_start']}\t{row['parent_end']}\t{p_strand}\n"
                            )
                            continue

                        # Chromosome names longer than 10 chars aren't allowed in the database.
                        # It's not technically unmapped, but may as well be
                        if len(new_parent_coords[0]) > 10:
                            self._handle_unmapped(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                            continue
                        new_row.update(
                            {
                                "parent_chrom": new_parent_coords[0],
                                "parent_start": new_parent_coords[1],
                                "parent_end": new_parent_coords[2],
                                "parent_strand": new_parent_coords[3] if row["parent_strand"] != "." else ".",
                            }
                        )
            else:
                no_parent = True

            if "chrom" in new_row and (no_parent or "parent_chrom" in new_row):
                self._write_mapped(new_row)

        self.write_summary()


class Analysis(Liftover):
    fieldnames = (
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
    )

    def __init__(self, liftover_file, output, unmapped_file, summary_file, mapper, **kwargs):
        super().__init__(liftover_file, output, unmapped_file, summary_file, mapper, **kwargs)
        self.hg38_genes = defaultdict(set)
        self.hg19_genes = {}
        for gene in kwargs["genes"].readlines():
            hg38_gene, ensembl, hg19_gene = gene.strip().split("\t")
            self.hg38_genes[hg38_gene].add(ensembl)
            if hg19_gene != "NULL":
                self.hg19_genes[(hg19_gene, ensembl)] = hg38_gene
        self.unmapped_genes = 0

    def liftover(self):
        for row in self.liftover_file:
            ensembl_id = row["gene_ensembl_id"]
            gene_name = row["gene_name"]
            if gene_name != "":
                if gene_name not in self.hg38_genes:
                    # Check if the hg19 gene name/ensembl id pair match an hg38 gene name.
                    # If so, use that hg38 gene name.
                    if (gene_name, ensembl_id) in self.hg19_genes:
                        gene_name = self.hg19_genes[(gene_name, ensembl_id)]
                    else:
                        self.unmapped_genes += 1
                        self._handle_unmapped(
                            f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['strand']}\t{gene_name}\t{row['gene_ensembl_id']}\n"
                        )
                        continue
                elif gene_name in self.hg38_genes:
                    ensembl_ids = self.hg38_genes[gene_name]
                    if ensembl_id not in ensembl_ids and len(ensembl_ids) == 1:
                        # In this case the gene has a new ensembl id
                        ensembl_id = list(ensembl_ids)[0]
                    elif len(ensembl_ids) > 1:
                        # In this case the gene has a new ensembl id, but in hg38 there
                        # are multiple ensembl ids to choose from, so we don't know which
                        # one is right and just don't map it

                        self.unmapped_genes += 1
                        self._handle_unmapped(
                            f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['strand']}\t{gene_name}\t{row['gene_ensembl_id']}\n"
                        )
                        continue

            strand = row["strand"]
            if strand == ".":
                strand = "+"

            new_coords = self.mapper.map_coordinates(row["chrom"], int(row["start"]), int(row["end"]), strand, "l")

            if new_coords is None:
                self._handle_unmapped(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                continue

            new_coords = new_coords[1::2]

            match len(new_coords):
                case 0:
                    self._handle_unmapped(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                case _:
                    try:
                        new_coords = merge_mapping(new_coords)
                    except:
                        self._handle_unmapped(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                        continue

                    # Chromosome names longer than 10 chars aren't allowed in the database.
                    # It's not technically unmapped, but may as well be
                    if len(new_coords[0]) > 10:
                        self._handle_unmapped(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{strand}\n")
                        continue

                    self._write_mapped(
                        {
                            "chrom": new_coords[0],
                            "start": new_coords[1],
                            "end": new_coords[2],
                            "strand": new_coords[3] if row["strand"] != "." else ".",
                            "gene_name": gene_name,
                            "gene_ensembl_id": ensembl_id,
                            "raw_p_val": row["raw_p_val"],
                            "adj_p_val": row["adj_p_val"],
                            "effect_size": row["effect_size"],
                            "facets": row["facets"],
                        }
                    )

        self.write_summary()


def get_args():
    parser = argparse.ArgumentParser("A tool for lifting over experiment and analysis data")
    parser.add_argument("-c", "--chain", required=True, help="The chain file containing the liftover data")
    parser.add_argument("-u", "--unmapped", type=argparse.FileType("w"), help="The file for logging unmapped regions")
    parser.add_argument("-s", "--summary", type=argparse.FileType("w"), help="A summary of the mapping results")
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
        regions = Experiment(args.experiment, args.output, args.unmapped, args.summary, mapper)
    elif args.analysis:
        genes_file = files("data_utilities").joinpath("hg38_genes.tsv")
        with open(genes_file, encoding="utf-8") as genes:
            regions = Analysis(args.analysis, args.output, args.unmapped, args.summary, mapper, genes=genes)

    regions.liftover()


def run_cli():
    args = get_args()
    run(args)


if __name__ == "__main__":
    run(get_args())
