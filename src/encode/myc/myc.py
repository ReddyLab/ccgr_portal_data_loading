import argparse
import csv
from typing import Optional
from dataclasses import dataclass


@dataclass
class EncodeData:
    chrom: str
    start: int
    end: int
    strand: Optional[str]
    target_gene_name: Optional[str]
    target_gene_ensembl_id: Optional[str]
    effect_size: float
    p_value: float
    adjusted_p_value: float
    significant: bool

    def element_data(self):
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": self.strand if self.strand is not None else ".",
            "bounds": "[)",
        }

    def observation_data(self):
        result = {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": self.strand if self.strand is not None else ".",
            "bounds": "[)",
            "raw_p_val": self.p_value,
            "adj_p_val": self.adjusted_p_value,
            "effect_size": self.effect_size,
        }

        if self.target_gene_name is not None:
            result["gene_name"] = self.target_gene_name
            result["gene_ensembl_id"] = self.target_gene_ensembl_id

        if not self.significant:
            result["facets"] = "Direction=Non-significant"
        else:
            if self.effect_size >= 0:
                result["facets"] = "Direction=Enriched Only"
            else:
                result["facets"] = "Direction=Depleted Only"


def read_data(input) -> list[EncodeData]:
    pass


def gen_metadata(screen_metadata, output_directory):
    pass


def gen_data(results, output_directory):
    pass


def get_args():
    parser = argparse.ArgumentParser(
        description="Converts encode element quantification bed files to portal tested element and observation files"
    )
    parser.add_argument("-i", "--input", required=True, type=argparse.FileType())
    parser.add_argument("--sm", "--screen-metdata", required=True, type=argparse.FileType())
    parser.add_argument("-o", "--output-directory", required=True)

    return parser.parse_args()


def run_cli():
    args = get_args()
    gen_metadata(args.sm, args.output_directory)
    results = read_data(args.input)
    metadata_paths = gen_data(results, args.output_directory)
