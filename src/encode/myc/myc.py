import argparse
import asyncio
import csv
import gzip
import json
import os
import ssl
import tempfile
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Optional

import httpx
import truststore

from data_utilities.experiment_validation import FIELD_NAMES as EXPERIMENT_FIELD_NAMES
from data_utilities.analysis_validation import FIELD_NAMES as ANALYSIS_FIELD_NAMES
from encode import ScreenMetadata, TESTED_ELEMENTS_FILE, OBSERVATIONS_FILE, first


class MycMetadata(ScreenMetadata):
    @property
    def source_type(self) -> str:
        return "Called Regulatory Element"

    @cached_property
    def parent_source_type(self):
        return None

    @cached_property
    def tested_elements_description(self):
        ref_files = self.screen_info["files"]
        guides_file = [f for f in ref_files if "element quantifications" == f["output_type"]][0]

        return guides_file["aliases"][0]

    @cached_property
    def results_description(self) -> str:
        af = self._get_analysis_files()
        q_file = first(af, lambda x: x["output_category"] == "quantification")
        header = "" if q_file is None else f"{q_file['output_type'].capitalize()}"

        return header

    @cached_property
    def guide_quantification_url(self) -> str:
        return [
            file["cloud_metadata"]["url"]
            for file in self.screen_info["files"]
            if file["output_type"] == "element quantifications"
        ][0]


class BingRenPeakTilingMetadata(MycMetadata):
    @cached_property
    def assay(self):
        return "Proliferation screen"

    @property
    def p_val_threshold(self):
        return 0.2


class BingRenFullTilingMetadata(MycMetadata):
    @cached_property
    def assay(self):
        return "Proliferation screen"


class PardisSabetiMetadata(MycMetadata):
    pass


@dataclass
class EncodeData:
    chrom: str
    start: int
    end: int
    effect_size: float
    p_value: float
    adjusted_p_value: float
    significant: bool
    strand: Optional[str] = None
    target_gene_name: Optional[str] = None
    target_gene_ensembl_id: Optional[str] = None

    def element_data(self):
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": self.strand if self.strand is not None else ".",
        }

    def observation_data(self):
        result = {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": self.strand if self.strand is not None else ".",
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

        return result


def gen_metadata(screen_metadata_file, output_directory):
    screen = json.load(screen_metadata_file)

    match screen["lab"]["title"]:
        case "Bing Ren, UCSD":
            match screen["elements_references"][0]["crispr_screen_tiling"]:
                case "peak tiling":
                    metadata = BingRenPeakTilingMetadata(screen)
                case "full tiling":
                    metadata = BingRenFullTilingMetadata(screen)
                case tiling:
                    raise ValueError(f"Unknown tiling: {tiling}")
        case "Pardis Sabeti, Broad":
            metadata = PardisSabetiMetadata(screen)
        case lab:
            raise ValueError(f"Unknown lab: {lab}")

    expr_metadata = metadata.build_experiment(output_directory)
    analysis_metadata = metadata.build_analysis(output_directory)

    expr_path = output_directory / "experiment.json"
    analysis_path = output_directory / "analysis001.json"
    if not expr_path.exists():
        expr_path.touch()
    if not analysis_path.exists():
        analysis_path.touch()
    expr_path.write_text(json.dumps(expr_metadata))
    analysis_path.write_text(json.dumps(analysis_metadata))

    return metadata


def line_to_data(line: list[str]) -> EncodeData:
    (
        chrom,
        chrom_start,
        chrom_end,
        _,
        effect_size,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        gene_name,
        gene_ensembl_id,
        _,
        _,
        significant,
        p_value,
        p_value_adjusted,
        _,
        _,
        _,
        _,
        _,
    ) = line

    significant = True if significant == "TRUE" else False

    try:
        p_value = float(p_value)
    except ValueError:
        p_value = 0.0 if significant else 1.0

    try:
        p_value_adjusted = float(p_value_adjusted)
    except ValueError:
        p_value_adjusted = 0.0 if significant else 1.0

    if gene_name.startswith("Peak_"):
        gene_name = None
        gene_ensembl_id = None

    return EncodeData(
        chrom=chrom,
        start=int(chrom_start),
        end=int(chrom_end),
        effect_size=float(effect_size),
        target_gene_name=gene_name,
        target_gene_ensembl_id=gene_ensembl_id,
        p_value=p_value,
        adjusted_p_value=p_value_adjusted,
        significant=significant,
    )


def download_file(url, output_path):
    ssl_context = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    with httpx.Client(verify=ssl_context) as client:
        # Don't re-download files we've already downloaded
        if Path(output_path).exists():
            return
        response = client.get(url, timeout=5)
        with open(output_path, "wb") as output:
            output.write(response.content)


def read_data(metadata, output_directory):
    with tempfile.TemporaryDirectory(dir=output_directory) as temp_dir:
        guides_url: str = metadata.guide_quantification_url
        guides_file = temp_dir / Path(os.path.basename(guides_url))

        download_file(guides_url, guides_file)

        with (
            gzip.open(guides_file, mode="rt", encoding="utf8")
            if guides_url.endswith(".gz")
            else open(guides_file, encoding="utf8")
        ) as input_file:
            reader = csv.reader(input_file, delimiter="\t")
            line = next(reader)
            if line[0] != "chrom":
                yield line_to_data(line)

            for line in reader:
                yield line_to_data(line)


def gen_data(results_gen, output_directory):
    tested_elements_file = open(output_directory / Path(TESTED_ELEMENTS_FILE), "w", encoding="utf-8")
    tested_elements_csv = csv.DictWriter(
        tested_elements_file, fieldnames=EXPERIMENT_FIELD_NAMES, delimiter="\t", quoting=csv.QUOTE_NONE
    )
    tested_elements_csv.writeheader()

    observations_file = open(output_directory / Path(OBSERVATIONS_FILE), "w", encoding="utf-8")
    observations_csv = csv.DictWriter(
        observations_file, fieldnames=ANALYSIS_FIELD_NAMES, delimiter="\t", quoting=csv.QUOTE_NONE
    )
    observations_csv.writeheader()

    for result in results_gen:
        tested_elements_csv.writerow(result.element_data())
        observations_csv.writerow(result.observation_data())


def get_args():
    parser = argparse.ArgumentParser(
        description="Converts encode element quantification bed files to portal tested element and observation files"
    )
    parser.add_argument("--sm", "--screen-metdata", required=True, type=argparse.FileType(encoding="utf8"))
    parser.add_argument("-o", "--output-directory", required=True, type=Path)

    return parser.parse_args()


def run_cli():
    args = get_args()
    args.output_directory.mkdir(exist_ok=True)
    metadata = gen_metadata(args.sm, args.output_directory)
    metadata_paths = gen_data(read_data(metadata, args.output_directory), args.output_directory)
