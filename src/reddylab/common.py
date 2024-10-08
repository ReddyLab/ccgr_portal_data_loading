import argparse
import csv
from enum import StrEnum
from pathlib import Path
from shutil import copyfile


from data_utilities import ANALYSIS_METADATA_FILE, EXPERIMENT_METADATA_FILE, TESTED_ELEMENTS_FILE, OBSERVATIONS_FILE
from data_utilities.analysis_validation import AnalysisField
from data_utilities.experiment_validation import ExperimentField


class Facet(StrEnum):
    GRNA_TYPE = "gRNA Type"
    PROMOTER = "Promoter Classification"
    DIRECTION = "Direction"


class GrnaFacet(StrEnum):
    POSITIVE_CONTROL = "Positive Control"
    NEGATIVE_CONTROL = "Negative Control"
    TARGETING = "Targeting"


class PromoterFacet(StrEnum):
    PROMOTER = "Promoter"
    NON_PROMOTER = "Non-promoter"


class DirectionFacet(StrEnum):
    DEPLETED = "Depleted Only"
    ENRICHED = "Enriched Only"
    NON_SIGNIFICANT = "Non-significant"
    BOTH = "Mixed"


class GenomeAssembly(StrEnum):
    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"
    HG19 = "hg19"
    HG38 = "hg38"


class ChromosomeStrands(StrEnum):
    POSITIVE = "+"
    NEGATIVE = "-"


def base_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("module")
    parser.add_argument("-t", "--tested-elements", required=True, type=argparse.FileType(encoding="utf-8"))
    parser.add_argument("-o", "--observations", required=True, type=argparse.FileType(encoding="utf-8"))
    parser.add_argument("-e", "--experiment-metadata", required=True)
    parser.add_argument("-a", "--analysis-metadata", required=True)
    parser.add_argument("--output-dir")
    return parser


def run(
    args,
    features,
    observations,
):
    output_directory = Path(args.output_dir)
    output_directory.mkdir(exist_ok=True)
    copyfile(args.experiment_metadata, output_directory / EXPERIMENT_METADATA_FILE)
    copyfile(args.analysis_metadata, output_directory / ANALYSIS_METADATA_FILE)

    with open(output_directory / Path(TESTED_ELEMENTS_FILE), "w", encoding="utf-8") as tested_elements_file:
        writer = csv.DictWriter(
            tested_elements_file, fieldnames=list(ExperimentField), delimiter="\t", quoting=csv.QUOTE_NONE
        )
        writer.writeheader()
        for feature in features:
            writer.writerow(feature)

    with open(output_directory / Path(OBSERVATIONS_FILE), "w", encoding="utf-8") as observations_files:
        writer = csv.DictWriter(
            observations_files, fieldnames=list(AnalysisField), delimiter="\t", quoting=csv.QUOTE_NONE
        )
        writer.writeheader()
        for observation in observations:
            writer.writerow(observation)
