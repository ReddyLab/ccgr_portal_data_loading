import argparse
import csv

from data_utilities.analysis_validation import AnalysisField
from data_utilities.experiment_validation import ExperimentField
from reddylab.common import *


def gene_ensembl_mapping(genes_file):
    gene_name_map = {}
    reader = csv.reader(genes_file, delimiter="\t")
    for row in reader:
        gene_name_map[row[0]] = row[1]

    return gene_name_map


def get_features(feature_tsv):
    reader = csv.DictReader(feature_tsv, quoting=csv.QUOTE_NONE)

    new_dhss = set()

    for line in reader:
        chrom_name = line["dhs_chrom"]
        dhs_start = int(line["dhs_start"])
        dhs_end = int(line["dhs_end"])

        dhs_name = f"{chrom_name}:[{dhs_start},{dhs_end})"

        if dhs_name not in new_dhss:
            new_dhss.add(dhs_name)
            yield {
                ExperimentField.CHROM: chrom_name,
                ExperimentField.START: dhs_start,
                ExperimentField.END: dhs_end,
                ExperimentField.STRAND: ".",
            }


def get_observations(results_file, gene_mapping_file):
    gene_name_map = gene_ensembl_mapping(gene_mapping_file)
    reader = csv.DictReader(results_file, quoting=csv.QUOTE_NONE)

    for line in reader:
        chrom_name = line["dhs_chrom"]
        dhs_start = int(line["dhs_start"])
        dhs_end = int(line["dhs_end"])

        significance = float(line["pval_empirical"])
        effect_size = float(line["avg_logFC"])
        if significance >= 0.01:
            direction = DirectionFacet.NON_SIGNIFICANT
        elif effect_size > 0:
            direction = DirectionFacet.ENRICHED
        elif effect_size < 0:
            direction = DirectionFacet.DEPLETED
        else:
            direction = DirectionFacet.NON_SIGNIFICANT

        yield {
            AnalysisField.CHROM: chrom_name,
            AnalysisField.START: dhs_start,
            AnalysisField.END: dhs_end,
            AnalysisField.STRAND: ".",
            AnalysisField.GENE_NAME: line["gene_symbol"],
            AnalysisField.GENE_ENSEMBL_ID: gene_name_map[line["gene_symbol"]],
            AnalysisField.RAW_P_VAL: float(line["p_val"]),
            AnalysisField.ADJ_P_VAL: significance,
            AnalysisField.EFFECT_SIZE: effect_size,
            AnalysisField.FACETS: f"{Facet.DIRECTION}={direction}",
        }


def get_args():
    parser = base_args()
    parser.add_argument(
        "-g",
        "--genes",
        required=True,
        type=argparse.FileType(encoding="utf-8"),
        help="gene name -> encode id mapping",
    )
    return parser.parse_args()


def run_cli():
    args = get_args()
    run(
        args,
        get_features(args.tested_elements),
        get_observations(args.observations, args.genes),
    )
