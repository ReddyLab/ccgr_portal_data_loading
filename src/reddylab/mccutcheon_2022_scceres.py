import argparse
import csv

from data_utilities.analysis_validation import AnalysisField
from data_utilities.experiment_validation import ExperimentField
from reddylab.common import *
from reddylab.mccutcheon_consts import DROP_GENE_NAMES, TRIM_GENE_NAMES


def gene_ensembl_mapping(genes_file):
    gene_name_map = {}
    reader = csv.reader(genes_file, delimiter="\t")
    for row in reader:
        if row[2] != "Gene Expression":
            continue
        gene_name_map[row[1]] = row[0]

    return gene_name_map


def read_dhss(dhs_file):
    dhss = {}
    reader = csv.DictReader(dhs_file, delimiter="\t", quoting=csv.QUOTE_NONE)
    for line in reader:
        dhss[line["grna_id"]] = (
            line["dhs_chrom"],
            int(line["dhs_start"]),
            int(line["dhs_end"]),
        )
    return dhss


def get_features(experiment_file, features_file):
    reader = csv.DictReader(experiment_file, delimiter="\t", quoting=csv.QUOTE_NONE)
    dhs_info = read_dhss(features_file)
    guide_set = set()
    for line in reader:
        grna_id = line["grna"]
        if grna_id in guide_set:
            continue

        guide_set.add(grna_id)
        dhs_chrom_name, dhs_start, dhs_end = dhs_info[line["grna"]]

        start = int(line["start"])
        end = int(line["end"])
        strand = line["Strand"]

        # Convert from "(]" to "[)" bounds
        if strand == "-":
            start += 1
            end += 1

        yield {
            ExperimentField.CHROM: line["chr"],
            ExperimentField.START: line["start"],
            ExperimentField.END: line["end"],
            ExperimentField.STRAND: line["Strand"],
            ExperimentField.PARENT_CHROM: dhs_chrom_name,
            ExperimentField.PARENT_START: dhs_start,
            ExperimentField.PARENT_END: dhs_end,
            ExperimentField.PARENT_STRAND: ".",
            ExperimentField.FACETS: f"{Facet.GRNA_TYPE}={GrnaFacet.TARGETING}",
            ExperimentField.MISC: f"grna={grna_id}",
        }


def get_observations(observation_file, genes_file):
    reader = csv.DictReader(observation_file, delimiter="\t", quoting=csv.QUOTE_NONE)
    gene_name_map = gene_ensembl_mapping(genes_file)

    for line in reader:
        target_gene = line["target_gene"]

        if target_gene in TRIM_GENE_NAMES:
            target_gene = target_gene[:-2]

        if gene_name_map[target_gene] in DROP_GENE_NAMES:
            continue

        strand = line["Strand"]
        chrom_name = line["chr"]
        grna_start = int(line["start"])
        grna_end = int(line["end"])

        # Convert from "(]" to "[)" bounds
        if strand == "-":
            grna_start += 1
            grna_end += 1

        significance = float(line["p_val_adj"])
        effect_size = float(line["avg_log2FC"])
        if significance >= 0.05:
            direction = DirectionFacet.NON_SIGNIFICANT
        elif effect_size > 0:
            direction = DirectionFacet.ENRICHED
        elif effect_size < 0:
            direction = DirectionFacet.DEPLETED
        else:
            direction = DirectionFacet.NON_SIGNIFICANT

        yield {
            AnalysisField.CHROM: chrom_name,
            AnalysisField.START: grna_start,
            AnalysisField.END: grna_end,
            AnalysisField.STRAND: strand,
            AnalysisField.GENE_NAME: target_gene,
            AnalysisField.GENE_ENSEMBL_ID: gene_name_map[target_gene],
            AnalysisField.RAW_P_VAL: line["p_val"],
            AnalysisField.ADJ_P_VAL: line["p_val_adj"],
            AnalysisField.EFFECT_SIZE: line["avg_log2FC"],
            AnalysisField.FACETS: f"{Facet.DIRECTION}={direction}",
        }


def get_args():
    parser = base_args()
    parser.add_argument("-f", "--features", required=True, type=argparse.FileType(encoding="utf-8"))
    parser.add_argument("-g", "--genes", required=True, type=argparse.FileType(encoding="utf-8"))
    return parser.parse_args()


def run_cli():
    args = get_args()
    run(
        args,
        get_features(args.tested_elements, args.features),
        get_observations(args.observations, args.genes),
    )
