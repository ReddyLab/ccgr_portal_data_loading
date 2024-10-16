import csv

from data_utilities.analysis_validation import AnalysisField
from data_utilities.experiment_validation import ExperimentField
from reddylab.common import *


def get_features(feature_file):
    reader = csv.DictReader(feature_file, delimiter="\t", quoting=csv.QUOTE_NONE)
    guide_set = set()
    for line in reader:
        grna_id = line["grna"]
        if grna_id in guide_set:
            continue

        guide_set.add(grna_id)

        grna_info = grna_id.split("-")

        # Skip non-targeting guides and guides with no assigned enhancer
        if not grna_info[0].startswith("chr") or line["dhs.chr"] == "NA":
            continue

        categorical_facets = []
        grna_type = line["type"]
        grna_promoter_class = line["annotation_manual"]

        if grna_type == "targeting":
            categorical_facets.append((Facet.GRNA_TYPE, GrnaFacet.TARGETING))
        elif grna_type.startswith("positive_control") == "":
            categorical_facets.append((Facet.GRNA_TYPE, GrnaFacet.POSITIVE_CONTROL))

        if grna_promoter_class == "promoter":
            categorical_facets.append((Facet.PROMOTER, PromoterFacet.PROMOTER))
        else:
            categorical_facets.append((Facet.PROMOTER, PromoterFacet.NON_PROMOTER))

        yield {
            ExperimentField.CHROM: line["grna.chr"],
            ExperimentField.START: line["grna.start"],
            ExperimentField.END: line["grna.end"],
            ExperimentField.STRAND: line["grna.strand"],
            ExperimentField.PARENT_CHROM: line["dhs.chr"],
            ExperimentField.PARENT_START: line["dhs.start"],
            ExperimentField.PARENT_END: line["dhs.end"],
            ExperimentField.PARENT_STRAND: ".",
            ExperimentField.FACETS: ";".join(f"{k}={v}" for k, v in categorical_facets),
        }


def get_observations(observation_file):
    reader = csv.DictReader(observation_file, delimiter="\t", quoting=csv.QUOTE_NONE)

    for i, line in enumerate(reader):
        # every other line in this file is basically a duplicate of the previous line
        if i % 2 == 0:
            continue

        significance = float(line["pval_fdr_corrected"])
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
            AnalysisField.CHROM: line["grna.chr"],
            AnalysisField.START: line["grna.start"],
            AnalysisField.END: line["grna.end"],
            AnalysisField.STRAND: line["grna.strand"],
            AnalysisField.GENE_NAME: line["gene_symbol"],
            AnalysisField.GENE_ENSEMBL_ID: line["gene_stable_id"],
            AnalysisField.RAW_P_VAL: line["p_val"],
            AnalysisField.ADJ_P_VAL: significance,
            AnalysisField.EFFECT_SIZE: effect_size,
            AnalysisField.FACETS: f"{Facet.DIRECTION}={direction}",
        }


def run_cli():
    args = base_args().parse_args()
    run(
        args,
        get_features(args.tested_elements),
        get_observations(args.observations),
    )
