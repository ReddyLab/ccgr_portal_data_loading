import csv

from data_utilities.analysis_validation import AnalysisField
from data_utilities.experiment_validation import ExperimentField
from reddylab.common import *


def get_features(feature_tsv):
    reader = csv.DictReader(feature_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    for line in reader:
        yield {
            ExperimentField.CHROM: line["seqnames"],
            ExperimentField.START: line["start"],
            ExperimentField.END: line["end"],
            ExperimentField.STRAND: ".",
        }


def get_observations(observation_file):
    reader = csv.DictReader(observation_file, delimiter="\t", quoting=csv.QUOTE_NONE)

    for line in reader:
        effect_size_field = line["logFC"].strip()
        if effect_size_field == "":
            effect_size = None
        else:
            effect_size = float(effect_size_field)

        log_significance = float(line["minusLog10PValue"])
        if log_significance < 2:  # p-value < 0.01, we're being stricter about significance with this data set
            direction = DirectionFacet.NON_SIGNIFICANT
        elif effect_size is not None and effect_size > 0:
            direction = DirectionFacet.ENRICHED
        elif effect_size is not None and effect_size < 0:
            direction = DirectionFacet.DEPLETED
        else:
            direction = None

        significance = pow(10, -float(log_significance))
        p_val = significance

        yield {
            AnalysisField.CHROM: line["seqnames"],
            AnalysisField.START: line["start"],
            AnalysisField.END: line["end"],
            AnalysisField.STRAND: ".",
            AnalysisField.RAW_P_VAL: significance,
            AnalysisField.ADJ_P_VAL: significance,
            AnalysisField.EFFECT_SIZE: effect_size_field,
            AnalysisField.FACETS: f"{Facet.DIRECTION}={direction}",
        }


def run_cli():
    args = base_args().parse_args()
    run(
        args,
        get_features(args.tested_elements),
        get_observations(args.observations),
    )
