import csv

from data_utilities.analysis_validation import AnalysisField
from data_utilities.experiment_validation import ExperimentField
from reddylab.common import *


def get_features(feature_tsv):
    reader = csv.DictReader(feature_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)

    new_dhss = set()

    for line in reader:
        chrom_name = line["chrom"]
        dhs_start = int(line["chromStart"])
        dhs_end = int(line["chromEnd"])

        dhs_name = f"{chrom_name}:[{dhs_start},{dhs_end})"

        if dhs_name not in new_dhss:
            new_dhss.add(dhs_name)
            yield {
                ExperimentField.CHROM: chrom_name,
                ExperimentField.START: dhs_start,
                ExperimentField.END: dhs_end,
                ExperimentField.STRAND: ".",
            }


def get_observations(results_file):
    reader = csv.DictReader(results_file, delimiter="\t", quoting=csv.QUOTE_NONE)

    for line in reader:
        chrom_name = line["chrom"]
        dhs_start = int(line["chromStart"])
        dhs_end = int(line["chromEnd"])

        effect_size_field = line["wgCERES_score_top3_wg"].strip()
        if effect_size_field == "":
            effect_size = None
        else:
            effect_size = float(effect_size_field)

        match line["direction_wg"]:
            case "non_sig":
                direction = DirectionFacet.NON_SIGNIFICANT
            case "enriched":
                direction = DirectionFacet.ENRICHED
            case "depleted":
                direction = DirectionFacet.DEPLETED
            case "both":
                direction = DirectionFacet.BOTH
            case _:
                direction = None

        yield {
            AnalysisField.CHROM: chrom_name,
            AnalysisField.START: dhs_start,
            AnalysisField.END: dhs_end,
            AnalysisField.STRAND: ".",
            # line[pValue] is -log10(actual p-value), so raw_p_value uses the inverse operation
            AnalysisField.RAW_P_VAL: pow(10, -float(line["pValue"])),
            # line[pValue] is -log10(actual p-value), but we want significance between 0 and 1
            # we perform the inverse operation.
            AnalysisField.ADJ_P_VAL: pow(10, -float(line["pValue"])),
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
