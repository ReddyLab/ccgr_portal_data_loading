#!/usr/bin/env python

import argparse
import csv
import re
import sys

from . import (
    ErrorSet,
    validate_bounds,
    validate_chrom,
    validate_facets,
    validate_fieldnames,
    validate_location,
    validate_strand,
)

FIELD_NAMES = [
    "chrom",
    "start",
    "end",
    "strand",
    "bounds",
    "gene_name",
    "gene_ensembl_id",
    "raw_p_val",
    "adj_p_val",
    "effect_size",
    "facets",
]

FIELD_NAMES_SET = set(FIELD_NAMES)


def validate_gene_name(gene_name: str):
    return isinstance(gene_name, str)


def validate_ensembl_id(ensembl_id: str):
    return re.match(r"ENSG\d{11}", ensembl_id) is not None


def validate_float(float_num: str):
    try:
        float(float_num)
    except ValueError:
        return False

    return True


def validate_data(reader):
    errors = ErrorSet()

    for i, line in enumerate(reader, start=1):
        if not validate_chrom(line["chrom"]):
            errors.add(f"Invalid chromosome: {line['chrom']}")

        if not validate_location(line["start"], line["end"]):
            errors.add(f"Invalid location: {line['start']}, {line['end']} (line {i})")

        if not validate_strand(line["strand"]):
            errors.add(f"Invalid strand: {line['strand']} (line {i})")

        if not validate_bounds(line["bounds"]):
            errors.add(f"Invalid bounds: {line['bounds']} (line {i})")

        if not validate_gene_name(line["gene_name"]):
            errors.add(f"Invalid gene name: {line['gene_name']} (line {i})")

        if not validate_ensembl_id(line["gene_ensembl_id"]):
            errors.add(f"Invalid gene ensembl id: {line['gene_ensembl_id']} (line {i})")

        if not validate_float(line["raw_p_val"]):
            errors.add(f"Invalid raw p value: {line['raw_p_val']} (line {i})")

        if not validate_float(line["adj_p_val"]):
            errors.add(f"Invalid adjusted p value: {line['adj_p_val']} (line {i})")

        if not validate_float(line["effect_size"]):
            errors.add(f"Invalid effect size: {line['effect_size']} (line {i})")

        if line["facets"] != "" and not validate_facets(line["facets"]):
            errors.add(f"Invalid facets: {line['facets']} (line {i})")

    return errors.error_list


def validate(experiment_data):
    errors = []

    reader = csv.DictReader(experiment_data, delimiter="\t")

    errors += validate_fieldnames(reader, FIELD_NAMES_SET)
    errors += validate_data(reader)

    return errors


def get_args():
    parser = argparse.ArgumentParser(description="Validator for CCGR Portal analysis data")
    parser.add_argument("analysis data", help="Analysis data file to validate")
    return parser.parse_args()


def run_cli():
    args = get_args()
    analysis_data_file = getattr(args, "analysis data")
    with open(analysis_data_file, encoding="utf-8") as experiment_data:
        errors = validate(experiment_data)

    if len(errors) == 0:
        print(f"File is valid analysis data: {analysis_data_file}")
    else:
        error_list = "\n".join(errors)
        print(f"The following issues were found:\n{error_list}")
        sys.exit(1)
