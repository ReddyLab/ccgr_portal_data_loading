#!/usr/bin/env python

import argparse
import csv
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
    "parent_chrom",
    "parent_start",
    "parent_end",
    "parent_strand",
    "parent_bounds",
    "facets",
    "misc",
]

FIELD_NAMES_SET = set(FIELD_NAMES)


def validate_misc(key_values: str):
    try:
        {k: v for k, v in [pair.split("=") for pair in key_values.split(";")]}
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

        if not (
            line["parent_chrom"]
            == line["parent_start"]
            == line["parent_end"]
            == line["parent_strand"]
            == line["parent_bounds"]
            == ""
        ):
            if not validate_chrom(line["parent_chrom"]):
                errors.add(f"Invalid chromosome: {line['chrom']}")

            if not validate_location(line["parent_start"], line["parent_end"]):
                errors.add(f"Invalid parent location: {line['parent_start']}, {line['parent_end']} (line {i})")

            if not validate_strand(line["parent_strand"]):
                errors.add(f"Invalid parent strand: {line['parent_strand']} (line {i})")

            if not validate_bounds(line["bounds"]):
                errors.add(f"Invalid parent bounds: {line['parent_bounds']} (line {i})")

        if line["facets"] != "" and not validate_facets(line["facets"]):
            errors.add(f"Invalid facets: {line['facets']} (line {i})")

        if line["misc"] != "" and not validate_misc(line["misc"]):
            errors.add(f"Invalid misc: {line['misc']} (line {i})")

    return errors.error_list


def validate(experiment_data):
    errors = []

    reader = csv.DictReader(experiment_data, delimiter="\t")

    errors += validate_fieldnames(reader, FIELD_NAMES_SET)
    errors += validate_data(reader)

    return errors


def get_args():
    parser = argparse.ArgumentParser(description="Validator for CCGR Portal experiment data")
    parser.add_argument("experiment data", help="Experiment data file to validate")
    return parser.parse_args()


def run_cli():
    args = get_args()
    experiment_data_file = getattr(args, "experiment data")
    with open(experiment_data_file, encoding="utf-8") as experiment_data:
        errors = validate(experiment_data)

    if len(errors) == 0:
        print(f"Experiment data is valid: {experiment_data_file}")
    else:
        error_list = "\n".join(errors)
        print(f"The following issues were found:\n{error_list}")
        sys.exit(1)
