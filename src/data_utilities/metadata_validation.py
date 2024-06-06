import argparse
import json
from importlib.resources import files

from jsonschema import validate


def get_args():
    parser = argparse.ArgumentParser(description="Validate experiment and analysis metadata")
    parser.add_argument("metadata_file", help="The metadata file to validate", type=argparse.FileType())

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-e", "--experiment", action="store_true", help="Validate experiment metadata")
    group.add_argument("-a", "--analysis", action="store_true", help="Validate analysis metadata")

    return parser.parse_args()


def validate_metadata(args):
    if args.experiment:
        validation_schema_file = files("data_utilities").joinpath("experiment_metadata.schema.json")
    elif args.analysis:
        validation_schema_file = files("data_utilities").joinpath("analysis_metadata.schema.json")

    with open(validation_schema_file, encoding="utf-8") as validation_schema:
        schema = json.load(validation_schema)

    metadata = json.load(args.metadata_file)

    validate(instance=metadata, schema=schema)


def run_cli():
    args = get_args()
    validate_metadata(args)
