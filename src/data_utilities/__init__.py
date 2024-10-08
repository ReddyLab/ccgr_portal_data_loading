import re

TESTED_ELEMENTS_FILE = "tested_elements.tsv"
OBSERVATIONS_FILE = "observations.tsv"
VALID_STRANDS = {".", "+", "-"}


class ErrorSet:
    def __init__(self):
        self.error_list = []
        self.error_set = set()

    def add(self, error):
        if error not in self.error_set:
            self.error_set.add(error)
            self.error_list.append(error)


def validate_fieldnames(reader, field_names):
    errors = []

    file_fieldnames_set = set(reader.fieldnames)

    diff = field_names - file_fieldnames_set

    if len(diff) > 0:
        errors.append(f"Data file is missing field(s): {diff}")

    diff = file_fieldnames_set - field_names
    if len(diff) > 0:
        errors.append(f"Data file includes extra field(s): {diff}")

    return errors


def validate_chrom(chrom: str):
    return re.match(r"chr([0-9]|1[0-9]|2[0-2]|[XYM])", chrom) is not None


def validate_location(start: str, end: str):
    if not (start.isdigit() and end.isdigit()):
        return False

    return 0 <= int(start) < int(end)


def validate_strand(strand: str):
    return strand in VALID_STRANDS


def validate_facets(key_values: str):
    return all(len(pair.split("=")) == 2 for pair in key_values.split(";"))
