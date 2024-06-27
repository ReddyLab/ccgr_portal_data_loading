import argparse
import asyncio
import csv
import json
import os.path
import re
import tempfile
from datetime import datetime
from collections import defaultdict
from os import SEEK_SET
from pathlib import Path
from typing import Optional

import httpx

from data_utilities.experiment_validation import FIELD_NAMES as EXPERIMENT_FIELD_NAMES
from data_utilities.analysis_validation import FIELD_NAMES as ANALYSIS_FIELD_NAMES

# These experiments have bad data -- the files are the wrong format
# and don't contain all the necessary information for
BAD_ENGREITZ_EXPERIMENTS = {
    "ENCSR108MPF",
    "ENCSR144TTP",
    "ENCSR194SMN",
    "ENCSR273LJW",
    "ENCSR376SMC",
    "ENCSR380BWA",
    "ENCSR416EUJ",
    "ENCSR473FYG",
    "ENCSR525DAZ",
    "ENCSR810OLG",
}
GOOD_ENGREITZ_EXPERIMENTS = {
    "ENCSR006WCB",
    "ENCSR019BVO",
    "ENCSR020TLT",
    "ENCSR029BRD",
    "ENCSR101JKC",
    "ENCSR114BIZ",
    "ENCSR116BPW",
    "ENCSR138CXY",
    "ENCSR154GRV",
    "ENCSR176XPW",
    "ENCSR245NIU",
    "ENCSR245ZZP",
    "ENCSR306ZVG",
    "ENCSR324QDY",
    "ENCSR336WSZ",
    "ENCSR532VIG",
    "ENCSR541HXR",
    "ENCSR617XYY",
    "ENCSR661PUY",
    "ENCSR752CCF",
    "ENCSR760TSA",
    "ENCSR858PBD",
    "ENCSR876DIL",
    "ENCSR884HBT",
    "ENCSR884XXE",
    "ENCSR898QMI",
    "ENCSR905FEH",
    "ENCSR905MSH",
    "ENCSR922PHL",
    "ENCSR954IYH",
}
GRCH37 = "GRCh37"
GRCH37_ASSEMBLIES = {GRCH37, "hg19"}
TISSUE_TYPES = {"K562": "Bone Marrow"}
P_VAL_THRESHOLD = 0.05
TESTED_ELEMENTS_FILE = "tested_elements.tsv"
OBSERVATIONS_FILE = "observations.tsv"


def first(iterable, test):
    for x in iterable:
        if test(x):
            return x


def gen_experiment_data(guides_file, dhs_file, results_file, strand_file):
    # The experiment data requires 4 files from the ENCODE data set:
    # 1) element quantifications aka results_file
    # 2) elements reference (guides) aka element_file
    # 3) elements reference (DHS peaks) aka parent_element_file
    # 4) a guide quantifications file aka guide_quant_file
    #
    # Loading the files requires compiling data from all four files. The element quantifications (1) file
    # Tells us which peaks DHS peaks from the elements reference (3) to include. The elements reference (3)
    # includes all the oligo ids for a given DHS peak which we can use to get the guids from elements references (2).
    # Unfortunately, elements reference (2) doesn't have the strand information! For this we need one of the guide
    # quantification files. We can match the guide in (2) to the guide information in (4) via the guide sequence.
    #
    # Once we have all the guide and DHS peak information we can add the guides and dhs peaks they are children of to the DB

    element_tsv = open(guides_file, encoding="utf-8")
    element_reader = csv.DictReader(element_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)

    parent_element_tsv = open(dhs_file, encoding="utf-8")
    parent_reader = csv.DictReader(parent_element_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)

    oligo_to_parents = {}

    results_tsv = open(results_file, encoding="utf-8")
    results_reader = csv.DictReader(results_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)

    #
    # Figure out which DHS peaks to include for this experiment
    #
    result_targets = {
        f'{line["chrPerturbationTarget"]}:{line["startPerturbationTarget"]}-{line["endPerturbationTarget"]}'
        for line in results_reader
    }

    #
    # Read guide strand and type information from the guide quantification file.
    # The type information is used for the GrnaType facet
    #
    strand_tsv = open(strand_file, encoding="utf-8")
    strand_reader = csv.reader(strand_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    chrom_strands = ["+", "-"]
    guide_strands = {}
    guide_types = {}
    for line in strand_reader:
        if line[5] in chrom_strands:
            guide_strands[line[14]] = line[5]
        else:
            guide_strands[line[14]] = None
        guide_types[line[14]] = line[15]

    #
    # Build parent (DHS Peak) features and create the oligo->peak mapping
    #
    for line in parent_reader:
        oligo_id, parent_string = line["OligoID"], line["target"]
        parent_string.strip()
        if parent_string not in result_targets:
            continue

        if (parent_info := re.match(r"(.+):(\d+)-(\d+)", parent_string)) is not None:
            parent_chrom, parent_start, parent_end = (
                parent_info.group(1),
                int(parent_info.group(2)),
                int(parent_info.group(3)),
            )
        else:
            continue

        oligo_to_parents[oligo_id] = (parent_chrom, parent_start, parent_end)

    #
    # Build the guide features
    #
    for line in element_reader:
        oligo_id = line["OligoID"]
        if oligo_id not in oligo_to_parents:
            continue

        guide_seq = line["GuideSequence"]

        #
        # We previously filtered out guides invalid strands
        # Here, we skip over those guides
        #
        if guide_seq not in guide_strands or guide_strands[guide_seq] is None:
            continue

        parent_chrom, parent_start, parent_end = oligo_to_parents[oligo_id]

        element_chrom, element_start, element_end = (line["chr"], int(line["start"]), int(line["end"]))
        strand = guide_strands[guide_seq]

        if guide_types[guide_seq] == "targeting":
            guide_type = "Targeting"
        elif guide_types[guide_seq] == "negative_control":
            guide_type = "Negative Control"
        else:
            raise ValueError(f"Guide Type: {guide_types[guide_seq]}")

        yield (
            element_chrom,
            element_start,
            element_end,
            strand,
            "[)",
            parent_chrom,
            parent_start,
            parent_end,
            ".",
            "[)",
            f"gRNA Type={guide_type}",
            f"grna={guide_seq}",
        )

    element_tsv.close()
    parent_element_tsv.close()
    results_tsv.close()
    strand_tsv.close()


def gen_analysis_data(guides_file, dhs_file, results_file, strand_file):
    # The analysis data requires 4 files from the ENCODE data set:
    # 1) element quantifications aka results_file
    # 2) elements reference (guides) aka guides_file
    # 3) elements reference (DHS peaks) aka elements_file
    # 4) a guide quantifications file aka strand_file
    #
    # Loading the files requires compiling data from all four files. The element quantifications (1) file
    # Tells us which peaks DHS peaks from the elements reference (3) to include. The elements reference (3)
    # includes all the oligo ids for a given DHS peak which we can use to get the guids from elements references (2).
    # Unfortunately, elements reference (2) doesn't have the strand information! For this we need one of the guide
    # quantification files. We can match the guide in (2) to the guide information in (4) via the guide sequence.
    #
    # Once we have all the guide information we can match the guide (source) to the observation and target information
    # which are in the element quantifications (1) file.

    #
    # Much like when loading the experiment we have to use the results file to figure out which DHS peaks to include
    #
    results_tsv = open(results_file, encoding="utf-8")
    results_reader = csv.DictReader(results_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    result_targets = {
        f'{line["chrPerturbationTarget"]}:{line["startPerturbationTarget"]}-{line["endPerturbationTarget"]}'
        for line in results_reader
    }

    #
    # Match DHS Peaks to oligo ids and create a set of all oligo ids
    #
    parent_element_tsv = open(dhs_file, encoding="utf-8")
    parent_reader = csv.DictReader(parent_element_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    parent_elements = defaultdict(set)
    parent_oligos = set()
    for line in parent_reader:
        if line["target"] in result_targets:
            parent_elements[line["target"]].add(line["OligoID"])
            parent_oligos.add(line["OligoID"])

    #
    # Get all guides associated with DHS peaks
    #
    element_tsv = open(guides_file, encoding="utf-8")
    element_reader = csv.DictReader(element_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    elements = {}
    for e in element_reader:
        if e["OligoID"] not in parent_oligos:
            continue

        elements[e["OligoID"]] = (e["chr"], int(e["start"]), int(e["end"]), e["GuideSequence"])

    #
    # Get strands for guides
    #
    strand_tsv = open(strand_file, encoding="utf-8")
    strand_reader = csv.reader(strand_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    chrom_strands = ["+", "-"]
    guide_strands = {}
    for line in strand_reader:
        if line[5] in chrom_strands:
            guide_strands[line[14]] = line[5]
        else:
            guide_strands[line[14]] = None

    #
    # Go back through the results, matching guides to observations using the parent_elements dictionary
    #
    results_tsv.seek(0, SEEK_SET)
    results_reader = csv.DictReader(results_tsv, delimiter="\t", quoting=csv.QUOTE_NONE)
    for line in results_reader:
        chrom_name, start, end = (
            line["chrPerturbationTarget"],
            int(line["startPerturbationTarget"]),
            int(line["endPerturbationTarget"]),
        )
        parent_element = f"{chrom_name}:{start}-{end}"
        oligos = parent_elements[parent_element]

        raw_p_value = line["pValue"]
        adjusted_p_value = line["pValueAdjusted"]

        # An explanation of the effect size values, from an email with Ben Doughty:
        #
        # For the effect size calculations, since we do a 6-bin sort, we don't actually compute a log2-fold change.
        # Instead, we use the data to compute an effect size in "gene expression" space, which we normalize to the
        # negative controls. The values are then scaled, so what an effect size of -0.2 means is that this guide
        # decreased the expression of the target gene by 20%. An effect size of 0 would be no change in expression,
        # and an effect size of +0.1 would mean a 10% increase in expression. The lowest we can go is -1 (which means
        # total elimination of signal), and technically the effect size is unbounded in the positive direction,
        # although we never see _super_ strong positive guides with CRISPRi.
        effect_size = line["EffectSize"]

        target_gene_name = line["measuredGeneSymbol"]
        target_ensembl_id = line["measuredEnsemblID"]

        if float(adjusted_p_value) <= P_VAL_THRESHOLD:
            if float(effect_size) > 0:
                cat_facet = "Direction=Enriched Only"
            else:
                cat_facet = "Direction=Depleted Only"
        else:
            cat_facet = "Direction=Non-significant"

        for oligo in oligos:
            guide_chrom, guide_start, guide_end, guide_seq = elements[oligo]

            #
            # We previously filtered out guides invalid strands
            # Here, we skip over those guides
            #
            if guide_seq not in guide_strands or guide_strands[guide_seq] is None:
                continue

            # We previously filtered out guides invalid strands
            # Here, we skip over those guides
            if guide_seq in guide_strands:
                yield (
                    guide_chrom,
                    guide_start,
                    guide_end,
                    guide_strands[guide_seq],
                    "[)",
                    target_gene_name,
                    target_ensembl_id,
                    raw_p_value,
                    adjusted_p_value,
                    effect_size,
                    cat_facet,
                )

    results_tsv.close()
    element_tsv.close()
    parent_element_tsv.close()
    strand_tsv.close()


def normalize_assembly(assembly):
    if assembly in GRCH37_ASSEMBLIES:
        return GRCH37

    raise ValueError(f"Invalid assembly {assembly}")


def get_source_type(element_reference):
    source_types = element_reference.get("elements_selection_method")
    if source_types is None:
        return "Tested Element"

    st_set = set(source_types)
    if "accessible genome regions" in st_set:
        return "Chromatin Accessible Region"

    if "candidate cis-regulatory elements" in st_set:
        return "cCRE"

    if "DNase hypersensitive sites" in st_set:
        return "DHS"

    return "Tested Element"


def get_assembly(screen_info):
    if len(screen_info["assembly"]) > 0:
        assembly = screen_info["assembly"][0]
    else:
        assembly = screen_info["elements_references"][0]["assembly"][0]

    return assembly


def get_analysis_files(screen_info):
    af_accessions = screen_info["analyses"][0]["files"]
    return [file for file in screen_info["files"] if file["@id"] in af_accessions]


def build_experiment(screen_info, output_dir: Path):
    gene_assembly = normalize_assembly(get_assembly(screen_info))

    cell_line = screen_info["biosample_ontology"][0]["term_name"]
    tissue_type = TISSUE_TYPES[cell_line]
    crispr = screen_info["related_datasets"][0]["perturbation_type"]

    create_date = datetime.fromisoformat(screen_info["date_created"])

    ref_files = screen_info["elements_references"][0]["files"]
    guides_file = [f for f in ref_files if "guide" in f["aliases"][0]][0]

    # CRISPRi Flow-FISH screen of multiple loci in K562 with PrimeFlow readout of PQBP1
    experiment = {
        "name": f"{crispr} {screen_info['assay_title'][0]} in {cell_line}",
        "description": screen_info.get("biosample_summary"),
        "biosamples": [{"cell_type": cell_line, "tissue_type": tissue_type}],
        "assay": screen_info["assay_title"][0],
        "source type": "gRNA",
        "parent source type": get_source_type(screen_info["elements_references"][0]),
        "year": str(create_date.year),
        "lab": screen_info["lab"]["title"],
        "tested_elements_file": {
            "description": guides_file["aliases"][0],
            "filename": TESTED_ELEMENTS_FILE,
            "file_location": str((output_dir / Path(TESTED_ELEMENTS_FILE)).absolute()),
            "genome_assembly": gene_assembly,
        },
    }

    return experiment


def build_analysis(screen_info, output_dir: Path):
    gene_assembly = normalize_assembly(get_assembly(screen_info))

    af = get_analysis_files(screen_info)
    q_file = first(af, lambda x: x["output_category"] == "quantification")
    crispr = screen_info["related_datasets"][0]["perturbation_type"]
    cell_line = screen_info["biosample_ontology"][0]["term_name"]

    analysis = {
        "name": f"{crispr} {screen_info['assay_title'][0]} in {cell_line}",
        "description": screen_info.get("biosample_summary"),
        "source type": "gRNA",
        "genome_assembly": gene_assembly,
        "p_val_adj_method": "Benjamini-Hochberg",
        "p_val_threshold": 0.05,
        "results": {
            "description": f"{q_file['output_type'].capitalize()}\nFrom Ben Doughty, via Email:\nFor the effect size calculations, since we do a 6-bin sort, we don't actually compute a log2-fold change. Instead, we use the data to compute an effect size in \"gene expression\" space, which we normalize to the negative controls. The values are then scaled, so what an effect size of -0.2 means is that this guide decreased the expression of the target gene by 20%. An effect size of 0 would be no change in expression, and an effect size of +0.1 would mean a 10% increase in expression. The lowest we can go is -1 (which means total elimination of signal), and technically the effect size is unbounded in the positive direction, although we never see _super_ strong positive guides with CRISPRi. ",
            "filename": OBSERVATIONS_FILE,
            "file_location": str((output_dir / Path(OBSERVATIONS_FILE)).absolute()),
        },
    }

    if (desc := screen_info.get("description")) is not None:
        analysis["description"] = desc

    return analysis


# Guides
def get_guides_file_url(screen) -> str:
    ref_files = screen["elements_references"][0]["files"]
    return [f for f in ref_files if "guide" in f["aliases"][0]][0]["cloud_metadata"]["url"]


# This file helps map DHSs from the results file to guides
def get_elements_file_url(screen) -> str:
    ref_files = screen["elements_references"][0]["files"]
    return [f for f in ref_files if "element" in f["aliases"][0]][0]["cloud_metadata"]["url"]


# This is the results file
def get_element_quantification_url(screen) -> list[str]:
    return [
        file["cloud_metadata"]["url"] for file in screen["files"] if file["output_type"] == "element quantifications"
    ][0]


# We only need one because the only reason we use this is to get the
# guide strand
def get_guide_quantification_url(screen) -> list[str]:
    return [
        file["cloud_metadata"]["url"]
        for file in screen["related_datasets"][0]["files"]
        if file["output_type"] == "guide quantifications"
    ][0]


async def download_file(client, url, output_path):
    # Don't re-download files we've already downloaded
    if Path(output_path).exists():
        return
    response = await client.get(url, timeout=5)
    with open(output_path, "w", encoding="utf-8") as output:
        output.write(response.content.decode())


async def gen_data(metadata_path, output_path) -> tuple[Optional[dict], Optional[dict]]:
    curr_dir = Path(metadata_path)
    json_files = list(curr_dir.glob("*.json"))

    output_dir = Path(output_path)
    output_dir.mkdir(exist_ok=True)

    metadata_paths = []

    for file in json_files:
        with open(file, encoding="utf-8") as f:
            screen = json.load(f)

        if screen["accession"] not in GOOD_ENGREITZ_EXPERIMENTS:
            continue

        expr_dir = output_dir / Path(screen["accession"])
        expr_dir.mkdir(exist_ok=True)

        tested_elements_file = open(expr_dir / Path(TESTED_ELEMENTS_FILE), "w", encoding="utf-8")
        tested_elements_csv = csv.writer(tested_elements_file, delimiter="\t", quoting=csv.QUOTE_NONE)
        tested_elements_csv.writerow(EXPERIMENT_FIELD_NAMES)

        observations_file = open(expr_dir / Path(OBSERVATIONS_FILE), "w", encoding="utf-8")
        observations_csv = csv.writer(observations_file, delimiter="\t", quoting=csv.QUOTE_NONE)
        observations_csv.writerow(ANALYSIS_FIELD_NAMES)

        with tempfile.TemporaryDirectory(dir=output_path) as temp_dir:
            downloads = []
            guides_url = get_guides_file_url(screen)
            guides_file = temp_dir / Path(os.path.basename(guides_url))
            downloads.append((guides_url, guides_file))

            elements_url = get_elements_file_url(screen)
            elements_file = temp_dir / Path(os.path.basename(elements_url))
            downloads.append((elements_url, elements_file))

            results_url = get_element_quantification_url(screen)
            results_file = temp_dir / Path(os.path.basename(results_url))
            downloads.append((results_url, results_file))

            guide_strand_url = get_guide_quantification_url(screen)
            guide_strand_file = temp_dir / Path(os.path.basename(guide_strand_url))
            downloads.append((guide_strand_url, guide_strand_file))

            async with httpx.AsyncClient() as client:
                tasks = [asyncio.create_task(download_file(client, url, file)) for url, file in downloads]
                await asyncio.wait(tasks)

            print(f"Working on {screen['accession']}")
            tested_elements_csv.writerows(
                gen_experiment_data(guides_file, elements_file, results_file, guide_strand_file)
            )
            observations_csv.writerows(gen_analysis_data(guides_file, elements_file, results_file, guide_strand_file))

        tested_elements_file.close()
        observations_file.close()

        expr_metadata = build_experiment(screen, expr_dir)
        analysis_metadata = build_analysis(screen, expr_dir)

        metadata_paths.append(write_experiment(expr_metadata, analysis_metadata, expr_dir))

    with open(output_dir / Path("upload_metadata.tsv"), "w", encoding="utf-8") as upload_metadata:
        for expr, analysis in metadata_paths:
            upload_metadata.write(f"\t{expr.absolute()}\t{analysis.absolute()}\n")


def get_args():
    parser = argparse.ArgumentParser(
        description="Generate an experiment/analysis with the Engreitz Flow-FISH  CRISPRi screens from Encode"
    )
    parser.add_argument("metadata_path", help="Encode screen metadata files")
    parser.add_argument("output_directory", help="Directory to put output files")

    return parser.parse_args()


def write_experiment(experiment_metadata, analysis_metadata, output_directory):
    output_dir = Path(output_directory)
    expr_path = output_dir / "experiment.json"
    analysis_path = output_dir / "analysis001.json"
    if not expr_path.exists():
        expr_path.touch()
    if not analysis_path.exists():
        analysis_path.touch()
    expr_path.write_text(json.dumps(experiment_metadata))
    analysis_path.write_text(json.dumps(analysis_metadata))
    return expr_path, analysis_path


def run_cli():
    args = get_args()
    metadata_paths = asyncio.run(gen_data(args.metadata_path, args.output_directory))
