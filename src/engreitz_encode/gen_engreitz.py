import argparse
import asyncio
import json
import os.path
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import httpx

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


def first(iterable, test):
    for x in iterable:
        if test(x):
            return x


def write_experiment_data():
    pass


def write_analysis_data():
    pass


def normalize_assembly(assembly):
    if assembly in GRCH38_ASSEMBLIES:
        return GRCH38
    elif assembly in GRCH37_ASSEMBLIES:
        return GRCH37

    raise ValueError(f"Invalid assembly {assembly}")


def get_source_type(element_reference):
    source_types = element_reference.get("elements_selection_method")
    if source_types is None:
        return "Tested Element"

    st_set = set(source_types)
    if "accessible genome regions" in st_set:
        return "Chromatin Accessible Region"
    elif "candidate cis-regulatory elements" in st_set:
        return "cCRE"
    elif "DNase hypersensitive sites" in st_set:
        return "DHS"
    else:
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


def build_experiment(screen_info):
    gene_assembly = normalize_assembly(get_assembly(screen_info))

    cell_line = screen_info["biosample_ontology"][0]["term_name"]
    tissue_type = TISSUE_TYPES[cell_line]
    crispr = screen_info["related_datasets"][0]["perturbation_type"]

    create_date = datetime.fromisoformat(screen_info["date_created"])

    af = get_analysis_files(screen_info)
    q_file = first(af, lambda x: x["output_category"] == "quantification")

    ref_files = screen_info["elements_references"][0]["files"]
    guides_file = [f for f in ref_files if "guide" in f["aliases"][0]][0]
    elements_file = [f for f in ref_files if "element" in f["aliases"][0]][0]

    guide_quant = first(
        screen_info["related_datasets"][0]["files"],
        lambda x: x["output_type"] == "guide quantifications",
    )

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
            "filename": os.path.basename(guides_file["cloud_metadata"]["url"]),
            "file_location": guides_file["cloud_metadata"]["url"],
            "genome_assembly": gene_assembly,
        },
        "misc_files": [
            {
                "description": elements_file["aliases"][0],
                "filename": os.path.basename(elements_file["cloud_metadata"]["url"]),
                "file_location": elements_file["cloud_metadata"]["url"],
                "genome_assembly": gene_assembly,
            },
            {
                "description": q_file["output_type"].capitalize(),
                "filename": os.path.basename(q_file["cloud_metadata"]["url"]),
                "file_location": q_file["cloud_metadata"]["url"],
                "genome_assembly": gene_assembly,
            },
            {
                "description": guide_quant["output_type"].capitalize(),
                "filename": os.path.basename(guide_quant["cloud_metadata"]["url"]),
                "file_location": guide_quant["cloud_metadata"]["url"],
                "genome_assembly": gene_assembly,
            },
        ],
        "data_format": "jesse-engreitz",
    }

    return experiment


def build_analysis(screen_info):
    gene_assembly = normalize_assembly(get_assembly(screen_info))

    af = get_analysis_files(screen_info)
    q_file = first(af, lambda x: x["output_category"] == "quantification")
    gene = screen_info["examined_loci"][0]["gene"]["symbol"]
    crispr = screen_info["related_datasets"][0]["perturbation_type"]
    cell_line = screen_info["biosample_ontology"][0]["term_name"]

    ref_files = screen_info["elements_references"][0]["files"]
    guides_file = [f for f in ref_files if "guide" in f["aliases"][0]][0]
    elements_file = [f for f in ref_files if "element" in f["aliases"][0]][0]

    guide_quant = first(
        screen_info["related_datasets"][0]["files"],
        lambda x: x["output_type"] == "guide quantifications",
    )

    analysis = {
        "name": f"{crispr} {screen_info['assay_title'][0]} in {cell_line} for {gene} ({screen_info['accession']})",
        "description": screen_info.get("biosample_summary"),
        "source type": "gRNA",
        "genome_assembly": gene_assembly,
        "p_val_adj_method": "Benjamini-Hochberg",
        "p_val_threshold": 0.05,
        "results": {
            "description": f"{q_file['output_type'].capitalize()}\nFrom Ben Doughty, via Email:\nFor the effect size calculations, since we do a 6-bin sort, we don't actually compute a log2-fold change. Instead, we use the data to compute an effect size in \"gene expression\" space, which we normalize to the negative controls. The values are then scaled, so what an effect size of -0.2 means is that this guide decreased the expression of the target gene by 20%. An effect size of 0 would be no change in expression, and an effect size of +0.1 would mean a 10% increase in expression. The lowest we can go is -1 (which means total elimination of signal), and technically the effect size is unbounded in the positive direction, although we never see _super_ strong positive guides with CRISPRi. ",
            "filename": os.path.basename(q_file["cloud_metadata"]["url"]),
            "file_location": q_file["cloud_metadata"]["url"],
        },
        "misc_files": [
            {
                "description": guides_file["aliases"][0],
                "filename": os.path.basename(guides_file["cloud_metadata"]["url"]),
                "file_location": guides_file["cloud_metadata"]["url"],
            },
            {
                "description": elements_file["aliases"][0],
                "filename": os.path.basename(elements_file["cloud_metadata"]["url"]),
                "file_location": elements_file["cloud_metadata"]["url"],
            },
            {
                "description": guide_quant["output_type"].capitalize(),
                "filename": os.path.basename(guide_quant["cloud_metadata"]["url"]),
                "file_location": guide_quant["cloud_metadata"]["url"],
            },
        ],
        "data_format": "jesse-engreitz",
    }

    if (desc := screen_info.get("description")) is not None:
        analysis["description"] = desc

    return analysis


# only one file
def get_guides_file_url(screens) -> str:
    screen = screens[0]
    ref_files = screen["elements_references"][0]["files"]
    return [f for f in ref_files if "guide" in f["aliases"][0]][0]["cloud_metadata"]["url"]


# only one file
def get_elements_file_url(screens) -> str:
    screen = screens[0]
    ref_files = screen["elements_references"][0]["files"]
    return [f for f in ref_files if "element" in f["aliases"][0]][0]["cloud_metadata"]["url"]


def get_element_quantification_urls(screens) -> list[str]:
    element_quants = []
    for screen in screens:
        af = get_analysis_files(screen)
        q_file = first(af, lambda x: x["output_category"] == "quantification")
        element_quants.append(q_file["cloud_metadata"]["url"])

    return element_quants


def get_guide_quantification_urls(screens) -> list[str]:
    return [
        first(
            screen["related_datasets"][0]["files"],
            lambda x: x["output_type"] == "guide quantifications",
        )[
            "cloud_metadata"
        ]["url"]
        for screen in screens
    ]


async def download_file(client, url, output_path):
    response = await client.get(url, timeout=5)
    with open(output_path, "w", encoding="utf-8") as output:
        output.write(response.content.decode())


async def gen_metadata(metadata_path, output_path) -> tuple[Optional[dict], Optional[dict]]:
    curr_dir = Path(metadata_path)
    json_files = list(curr_dir.glob("*.json"))

    output_dir = Path(output_path)
    output_dir.mkdir(exist_ok=True)

    experiments = []
    for file in json_files:
        with open(file, encoding="utf-8") as f:
            screen = json.load(f)

        if screen["accession"] not in GOOD_ENGREITZ_EXPERIMENTS:
            continue

        experiments.append((file, screen))

    if len(experiments) == 0:
        return None, None

    screens = [expr[1] for expr in experiments]

    with tempfile.TemporaryDirectory(dir=output_path) as temp_dir:
        downloads = []
        guides_url = get_guides_file_url(screens)
        guides_file = temp_dir / Path(os.path.basename(guides_url))
        downloads.append((guides_url, guides_file))

        elements_url = get_elements_file_url(screens)
        elements_file = temp_dir / Path(os.path.basename(elements_url))
        downloads.append((elements_url, elements_file))

        element_quant_urls = get_element_quantification_urls(screens)
        element_quant_files = [temp_dir / Path(os.path.basename(qu)) for qu in element_quant_urls]
        downloads.extend(zip(element_quant_urls, element_quant_files))

        guide_quant_urls = get_guide_quantification_urls(screens)
        guide_quant_files = [temp_dir / Path(os.path.basename(qu)) for qu in guide_quant_urls]
        downloads.extend(zip(guide_quant_urls, guide_quant_files))

        async with httpx.AsyncClient() as client:
            tasks = [asyncio.create_task(download_file(client, url, file)) for url, file in downloads]
            await asyncio.wait(tasks)

        # for file, screen in experiments:
        #     new_dir = temp_dir / Path(file.stem)
        #     new_dir.mkdir(exist_ok=True)
        #     build_experiment(screen)
        #     build_analysis(screen)

    return {}, {}


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


def run_cli():
    args = get_args()
    experiment_metadata, analysis_metadata = asyncio.run(gen_metadata(args.metadata_path, args.output_directory))
    if experiment_metadata is None or analysis_metadata is None:
        return
    write_experiment(experiment_metadata, analysis_metadata, args.output_directory)
