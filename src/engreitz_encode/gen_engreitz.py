import argparse
import json
import os.path
from datetime import datetime
from pathlib import Path

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


def gen_metadata():
    curr_dir = Path("./temp_data")
    json_files = list(curr_dir.glob("*.json"))

    experiments = []
    for file in json_files:
        with open(file, encoding="utf-8") as f:
            screen = json.load(f)

        if screen["accession"] not in GOOD_ENGREITZ_EXPERIMENTS:
            continue

        experiments.append((file, screen))

    screens = set()
    biosamples = set()
    classifications = set()

    for file3, screen3 in experiments:
        classifications.add(screen3["biosample_ontology"][0]["classification"])
        biosamples.add(screen3["biosample_ontology"][0]["term_name"])
        screens.add(screen3["assay_term_name"][0])
        new_dir = curr_dir / Path(file3.stem)
        new_dir.mkdir(exist_ok=True)
        expr_path = new_dir / "experiment.json"
        analysis_path = new_dir / "analysis001.json"
        if not expr_path.exists():
            expr_path.touch()
        if not analysis_path.exists():
            analysis_path.touch()
        expr_path.write_text(json.dumps(build_experiment(screen3)))
        analysis_path.write_text(json.dumps(build_analysis(screen3)))

    print(biosamples)
    print(classifications)
    print(screens)


def get_args():
    parser = argparse.ArgumentParser(description="")

    return parser.parse_args()


def run_cli():
    args = get_args
    gen_metadata()
