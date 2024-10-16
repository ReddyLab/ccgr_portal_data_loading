from datetime import datetime
from functools import cached_property
from pathlib import Path
from typing import Optional


GRCH37 = "GRCh37"
GRCH38 = "GRCh38"
HG19 = "hg19"
HG38 = "hg38"
GRCH37_ASSEMBLIES = {GRCH37, HG19}
GRCH38_ASSEMBLIES = {GRCH38, HG38}
TISSUE_TYPES = {
    "K562": "Bone Marrow",
    "MCF-7": "Breast",
    "HCT116": "Colon",
    "DU 145": "Prostate",
    "A549": "Lung",
    "NCI-H460": "Lung",
    "PC-3": "Prostate",
    "MDA-MB-231": "Breast",
    "HepG2": "Liver",
    "SW620": "Colon",
}


def first(iterable, test):
    for x in iterable:
        if test(x):
            return x


class ScreenMetadata:
    def __init__(self, screen_info):
        self.screen_info = screen_info

    @cached_property
    def name(self):
        return f"{self.functional_characterization} {self.screen_info['assay_title'][0]} in {self.cell_line}"

    @cached_property
    def description(self):
        return self.screen_info.get("description")

    @cached_property
    def summary(self):
        return self.screen_info.get("biosample_summary")

    @cached_property
    def cell_line(self):
        return self.screen_info["biosample_ontology"][0]["term_name"]

    @cached_property
    def tissue_type(self):
        return TISSUE_TYPES[self.cell_line]

    @cached_property
    def assembly(self):
        if len(self.screen_info["assembly"]) > 0:
            assembly = self.screen_info["assembly"][0]
        else:
            assembly = self.screen_info["elements_references"][0]["assembly"][0]

        if assembly in GRCH37_ASSEMBLIES:
            return HG19

        if assembly in GRCH38_ASSEMBLIES:
            return HG38

        raise ValueError(f"Invalid assembly {assembly}")

        return assembly

    @cached_property
    def assay(self):
        return self.screen_info["assay_title"][0]

    @cached_property
    def lab(self):
        return self.screen_info["lab"]["title"]

    @cached_property
    def functional_characterization(self):
        try:
            return self.screen_info["related_datasets"][0]["perturbation_type"]
        except:
            return None

    @property
    def source_type(self) -> str:
        return "gRNA"

    @cached_property
    def parent_source_type(self) -> Optional[str]:
        return self._get_source_type(self.screen_info["elements_references"][0])

    @cached_property
    def tested_elements_description(self):
        ref_files = self.screen_info["elements_references"][0]["files"]
        guides_file = [f for f in ref_files if "guide" in f["aliases"][0]][0]

        return guides_file["aliases"][0]

    @cached_property
    def results_description(self) -> str:
        af = self._get_analysis_files()
        q_file = first(af, lambda x: x["output_category"] == "quantification")
        header = "" if q_file is None else f"{q_file['output_type'].capitalize()}\n"

        return f'{header}From Ben Doughty, via Email:\nFor the effect size calculations, since we do a 6-bin sort, we don\'t actually compute a log2-fold change. Instead, we use the data to compute an effect size in "gene expression" space, which we normalize to the negative controls. The values are then scaled, so what an effect size of -0.2 means is that this guide decreased the expression of the target gene by 20%. An effect size of 0 would be no change in expression, and an effect size of +0.1 would mean a 10% increase in expression. The lowest we can go is -1 (which means total elimination of signal), and technically the effect size is unbounded in the positive direction, although we never see _super_ strong positive guides with CRISPRi. '

    @property
    def p_val_adj_method(self):
        return "Benjamini-Hochberg"

    @property
    def p_val_threshold(self):
        return 0.05

    def _get_source_type(self, element_reference):
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

    def _get_analysis_files(self):
        af_accessions = self.screen_info["analyses"][0]["files"]
        return [file for file in self.screen_info["files"] if file["@id"] in af_accessions]

    def build_experiment(self, output_dir: Path):
        create_date = datetime.fromisoformat(self.screen_info["date_created"])

        experiment = {
            "name": self.name,
            "description": self.summary,
            "biosamples": [{"cell_type": self.cell_line, "tissue_type": self.tissue_type}],
            "assay": self.assay,
            "source type": self.source_type,
            "year": str(create_date.year),
            "lab": self.lab,
            "functional_characterization_modality": self.functional_characterization,
            "tested_elements_file": {
                "description": self.tested_elements_description,
                "filename": TESTED_ELEMENTS_FILE,
                "file_location": str((output_dir / Path(TESTED_ELEMENTS_FILE)).absolute()),
                "genome_assembly": self.assembly,
            },
        }

        if self.parent_source_type is not None:
            experiment["parent source type"] = self.parent_source_type

        return experiment

    def build_analysis(self, output_dir: Path):
        analysis = {
            "name": self.name,
            "description": self.summary,
            "source type": self.source_type,
            "genome_assembly": self.assembly,
            "p_val_adj_method": self.p_val_adj_method,
            "p_val_threshold": self.p_val_threshold,
            "results": {
                "description": self.results_description,
                "filename": OBSERVATIONS_FILE,
                "file_location": str((output_dir / Path(OBSERVATIONS_FILE)).absolute()),
            },
        }

        if (desc := self.description) is not None:
            analysis["description"] = desc

        return analysis
