import sys
from enum import StrEnum


class ReddylabModule(StrEnum):
    KLANN_SCCERES = "klann_2021_scceres"
    KLANN_WGCERES = "klann_2021_wgceres"
    BOUNDS_SCCERES = "bounds_2021_scceres"
    SIKLENKA_STARRSEQ = "siklinka_2022_atacstarrseq"
    MCCUTCHEON_SCCERES = "mccutcheon_2022_scceres"


def run_cli():
    if len(sys.argv) <= 2:
        sys.exit()

    match sys.argv[1]:
        case ReddylabModule.KLANN_SCCERES:
            from reddylab.klann_2021_scceres import run_cli as rc
        case ReddylabModule.KLANN_WGCERES:
            from reddylab.klann_2021_wgceres import run_cli as rc
        case ReddylabModule.BOUNDS_SCCERES:
            from reddylab.bounds_2021_scceres import run_cli as rc
        case ReddylabModule.SIKLENKA_STARRSEQ:
            from reddylab.siklinka_2022_atacstarrseq import run_cli as rc
        case ReddylabModule.MCCUTCHEON_SCCERES:
            from reddylab.mccutcheon_2022_scceres import run_cli as rc
        case _:
            modules = "".join([f"\n  {module}" for module in list(ReddylabModule)])
            sys.exit(f"Invalid module. Please use one of {modules}")
    print(sys.argv[1])
    rc()
