"""Small NHANES laboratory file catalog.

This starts with a practical 2017-2018 multi-class catalog and PAH cycles from
2003-2004 onward where file keys are available. Expand this dictionary as new
classes are validated.
"""

from __future__ import annotations

NHANES_BASE_URL = "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle_suffix}/{file_key}.XPT"

CYCLE_SUFFIX = {
    "2003-2004": "2003-2004",
    "2005-2006": "2005-2006",
    "2007-2008": "2007-2008",
    "2009-2010": "2009-2010",
    "2011-2012": "2011-2012",
    "2013-2014": "2013-2014",
    "2015-2016": "2015-2016",
    "2017-2018": "2017-2018",
}

NHANES_FILE_CATALOG = {
    "2003-2004": {"PAH": "PAH_C", "DEMO": "DEMO_C", "SMQ": "SMQ_C", "UCREAT": "ALB_CR_C"},
    "2005-2006": {"PAH": "PAH_D", "DEMO": "DEMO_D", "SMQ": "SMQ_D", "UCREAT": "ALB_CR_D"},
    "2007-2008": {"PAH": "PAH_E", "DEMO": "DEMO_E", "SMQ": "SMQ_E", "UCREAT": "ALB_CR_E"},
    "2009-2010": {"PAH": "PAH_F", "DEMO": "DEMO_F", "SMQ": "SMQ_F", "UCREAT": "ALB_CR_F"},
    "2011-2012": {"PAH": "PAH_G", "DEMO": "DEMO_G", "SMQ": "SMQ_G", "UCREAT": "ALB_CR_G"},
    "2013-2014": {"PAH": "PAH_H", "DEMO": "DEMO_H", "SMQ": "SMQ_H", "UCREAT": "ALB_CR_H"},
    "2015-2016": {"PAH": "PAH_I", "DEMO": "DEMO_I", "SMQ": "SMQ_I", "UCREAT": "ALB_CR_I"},
    "2017-2018": {
        "PAH": "PAH_J",
        "COT": "COT_J",
        "PFAS": "PFAS_J",
        "PBCD": "PBCD_J",
        "UM": "UM_J",
        "PHTHTE": "PHTHTE_J",
        "EPHPP": "EPHPP_J",
        "UVOC": "UVOC_J",
        "VOCWB": "VOCWB_J",
        "OPD": "OPD_J",
        "UPHOPM": "UPHOPM_J",
        "DEMO": "DEMO_J",
        "SMQ": "SMQ_J",
        "UCREAT": "ALB_CR_J",
    },
}


def available_cycles() -> list[str]:
    return sorted(NHANES_FILE_CATALOG)


def get_cycle_files(cycle: str) -> dict[str, str]:
    try:
        return dict(NHANES_FILE_CATALOG[cycle])
    except KeyError as exc:
        raise KeyError(f"Unsupported NHANES cycle: {cycle}") from exc


def get_file_url(cycle: str, file_key: str) -> str:
    files = get_cycle_files(cycle)
    if file_key not in files:
        raise KeyError(f"File key {file_key!r} is not available for cycle {cycle}")
    return NHANES_BASE_URL.format(cycle_suffix=CYCLE_SUFFIX[cycle], file_key=files[file_key])


def class_available(cycle: str, class_name: str) -> bool:
    from .class_registry import NHANES_CLASS_REGISTRY

    info = NHANES_CLASS_REGISTRY[class_name]
    return info["default_file_key"] in NHANES_FILE_CATALOG.get(cycle, {})
