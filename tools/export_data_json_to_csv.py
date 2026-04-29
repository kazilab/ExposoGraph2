#!/usr/bin/env python3
"""Export ExposoGraph data JSON files to CSV files.

The exporter writes two types of CSVs:
- table-oriented files for common manuscript/supplementary-table use cases
- one long-form path/value CSV per JSON file so nested content is not lost
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any, Iterable


DEFAULT_DATA_DIR = Path("ExposoGraph/data")
DEFAULT_OUTPUT_DIR = DEFAULT_DATA_DIR / "csv_exports"


def scalar(value: Any) -> bool:
    return value is None or isinstance(value, (str, int, float, bool))


def stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        if all(scalar(item) for item in value):
            return "; ".join(stringify(item) for item in value)
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def flatten_record(record: dict[str, Any], prefix: str = "") -> dict[str, str]:
    flat: dict[str, str] = {}
    for key, value in record.items():
        col = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            if value and all(scalar(v) or isinstance(v, list) for v in value.values()):
                for child_key, child_value in value.items():
                    flat[f"{col}.{child_key}"] = stringify(child_value)
            else:
                flat[col] = stringify(value)
        else:
            flat[col] = stringify(value)
    return flat


def write_csv(path: Path, rows: Iterable[dict[str, Any]], preferred: Iterable[str] = ()) -> None:
    normalized = [{key: stringify(value) for key, value in row.items()} for row in rows]
    preferred_cols = list(dict.fromkeys(preferred))
    seen = set(preferred_cols)
    other_cols: list[str] = []
    for row in normalized:
        for key in row:
            if key not in seen:
                seen.add(key)
                other_cols.append(key)
    columns = preferred_cols + sorted(other_cols)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(normalized)


def walk_long(value: Any, source_file: str, path: str = "$") -> list[dict[str, str]]:
    if isinstance(value, dict):
        if not value:
            return [{"source_file": source_file, "json_path": path, "value_type": "object", "value": "{}"}]
        rows: list[dict[str, str]] = []
        for key, child in value.items():
            rows.extend(walk_long(child, source_file, f"{path}.{key}"))
        return rows
    if isinstance(value, list):
        if not value:
            return [{"source_file": source_file, "json_path": path, "value_type": "array", "value": "[]"}]
        rows = []
        for index, child in enumerate(value):
            rows.extend(walk_long(child, source_file, f"{path}[{index}]"))
        return rows
    value_type = "null" if value is None else type(value).__name__
    return [{"source_file": source_file, "json_path": path, "value_type": value_type, "value": stringify(value)}]


def export_biomarker_mapping(data: dict[str, Any], out: Path) -> list[Path]:
    path = out / "biomarker_mapping__entries.csv"
    rows = [flatten_record(row) for row in data.get("entries", [])]
    write_csv(
        path,
        rows,
        [
            "lifestyle_factor",
            "biomarker",
            "matrix",
            "reference_range",
            "reference_units",
            "partition_coefficient",
            "target_tissue",
            "target_enzyme",
            "carcinogen_class",
            "Km_uM",
            "S_over_Km_central",
            "S_over_Km_range",
            "tier_multiplier",
            "central_tissue_S_uM",
            "source_status",
            "provenance_note",
            "references",
        ],
    )
    return [path]


def export_exposure_database(data: dict[str, Any], out: Path) -> list[Path]:
    written: list[Path] = []
    classes = data.get("carcinogen_classes", {})

    summary_rows = []
    scenario_rows = []
    for class_key, class_record in classes.items():
        summary = {"carcinogen_class": class_key}
        summary.update(
            flatten_record(
                {k: v for k, v in class_record.items() if k != "exposure_scenarios"}
            )
        )
        summary_rows.append(summary)
        for scenario_id, scenario in class_record.get("exposure_scenarios", {}).items():
            row = {"carcinogen_class": class_key, "scenario_id": scenario_id}
            row.update(flatten_record(scenario))
            scenario_rows.append(row)

    path = out / "exposure_database__carcinogen_classes.csv"
    write_csv(path, summary_rows, ["carcinogen_class", "class_id", "class_label", "index_carcinogen"])
    written.append(path)

    path = out / "exposure_database__exposure_scenarios.csv"
    write_csv(path, scenario_rows, ["carcinogen_class", "scenario_id", "label", "multiplier_vs_baseline"])
    written.append(path)

    question_rows = []
    for section_id, section in data.get("exposure_questionnaire_schema", {}).get("sections", {}).items():
        for question in section.get("questions", []):
            row = {"section_id": section_id, "section_label": section.get("label", "")}
            row.update(flatten_record(question))
            question_rows.append(row)
    path = out / "exposure_database__questionnaire_questions.csv"
    write_csv(path, question_rows, ["section_id", "section_label", "id", "type", "text", "unit", "options"])
    written.append(path)

    flux = data.get("flux_model_integration", {})
    risk_rows = []
    for threshold_id, threshold in flux.get("risk_thresholds", {}).items():
        row = {"threshold_id": threshold_id}
        row.update(flatten_record(threshold))
        risk_rows.append(row)
    path = out / "exposure_database__risk_thresholds.csv"
    write_csv(path, risk_rows, ["threshold_id", "label", "min", "max", "color"])
    written.append(path)

    modifier_rows = []
    for modifier_id, modifier in flux.get("genotype_modifiers", {}).items():
        row = {"modifier_id": modifier_id}
        if isinstance(modifier, dict):
            row.update(flatten_record(modifier))
        else:
            row["value"] = modifier
        modifier_rows.append(row)
    path = out / "exposure_database__genotype_modifiers.csv"
    write_csv(path, modifier_rows, ["modifier_id"])
    written.append(path)

    return written


def export_interaction_parameters(data: dict[str, Any], out: Path) -> list[Path]:
    written: list[Path] = []

    term_rows = []
    for section_name in ("competitive_inhibition", "phase2_conjugation"):
        for enzyme, enzyme_record in data.get(section_name, {}).items():
            if enzyme.startswith("_") or not isinstance(enzyme_record, dict):
                continue
            for term_type in ("substrates", "inhibitors", "induction"):
                terms = enzyme_record.get(term_type, {})
                if not isinstance(terms, dict):
                    continue
                for term_id, term_record in terms.items():
                    row = {"section": section_name, "enzyme": enzyme, "term_type": term_type, "term_id": term_id}
                    if isinstance(term_record, dict):
                        row.update(flatten_record(term_record))
                    else:
                        row["value"] = term_record
                    term_rows.append(row)
    path = out / "interaction_parameters__enzyme_terms.csv"
    write_csv(path, term_rows, ["section", "enzyme", "term_type", "term_id", "Km_uM", "Ki_uM", "Vmax_relative"])
    written.append(path)

    induction_rows = []
    for induction_id, induction in data.get("enzyme_induction", {}).items():
        if induction_id.startswith("_"):
            continue
        row = {"induction_id": induction_id}
        if isinstance(induction, dict):
            row.update(flatten_record(induction))
        else:
            row["value"] = induction
        induction_rows.append(row)
    path = out / "interaction_parameters__enzyme_induction.csv"
    write_csv(path, induction_rows, ["induction_id", "mechanism", "fold_induction", "fold"])
    written.append(path)

    consumer_rows = []
    for consumer_id, consumer in data.get("gsh_depletion", {}).get("consumers", {}).items():
        row = {"consumer_id": consumer_id}
        row.update(flatten_record(consumer))
        consumer_rows.append(row)
    path = out / "interaction_parameters__gsh_consumers.csv"
    write_csv(path, consumer_rows, ["consumer_id", "substrate_class", "enzyme", "gsh_per_umol_substrate"])
    written.append(path)

    genotype_rows = []
    for gene, variants in data.get("genotype_modifiers", {}).items():
        if not isinstance(variants, dict):
            continue
        for variant_id, variant in variants.items():
            row = {"gene": gene, "variant_id": variant_id}
            if isinstance(variant, dict):
                row.update(flatten_record(variant))
            else:
                row["value"] = variant
            genotype_rows.append(row)
    path = out / "interaction_parameters__genotype_modifiers.csv"
    write_csv(path, genotype_rows, ["gene", "variant_id", "activity_multiplier", "frequency", "alleles"])
    written.append(path)

    rule_rows = []
    for rule_type, rules in data.get("interaction_rules", {}).items():
        if isinstance(rules, list):
            for index, rule in enumerate(rules, start=1):
                row = {"rule_type": rule_type, "rule_index": index}
                row.update(flatten_record(rule))
                rule_rows.append(row)
    path = out / "interaction_parameters__interaction_rules.csv"
    write_csv(path, rule_rows, ["rule_type", "rule_index", "name", "mechanism", "synergy_score", "antagonism_score"])
    written.append(path)

    return written


def export_kinetic_parameters(data: dict[str, Any], out: Path) -> list[Path]:
    written: list[Path] = []
    classes = data.get("carcinogen_classes", {})

    summary_rows = []
    pathway_rows = []
    validation_rows = []
    for class_key, class_record in classes.items():
        summary = {"carcinogen_class": class_key}
        summary.update(
            flatten_record(
                {
                    k: v
                    for k, v in class_record.items()
                    if k not in {"pathways", "epidemiological_validation"}
                }
            )
        )
        summary_rows.append(summary)

        for pathway_type, pathway_record in class_record.get("pathways", {}).items():
            if isinstance(pathway_record, dict):
                for term_id, term_record in pathway_record.items():
                    row = {"carcinogen_class": class_key, "pathway_type": pathway_type, "term_id": term_id}
                    if isinstance(term_record, dict):
                        row.update(flatten_record(term_record))
                    else:
                        row["value"] = term_record
                    pathway_rows.append(row)
            else:
                pathway_rows.append(
                    {"carcinogen_class": class_key, "pathway_type": pathway_type, "term_id": "", "value": pathway_record}
                )

        for validation_id, validation in class_record.get("epidemiological_validation", {}).items():
            row = {"carcinogen_class": class_key, "validation_id": validation_id}
            if isinstance(validation, dict):
                row.update(flatten_record(validation))
            else:
                row["value"] = validation
            validation_rows.append(row)

    path = out / "kinetic_parameters__carcinogen_classes.csv"
    write_csv(path, summary_rows, ["carcinogen_class", "index_carcinogen", "description", "confidence_overall"])
    written.append(path)

    path = out / "kinetic_parameters__pathways.csv"
    write_csv(path, pathway_rows, ["carcinogen_class", "pathway_type", "term_id", "reaction", "Km_uM", "Vmax_relative", "confidence"])
    written.append(path)

    path = out / "kinetic_parameters__epidemiological_validation.csv"
    write_csv(path, validation_rows, ["carcinogen_class", "validation_id", "value", "source"])
    written.append(path)

    special_rows = []
    for modifier_id, modifier in data.get("genotype_modifiers", {}).get("special_cases", {}).items():
        row = {"modifier_id": modifier_id}
        row.update(flatten_record(modifier))
        special_rows.append(row)
    path = out / "kinetic_parameters__genotype_special_cases.csv"
    write_csv(path, special_rows, ["modifier_id", "phenotype", "activity_fraction", "mechanism", "source"])
    written.append(path)

    weight_rows = []
    for gene, weights in data.get("tissue_expression_weights", {}).items():
        if gene == "description":
            continue
        row = {"gene": gene}
        row.update(flatten_record(weights) if isinstance(weights, dict) else {"value": weights})
        weight_rows.append(row)
    path = out / "kinetic_parameters__tissue_expression_weights.csv"
    write_csv(path, weight_rows, ["gene"])
    written.append(path)

    return written


def export_mutational_signatures(data: dict[str, Any], out: Path) -> list[Path]:
    written: list[Path] = []
    rows = []
    for signature_id, signature in data.get("signatures", {}).items():
        row = {"signature_id": signature_id}
        row.update(flatten_record(signature))
        rows.append(row)
    path = out / "mutational_signatures__signatures.csv"
    write_csv(path, rows, ["signature_id", "description", "aetiology", "dominant_mutation", "confidence", "references"])
    written.append(path)

    rows = []
    for class_key, mapping in data.get("carcinogen_class_map", {}).items():
        row = {"carcinogen_class": class_key}
        row.update(flatten_record(mapping))
        rows.append(row)
    path = out / "mutational_signatures__carcinogen_class_map.csv"
    write_csv(path, rows, ["carcinogen_class", "primary", "secondary", "notes"])
    written.append(path)
    return written


def export_parameter_provenance(data: dict[str, Any], out: Path) -> list[Path]:
    rows = []
    for enzyme, substrates in data.get("pairs", {}).items():
        for substrate, provenance in substrates.items():
            row = {"enzyme": enzyme, "substrate": substrate}
            row.update(flatten_record(provenance))
            rows.append(row)
    path = out / "parameter_provenance__pairs.csv"
    write_csv(path, rows, ["enzyme", "substrate", "km_source", "km_confidence", "vmax_source", "vmax_confidence", "ki_status", "ki_value_uM", "ki_reference"])
    return [path]


def export_proxy_flux(data: dict[str, Any], out: Path, stem: str) -> list[Path]:
    written: list[Path] = []
    classes = data.get("classes", {})

    summary_rows = []
    term_rows = []
    for class_key, class_record in classes.items():
        summary = {"class": class_key}
        summary.update(
            flatten_record(
                {
                    k: v
                    for k, v in class_record.items()
                    if k not in {"activation_terms", "detox_terms", "repair_terms"}
                }
            )
        )
        summary_rows.append(summary)
        for term_type in ("activation_terms", "detox_terms", "repair_terms"):
            terms = class_record.get(term_type, {})
            if not isinstance(terms, dict):
                continue
            for term_id, term_record in terms.items():
                row = {"class": class_key, "term_type": term_type, "term_id": term_id}
                if isinstance(term_record, dict):
                    row.update(flatten_record(term_record))
                else:
                    row["value"] = term_record
                term_rows.append(row)

    path = out / f"{stem}__classes.csv"
    write_csv(path, summary_rows, ["class", "model_kind", "signal", "exposure_default", "unit_note"])
    written.append(path)

    path = out / f"{stem}__terms.csv"
    write_csv(path, term_rows, ["class", "term_type", "term_id", "gene", "equation", "km", "vmax", "confidence", "provenance_ref", "parameter_basis", "sources"])
    written.append(path)
    return written


def export_tissue_expression(data: dict[str, Any], out: Path) -> list[Path]:
    written: list[Path] = []
    for section in ("expression", "weights"):
        rows = []
        for gene, values in data.get(section, {}).items():
            row = {"gene": gene}
            row.update(flatten_record(values) if isinstance(values, dict) else {"value": values})
            rows.append(row)
        path = out / f"tissue_expression_data__{section}.csv"
        write_csv(path, rows, ["gene", "Liver", "Lung", "Prostate", "Bladder", "Colon", "Breast", "Kidney", "Esophagus"])
        written.append(path)
    return written


SPECIAL_EXPORTERS = {
    "biomarker_mapping": export_biomarker_mapping,
    "exposure_database": export_exposure_database,
    "interaction_parameters": export_interaction_parameters,
    "kinetic_parameters": export_kinetic_parameters,
    "mutational_signatures": export_mutational_signatures,
    "parameter_provenance": export_parameter_provenance,
    "tissue_expression_data": export_tissue_expression,
}


def export_all(data_dir: Path, output_dir: Path) -> list[Path]:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    written: list[Path] = []
    manifest_rows = []
    for json_path in sorted(data_dir.glob("*.json")):
        data = json.loads(json_path.read_text(encoding="utf-8"))
        stem = json_path.stem

        long_path = output_dir / f"{stem}__long.csv"
        write_csv(long_path, walk_long(data, json_path.name), ["source_file", "json_path", "value_type", "value"])
        written.append(long_path)
        manifest_rows.append({"source_json": json_path.name, "csv_file": long_path.name, "export_type": "long_path_value"})

        if stem in SPECIAL_EXPORTERS:
            paths = SPECIAL_EXPORTERS[stem](data, output_dir)
        elif stem in {"proxy_flux_parameters", "proxy_flux_provenance"}:
            paths = export_proxy_flux(data, output_dir, stem)
        else:
            paths = []
        written.extend(paths)
        for path in paths:
            manifest_rows.append({"source_json": json_path.name, "csv_file": path.name, "export_type": "table"})

    manifest_path = output_dir / "manifest.csv"
    write_csv(manifest_path, manifest_rows, ["source_json", "csv_file", "export_type"])
    written.append(manifest_path)

    readme_path = output_dir / "README.md"
    readme_path.write_text(
        "# ExposoGraph Data CSV Exports\n\n"
        "Generated from `ExposoGraph/data/*.json` by `tools/export_data_json_to_csv.py`.\n\n"
        "- `*__long.csv` files preserve every scalar JSON value as `json_path,value_type,value` rows.\n"
        "- Other `*__.csv` files are table-oriented exports for manuscript and supplementary-table use.\n"
        "- Regenerate with: `python tools/export_data_json_to_csv.py`.\n",
        encoding="utf-8",
    )
    written.append(readme_path)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    written = export_all(args.data_dir, args.output_dir)
    print(f"Wrote {len(written)} files to {args.output_dir}")


if __name__ == "__main__":
    main()
