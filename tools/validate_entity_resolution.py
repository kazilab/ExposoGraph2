#!/usr/bin/env python
"""Entity-resolution contract for the ExposoGraph JSON migration.

Section 6 of the entity-resolution proposal. Unresolved references are
errors. The original exposure_database.json is historical and is not read.

Two readings are narrower than the proposal sentence, matching decisions
already shipped:

- index_carcinogen_id was removed and is not required. If it reappears, the
  value must be a Carcinogen node. The free-text index_carcinogen field and
  the kinetic/proxy graph_group label stay on the deprecation list until
  step 10.
- exposure_source is required only when one carcinogen_group_id and
  base_scenario match more than one scenario. Unsplit families omit it.
  The two nitrosamine general_population ids must both exist.

Run from the repository root:

    python tools/validate_entity_resolution.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "ExposoGraph" / "data"
GRAPH_PATH = ROOT / "ExposoGraph" / "map" / "graph-data.json"

REACTION_EDGE_TYPES = frozenset({"SUBSTRATE_OF", "DETOXIFIED_BY"})
FORBIDDEN_ID_FIELDS = frozenset(
    {"substrate_node_id", "enzyme_node_id", "graph_node_id"}
)
ID_FILES = (
    "kinetic_parameters.json",
    "proxy_flux_parameters.json",
    "interaction_parameters.json",
    "parameter_provenance.json",
)
NITROSAMINE_BASELINES = (
    "group_nitrosamines_tobacco_general_population",
    "group_nitrosamines_dietary_general_population",
)


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def _walk(obj: object):
    if isinstance(obj, dict):
        yield obj
        for value in obj.values():
            yield from _walk(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from _walk(value)


def _nodes_by_type(graph: dict) -> dict[str, dict[str, dict]]:
    by_type: dict[str, dict[str, dict]] = {}
    for node in graph["nodes"]:
        by_type.setdefault(node.get("type"), {})[node["id"]] = node
    return by_type


def _reaction_edges(graph: dict) -> set[tuple[str, str]]:
    return {
        (edge["source"], edge["target"])
        for edge in graph["edges"]
        if edge.get("type") in REACTION_EDGE_TYPES
    }


def _includes(graph: dict) -> set[tuple[str, str]]:
    return {
        (edge["source"], edge["target"])
        for edge in graph["edges"]
        if edge.get("type") == "INCLUDES"
    }


def _member_of(group_id: str, member_id: str, carcinogens: dict[str, dict], includes: set[tuple[str, str]]) -> bool:
    node = carcinogens.get(member_id)
    if node is None:
        return False
    if node.get("carcinogen_group_id") == group_id:
        return True
    return (group_id, member_id) in includes or (member_id, group_id) in includes


def check_id_fields(errors: list[str]) -> None:
    for name in ID_FILES:
        payload = _load(DATA / name)
        for block in _walk(payload):
            for field in FORBIDDEN_ID_FIELDS:
                if field in block:
                    errors.append(f"{name}: forbidden field {field}")
            if "index_carcinogen_id" in block:
                errors.append(
                    f"{name}: index_carcinogen_id was removed and must not be reintroduced"
                )


def check_group_and_member_ids(
    errors: list[str],
    groups: dict[str, dict],
    carcinogens: dict[str, dict],
    includes: set[tuple[str, str]],
) -> None:
    kinetic = _load(DATA / "kinetic_parameters.json")["carcinogen_classes"]
    proxy = _load(DATA / "proxy_flux_parameters.json")["classes"]
    for filename, classes in (
        ("kinetic_parameters.json", kinetic),
        ("proxy_flux_parameters.json", proxy),
    ):
        for class_id, block in classes.items():
            group_id = block.get("carcinogen_group_id")
            if group_id not in groups:
                errors.append(
                    f"{filename}:{class_id}: carcinogen_group_id {group_id!r} "
                    "is not a CarcinogenGroup"
                )
            for member_id in block.get("member_carcinogen_ids") or []:
                if member_id not in carcinogens:
                    errors.append(
                        f"{filename}:{class_id}: member {member_id!r} is not a Carcinogen"
                    )
                elif group_id in groups and not _member_of(group_id, member_id, carcinogens, includes):
                    errors.append(
                        f"{filename}:{class_id}: member {member_id!r} is not in {group_id}"
                    )

    biomarkers = _load(DATA / "biomarker_mapping.json")["entries"]
    for entry in biomarkers:
        entry_id = entry.get("entry_id") or entry.get("biomarker")
        group_id = entry.get("carcinogen_group_id")
        carcinogen_id = entry.get("carcinogen_id")
        if not group_id or group_id not in groups:
            errors.append(
                f"biomarker_mapping.json:{entry_id}: carcinogen_group_id {group_id!r} "
                "is not a CarcinogenGroup"
            )
        if carcinogen_id is not None and carcinogen_id not in carcinogens:
            errors.append(
                f"biomarker_mapping.json:{entry_id}: carcinogen_id {carcinogen_id!r} "
                "is not a Carcinogen"
            )
        elif (
            carcinogen_id
            and group_id in groups
            and not _member_of(group_id, carcinogen_id, carcinogens, includes)
        ):
            errors.append(
                f"biomarker_mapping.json:{entry_id}: {carcinogen_id} is not in {group_id}"
            )


def check_reaction_ids(errors: list[str], nodes: dict[str, str], reaction_edges: set[tuple[str, str]]) -> None:
    for name in ("kinetic_parameters.json", "proxy_flux_parameters.json"):
        payload = _load(DATA / name)
        classes = payload.get("carcinogen_classes") or payload.get("classes") or {}
        for block in _walk(classes):
            source = block.get("source_node_id")
            target = block.get("target_node_id")
            for field, value in (("source_node_id", source), ("target_node_id", target)):
                if isinstance(value, str) and value not in nodes:
                    errors.append(f"{name}: {field} {value!r} does not resolve")

    interaction = _load(DATA / "interaction_parameters.json")
    for section in ("competitive_inhibition", "phase2_conjugation", "gsh_depletion"):
        for enzyme_key, block in interaction.get(section, {}).items():
            if not isinstance(block, dict) or not isinstance(block.get("substrates"), dict):
                continue
            target = block.get("target_node_id")
            if not isinstance(target, str) or target not in nodes:
                errors.append(
                    f"interaction_parameters.json:{section}.{enzyme_key}: "
                    f"target_node_id {target!r} does not resolve"
                )
                continue
            for substrate_key, substrate in block["substrates"].items():
                if not isinstance(substrate, dict):
                    continue
                source = substrate.get("source_node_id")
                if not isinstance(source, str) or source not in nodes:
                    errors.append(
                        f"interaction_parameters.json:{section}.{enzyme_key}.{substrate_key}: "
                        f"source_node_id {source!r} does not resolve"
                    )
                    continue
                if (source, target) not in reaction_edges:
                    errors.append(
                        f"interaction_parameters.json:{section}.{enzyme_key}.{substrate_key}: "
                        f"{source} -> {target} is not SUBSTRATE_OF or DETOXIFIED_BY"
                    )


def check_profile_entities(errors: list[str], nodes: dict[str, str]) -> None:
    profiles = _load(DATA / "exposure_database_revised.json")["exposure_profiles"]
    allowed = {"Carcinogen", "CarcinogenGroup"}
    for profile_id, profile in profiles.items():
        for index, component in enumerate(profile.get("components") or []):
            entity_id = component.get("entity_id")
            node_type = nodes.get(entity_id)
            if node_type not in allowed:
                errors.append(
                    f"exposure_profiles.{profile_id}[{index}]: entity_id {entity_id!r} "
                    f"is {node_type or 'unresolved'}, not a Carcinogen or CarcinogenGroup"
                )


def check_scenarios(
    errors: list[str],
    groups: dict[str, dict],
    carcinogens: dict[str, dict],
    includes: set[tuple[str, str]],
) -> None:
    scenarios = _load(DATA / "exposure_database_revised.json")["exposure_scenarios"]
    if len(scenarios) != len(set(scenarios)):
        errors.append("exposure scenario ids are not unique")
    by_identity: dict[tuple[str, str, str | None], list[str]] = {}
    members_by_source: dict[tuple[str, str], set[str]] = {}
    covered: dict[str, set[str]] = {}
    for scenario_id, scenario in scenarios.items():
        group_id = scenario.get("carcinogen_group_id")
        if group_id not in groups:
            errors.append(
                f"{scenario_id}: carcinogen_group_id {group_id!r} is not a CarcinogenGroup"
            )
        base = scenario.get("base_scenario")
        if not base:
            errors.append(f"{scenario_id}: missing base_scenario")
        source = scenario.get("exposure_source")
        by_identity.setdefault((group_id, base, source), []).append(scenario_id)
        members = scenario.get("member_carcinogen_ids") or []
        if not members:
            errors.append(f"{scenario_id}: missing member_carcinogen_ids")
        covered.setdefault(group_id, set()).update(members)
        if source:
            members_by_source.setdefault((group_id, source), set()).update(members)
        for member_id in members:
            if member_id not in carcinogens:
                errors.append(f"{scenario_id}: member {member_id!r} is not a Carcinogen")
            elif group_id in groups and not _member_of(group_id, member_id, carcinogens, includes):
                errors.append(f"{scenario_id}: member {member_id!r} is not in {group_id}")
    for identity, scenario_ids in by_identity.items():
        if len(scenario_ids) > 1:
            errors.append(
                f"duplicate scenario identity {identity}: {', '.join(scenario_ids)}"
            )
    for (group_id, base, source), scenario_ids in by_identity.items():
        siblings = [
            other
            for (other_group, other_base, other_source), other in by_identity.items()
            if other_group == group_id and other_base == base and other_source != source
        ]
        if siblings and source in (None, ""):
            errors.append(
                f"{scenario_ids[0]}: exposure_source is required because "
                f"{group_id} / {base} is split"
            )
    for scenario_id in NITROSAMINE_BASELINES:
        if scenario_id not in scenarios:
            errors.append(f"missing nitrosamine baseline {scenario_id}")
    if len(set(NITROSAMINE_BASELINES) & set(scenarios)) == 2:
        left = scenarios[NITROSAMINE_BASELINES[0]].get("exposure_source")
        right = scenarios[NITROSAMINE_BASELINES[1]].get("exposure_source")
        if left == right:
            errors.append("nitrosamine general_population baselines do not have distinct sources")
    for group_id, members in covered.items():
        for source, source_members in (
            (source, claimed)
            for (claimed_group, source), claimed in members_by_source.items()
            if claimed_group == group_id
        ):
            other_sources = {
                other
                for (other_group, other), claimed in members_by_source.items()
                if other_group == group_id and other != source
            }
            overlap = set()
            for other in other_sources:
                overlap |= source_members & members_by_source[(group_id, other)]
            if overlap:
                errors.append(
                    f"{group_id}: members claimed under more than one exposure_source: "
                    f"{sorted(overlap)}"
                )
                break
    for group_id, node in groups.items():
        if group_id not in covered:
            continue
        expected = {
            carcinogen_id
            for carcinogen_id, carcinogen in carcinogens.items()
            if carcinogen.get("carcinogen_group_id") == group_id
        }
        missing = expected - covered[group_id]
        if missing:
            errors.append(
                f"{group_id}: exposure scenarios omit members {sorted(missing)}"
            )
        _ = node


def check_quantitative_tissues(errors: list[str], vocabulary: set[str], graph: dict) -> None:
    raw = _load(DATA / "tissue_expression_data_raw.json")
    declared = raw.get("metadata", {}).get("tissues") or []
    for key in declared:
        if key not in vocabulary:
            errors.append(f"tissue_expression_data_raw.json metadata.tissues: {key!r} is not in tissue_vocabulary")
    for gene, weights in (raw.get("expression") or {}).items():
        if not isinstance(weights, dict):
            continue
        for key in weights:
            if key not in vocabulary:
                errors.append(
                    f"tissue_expression_data_raw.json:{gene}: quantitative tissue key {key!r} "
                    "is not in tissue_vocabulary"
                )
    source_map = (raw.get("metadata") or {}).get("source_tissue_map") or {}
    for key in source_map:
        if key not in vocabulary:
            errors.append(
                f"tissue_expression_data_raw.json source_tissue_map: {key!r} "
                "is not in tissue_vocabulary"
            )
    for node in graph["nodes"]:
        weights = node.get("tissue_weights")
        if not isinstance(weights, dict):
            continue
        for key in weights:
            if key not in vocabulary:
                errors.append(
                    f"graph node {node['id']}: tissue_weights key {key!r} "
                    "is not in tissue_vocabulary"
                )


def check_rate_tables(errors: list[str], nodes: dict[str, str]) -> None:
    interaction = _load(DATA / "interaction_parameters.json")
    for table in ("reference_concentrations_uM", "gsh_reference_rates"):
        for key in interaction.get(table, {}):
            if key.startswith("_"):
                continue
            if key not in nodes:
                errors.append(f"interaction_parameters.json:{table}: {key!r} does not resolve")
    for gene, block in _load(DATA / "genotypes.json")["genes"].items():
        for phenotype_id, phenotype in block.get("phenotypes", {}).items():
            folds = phenotype.get("risk_folds") or {}
            for member_id in folds:
                if member_id not in nodes or nodes[member_id] != "Carcinogen":
                    errors.append(
                        f"genotypes.json:{gene}.{phenotype_id}: {member_id!r} is not a Carcinogen"
                    )


def main() -> int:
    graph = _load(GRAPH_PATH)
    by_type = _nodes_by_type(graph)
    nodes = {node["id"]: node.get("type") for node in graph["nodes"]}
    vocabulary = set(graph.get("tissue_vocabulary") or {})
    if len(vocabulary) != 10:
        print(f"tissue_vocabulary has {len(vocabulary)} keys, expected 10")
        return 1
    errors: list[str] = []
    check_id_fields(errors)
    check_group_and_member_ids(
        errors, by_type.get("CarcinogenGroup", {}), by_type.get("Carcinogen", {}), _includes(graph)
    )
    check_reaction_ids(errors, nodes, _reaction_edges(graph))
    check_profile_entities(errors, nodes)
    check_scenarios(
        errors, by_type.get("CarcinogenGroup", {}), by_type.get("Carcinogen", {}), _includes(graph)
    )
    check_quantitative_tissues(errors, vocabulary, graph)
    check_rate_tables(errors, nodes)
    if errors:
        for error in errors:
            print(error)
        print(f"{len(errors)} error(s)")
        return 1
    print("entity resolution contract: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
