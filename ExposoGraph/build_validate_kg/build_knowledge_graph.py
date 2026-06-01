import json
import build_helpers as bh
from pydantic import TypeAdapter, ValidationError, Discriminator
from typing import Dict, List, Annotated, Union
import networkx as nx
import validation_classes as valid_classes


def build_node_dict(node_metadata_list, node_adapter):
    validated_nodes = node_adapter.validate_python(node_metadata_list)
    node_dict = {node.id: node for node in validated_nodes}
    return node_dict


def build_edge_dict(edge_metadata_list, edge_adapter):
    validated_edges = edge_adapter.validate_python(edge_metadata_list)
    edgelist = {(edge.source, edge.target): edge for edge in validated_edges}
    return edgelist


graph_data = "../map/graph-data.js"
inter_path = "../data/interaction_parameters.json"
raw_data = bh.parse_static_js_const(graph_data, "GRAPH_DATA")


node_dict = build_node_dict(raw_data["nodes"], valid_classes.node_adapter)
edge_dict = build_edge_dict(raw_data["edges"], valid_classes.edge_adapter)

net = nx.from_edgelist(edge_dict.keys())

with open(inter_path) as f:
    interaction_params = json.load(f)

genotype_modifiers_metadata = interaction_params["genotype_modifiers"]
genotype_modifiers_metadata.pop("_description")

# fix later, dumb to run through twice
genotype_class = valid_classes.GenotypeModifiersContainer(
    **genotype_modifiers_metadata
).model_dump()

for enzyme, gene_mod in genotype_class.items():
    node_dict[enzyme].genotype_modifiers = gene_mod


competitive_inhibition_metadata = interaction_params["competitive_inhibition"]
competitive_inhibition_metadata.pop("_description")
enzyme_substrate_comp_inhib_dict = {}
for enzyme_name, inter_metadata in competitive_inhibition_metadata.items():
    inter_metadata.pop("_description")
    for substrate_name, metadata_dict in inter_metadata["substrates"].items():
        enzyme_substrate_comp_inhib_dict[(enzyme_name, substrate_name)] = (
            valid_classes.CompetitiveInhibitionSubstrate(**metadata_dict)
        )
    if "inhibitors" in inter_metadata:
        pass  # come back to
    if "substrate_inhibition" in inter_metadata:
        pass  # come back to
    if "induction" in inter_metadata:
        pass  # come back to

# need to add this to edges, but as invisible
# we also need to some kind of type to this edge

phase2_conjugation_metadata = interaction_params["phase2_conjugation"]
phase2_conjugation_metadata.pop("_description")
enzyme_substrate_phase2_conj_dict = {}
for enzyme_name, metadata in phase2_conjugation_metadata.items():
    metadata.pop("_description")
    for substrate_name, metadata_dict in metadata["substrates"].items():
        enzyme_substrate_phase2_conj_dict[(enzyme_name, substrate_name)] = (
            valid_classes.Phase2Conjugation(**metadata_dict)
        )

print(enzyme_substrate_phase2_conj_dict)
