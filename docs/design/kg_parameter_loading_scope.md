# Scoping note: loading kinetic/interaction/tissue parameters into the knowledge graph

Status: **scoping only — no implementation in this commit.** This documents
options for a follow-up change; it intentionally does not modify any code.

## Problem

`flux_engine.py`, `exposure_engine.py`, and `interaction_engine.py` — the
three main computational modules — never import `GraphEngine` or
`KnowledgeGraph` at all. Each module independently lazy-loads and caches its
own copy of `data/kinetic_parameters.json`, `data/interaction_parameters.json`,
`data/exposure_database.json`, `data/proxy_flux_parameters.json`, and
`data/tissue_expression_data.json` at module scope (`_KINETIC_CACHE`,
`_INTERACTION_CACHE`, `_EXPOSURE_DB_CACHE`, ...). The bundled reference graph
(`GraphEngine` / `map/graph-data.json`, see the JSON-migration commit) is a
separate object that these modules never touch. The one exception is
`tissue_expression_data.json`, which is *sometimes* draped onto graph nodes
after the fact via `tissue_subgraphs.annotate_graph_with_tissue_weights`, but
that is an opt-in post-processing step, not something the graph carries by
default.

Net effect: parameters that describe the same real-world entities (a gene, an
enzyme, a carcinogen, a tissue) live in up to four disconnected places, and
each computational module re-implements its own lookup/caching logic instead
of asking one graph object "what do you know about this edge/node."

## What the schema already gives us for free

Before inventing new structure, it's worth noting the existing `models.py`
schema already anticipates most of this:

- `Edge.kinetics: Optional[dict[str, Any]]` — an open-ended per-edge payload
  already exists and is currently unused by any bundled edge. This is a
  natural home for Km/Vmax/CLint/confidence (`kinetic_parameters.json`) and
  for Ki / competitive-inhibition terms (`interaction_parameters.json`).
- `Edge.provenance: list[ProvenanceRecord]` already exists on every edge, so
  citation-style entries from `parameter_provenance.json` /
  `proxy_flux_provenance.json` can ride along on the same edge instead of
  living in a separately-loaded provenance JSON.
- `Edge.carcinogen: Optional[str]` is already used as disambiguating context
  on non-carcinogen edges (e.g. which carcinogen a Gene→Metabolite edge's
  activation pertains to). The same pattern extends naturally to "which
  competing substrate does this Ki apply against."
- `NodeType.TISSUE` and `EdgeType.EXPRESSED_IN` (and `ENCODES`, `CUSTOM`) are
  already defined in the enums but **unused** by the bundled reference graph
  — zero nodes of type `Tissue`, zero `EXPRESSED_IN` edges exist today. Tissue
  integration is a data-population problem, not a schema change.
- `Node.tissue_weights: Optional[dict[str, float]]` already exists on
  `Node`, currently populated only by the opt-in
  `annotate_graph_with_tissue_weights` helper.

## Proposed mapping, JSON → graph structure

| Source JSON | Proposed home | New schema needed? |
|---|---|---|
| `tissue_expression_data.json` | New `Tissue` nodes (one per GTEx tissue) + `EXPRESSED_IN` edges from Gene/Enzyme nodes, weight carried on the edge (or continue using `Node.tissue_weights` as a denormalized cache) | No — `NodeType.TISSUE`/`EdgeType.EXPRESSED_IN` already exist, unused |
| `kinetic_parameters.json` (`carcinogen_classes.*.pathways.{activation,detox}.<enzyme>`) | `Edge.kinetics` dict on the existing `ACTIVATES`/`DETOXIFIES` edges between a carcinogen (or a new class-level grouping node — see below) and the enzyme | Possibly a new `EdgeType.MEMBER_OF`/node type if class-level defaults need to be shared across multiple specific carcinogens (see open question below) |
| `interaction_parameters.json` (`competitive_inhibition`, `enzyme_induction`) | `Edge.kinetics` dict on `INHIBITS`/`INDUCES` edges (both already used in the bundled graph), using `Edge.carcinogen` as the "competing substrate" context | No |
| `interaction_parameters.json` (`gsh_depletion`) | `Edge.kinetics` on a carcinogen/metabolite → `Pathway` or new small-molecule node (GSH already likely exists as a `Metabolite` node — needs verification) | No, if a GSH node already exists |
| `interaction_parameters.json` (`genotype_modifiers`) / `kinetic_parameters.json` genotype terms | New auxiliary node type per enzyme allele/genotype, linked to the Enzyme node with a new edge carrying the activity multiplier | Yes — new node type + new edge type (naming TBD, e.g. `Genotype`/`Allele` and `MODIFIES_ACTIVITY`; not committing to these names here) |
| `parameter_provenance.json` / `proxy_flux_provenance.json` | Folded into the relevant edge's existing `provenance` list instead of a separately loaded file | No |

The genotype-modifier row is the one case that likely needs a genuinely new
node type; everything else fits into fields the schema already has. Whatever
new node/edge types end up being needed, the point raised in the original
ask stands: they don't have to be "Lifestyle" nodes specifically — the exact
set should be decided when this is implemented, driven by what the three
JSON files actually require, not fixed in advance here.

## Keeping new auxiliary nodes invisible to the existing UI

The requirement is: only the node/edge types the D3/Streamlit viewer already
renders today stay visible there; any new parameter-carrier nodes should not
show up in the bundled map by default. Two ways to get this, in order of
preference:

1. **Type-allowlist at export time.** `exporter.to_graph_data_js` (the
   function that regenerates `map/graph-data.js` for the viewer) filters
   nodes/edges to the set of types already rendered today (`Carcinogen`,
   `Enzyme`, `Gene`, `Metabolite`, `DNA_Adduct`, `Pathway`, plus whichever
   edge types are already used). Any new auxiliary type is excluded
   automatically, with no per-node flag to maintain. This also lines up with
   the JSON/JS split from the graph-data migration commit: the full,
   parameter-enriched graph lives in `map/graph-data.json` and is what
   `flux_engine`/`exposure_engine`/`interaction_engine` query; the UI keeps
   consuming a filtered `graph-data.js` derived from it.
2. **Explicit visibility flag**, e.g. `Node.role: Literal["domain",
   "auxiliary"] = "domain"`, checked by `filter_knowledge_graph`
   (`graph_filters.py`) alongside the existing `GraphVisibility` /
   `match_status` filtering. More flexible (lets a future UI opt into
   showing auxiliary nodes), but is new schema surface the allowlist
   approach avoids.

Recommendation: start with (1) since it requires no schema change and
reuses the export path that already exists; revisit (2) only if a future UI
needs to selectively surface specific auxiliary nodes.

## Ingestion mechanism (sketch, not implemented here)

A single ingestion step — likely a new function alongside
`reference_data.build_reference_graph()`/`build_reference_engine()`, or a
dedicated `graph_parameters.py` — would read the three JSON files once,
build the corresponding `Node`/`Edge` objects, and `engine.merge(...)` them
into the reference graph at build time (reusing the validated merge path
`GraphEngine` already has). `flux_engine`, `exposure_engine`, and
`interaction_engine` would then be refactored (separately, in their own
discrete commits) to call `engine.get_node(...)`/`engine.get_edge(...)`
(added in the prior commit) instead of their private `_load_*_params()`
JSON caches.

## Open questions to resolve before implementation

- **Class vs. specific-carcinogen granularity.** `kinetic_parameters.json`
  keys parameters by `carcinogen_classes.<class>` (e.g. `PAH`, `Aflatoxin`)
  with a single `index_carcinogen`, while the graph has 60 specific
  carcinogen nodes each tagged with a `group` attribute (e.g. `BaP` →
  `group: "PAH"`). The `group` values don't necessarily match the JSON's
  class keys one-to-one (e.g. graph uses `Mycotoxin`/`Alkylating` naming in
  some places) — a resolver/mapping table will be needed either way, whether
  parameters land directly on each specific carcinogen's edges or on a
  shared class-level grouping node that specific carcinogens point to.
- **Multi-edge disambiguation.** An enzyme can have distinct Km/Vmax for
  multiple carcinogens and distinct Ki for multiple competing substrates.
  The `MultiDiGraph` already supports parallel edges, and the new
  `get_edge_keys`/`edge_key` parameter (previous commit) gives a way to
  disambiguate them — but the ingestion step needs a stable, deterministic
  scheme for assigning edge keys so lookups are reproducible across rebuilds.
- **Backward compatibility during migration.** Per the discrete-commit
  constraint, `flux_engine`/`exposure_engine`/`interaction_engine` should
  probably keep reading their JSON caches until the graph-backed path is
  validated to produce identical results (e.g. a temporary dual-path with a
  test asserting parity), then have the JSON path removed in a later commit.
- **GSH and other small-molecule nodes.** A `GSH` node (plus conjugates
  like `AFB1_GSH`, `BPDE_GSH`) already exists in the bundled graph, so
  `gsh_depletion` parameters likely attach to existing nodes/edges without a
  new node type — this should be confirmed against the full node type of
  `GSH` and its existing edges before implementation.
- **Class-key naming is inconsistent today, beyond granularity.**
  Spot-checking confirms this isn't just an occasional mismatch:
  `kinetic_parameters.json` uses `Aflatoxin`, `ChlorinatedSolvent`,
  `HeavyMetal`, `Benzene`, `NDMA`; the graph's `group` attribute uses
  `Mycotoxin`, `Chlorinated_Solvent`, `Heavy_Metal`, and has no `Benzene`
  or `NDMA` group at all (those are presumably specific carcinogens filed
  under a different group). A resolver/mapping table is a prerequisite for
  any class-level ingestion approach, not a nice-to-have.
