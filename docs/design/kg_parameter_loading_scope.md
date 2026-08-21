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

## Addendum: `interaction_parameters.json` substrate-key coverage (verified)

The proposed mapping table above says `interaction_parameters.json`
substrates can ride on existing edges "using `Edge.carcinogen` as the
competing substrate context." That assumed the substrates already have
graph representation. A full inventory of `competitive_inhibition` shows
that assumption holds for only a small minority of the data:

- `competitive_inhibition` has 74 `(enzyme, substrate)` entries across 12
  enzymes, using 54 unique substrate keys.
- **8 of 54** substrate keys exact-match an existing `Carcinogen` node:
  `AFB1, BaP, DMBA, MeIQx, NDEA, NDMA, NNK, PhIP`.
- **5 of 54** are pure case/naming-convention mismatches, safe to
  auto-resolve via a small alias table: `benzene`→`Benzene`,
  `ethanol`→`Ethanol`, `vinyl_chloride`→`VinylChloride`,
  `cyclophosphamide`→`Cyclophosphamide`, `testosterone`→`Testosterone`.
- **41 of 54** have no corresponding node under any casing. These are not
  one problem:
  - Real IARC-relevant carcinogens/toxicants simply absent as graph nodes
    today: `naphthalene, styrene, chloroform, trichloroethylene (graph's
    "TCE"), nicotine, chrysene, benzo_a_anthracene, dibenz_ah_anthracene,
    1_nitropyrene, 3_methylindole, 4_aminobiphenyl (graph's "4ABP"),
    6_aminochrysene, AalphaC, IQ, sterigmatocystin, methoxsalen, NNAL,
    cotinine`. A content gap, not an aliasing problem — **deferred** to a
    separate content-curation effort, not in scope for the ingestion
    mechanism itself.
  - Non-carcinogenic CYP "probe substrates" from the pharmacology
    literature, used only to characterize enzyme kinetics, not carcinogens
    by any classification: `caffeine, midazolam, dextromethorphan,
    debrisoquine, S_warfarin, S_mephenytoin, tolbutamide, coumarin,
    bufuralol, bupropion, nifedipine, erythromycin, diclofenac,
    cyclosporine, omeprazole, p_nitrophenol, 7_ethoxyresorufin,
    phenacetin, theophylline, acetaminophen, resveratrol`. These do not
    obviously belong under `NodeType.CARCINOGEN`. **Deferred**: no new
    node type introduced yet; the ingestion mechanism should skip these
    (leave them JSON-only) rather than force-fit them, and
    `interaction_engine` keeps a narrow non-graph fallback for exactly
    these keys until/unless a `Xenobiotic`-style node type is deliberately
    added.
  - Two ambiguous cases needing SME confirmation, not string matching:
    `estradiol_2_OH`/`estradiol_4_OH` plausibly correspond to the existing
    Metabolite nodes `2OHE2`/`HydroxyE2` (the graph already has an
    `E2`→`2OHE2` `ACTIVATES` edge), but the naming conventions are
    unrelated enough that no mechanical rule should auto-resolve this.
    **Left unresolved** pending confirmation.
- The identical substrate-coverage problem exists in the sibling block
  `phase2_conjugation`: its outer enzyme keys match `Enzyme` nodes exactly,
  but its own nested `substrates` keys have the same gap (e.g.
  `UGT1A1.substrates.bilirubin` — an endogenous compound, no graph node).
- **Edge coverage gap, even for matched carcinogens.** Checking whether an
  edge already exists tagged with the right `Edge.carcinogen` value for
  each of the 21 `(enzyme, carcinogen)` pairs among the 8 exact matches:
  only **10 of 21** have one. The other 11 — including `CYP1A1`/`PhIP` and
  `NNK` against five separate enzymes — would need a brand-new edge
  created before any Km/Vmax could attach to it. **Scope for the first
  ingestion commit: only the 10 pairs with an existing edge**; the other
  11 are reported as warnings (mirroring the tissue-expression warning
  pattern) and deferred to a follow-up commit once product-node identity
  (see below) is settled.
- **`product` values are worse still.** Across all 74 substrate entries
  there are 57 unique `product` values (what the reaction yields); only 1
  (`BPDE`) matches an existing `Metabolite` node. For the first ingestion
  pass, `product`/`product_carcinogenic` stay as plain string/bool fields
  inside `Edge.kinetics` rather than becoming node references — minting
  ~56 new `Metabolite` nodes is a separate, larger content decision.
- **A third naming layer sits upstream of all of this.**
  `interaction_engine.competitive_inhibition_flux()` takes an arbitrary
  `substrates: dict[str, float]`, but that dict is normally built by
  `_build_competitive_substrates()`, which hardcodes a small mapping from
  exposure categories (`"PAH"`, `"HCA"`) to only a handful of the 54
  substrate names. Most `competitive_inhibition` substrate entries are
  never reached through the normal exposure-profile pathway today — only
  through direct calls with an explicit substrate name (e.g. validation
  code). Making `interaction_engine` graph-first eventually requires this
  bridge to become graph-driven too (`CARCINOGEN_ENZYME_MAP` is a fourth,
  separate hardcoded vocabulary covering the same ground) — out of scope
  for the initial ingestion commit, called out here so it isn't lost.

## Addendum 2: correction — non-carcinogen substrates are a designed category, not a gap

The "41 of 54 unmatched substrates" framing above is too pessimistic for a
subset of that group. `interaction_schema.ReactionRole` already defines
`PROBE_ONLY` alongside `BIOACTIVATION`/`DETOXIFICATION`/`CLEARANCE`, and
`endpoint_toxic_flux._interpretation_channel()` already routes
`PROBE_ONLY` (and `UNKNOWN`, `DUAL_ROLE`) substrates to a neutral outcome
that contributes nothing to `activation_burden_ratio` or
`detox_failure_ratio`. Separately, `kinetic_resolver.get_ki(enzyme,
inhibitor, target_substrate, ...)` treats "inhibitor" and "target
substrate" as distinct roles, and `product_carcinogenic` on every JSON
substrate entry already zeroes `activated_product_flux` when the reaction
product isn't carcinogenic. **The pipeline already has a mechanism for
correctly excluding non-carcinogenic substrates from toxic burden that
has nothing to do with graph-node identity.** A substrate does not need a
`Carcinogen` node to be handled correctly.

This resolves three things that were previously open questions:

1. **Design principle: node identity and reaction-role gating are
   decoupled.** Whether a substrate gets a graph node is a
   provenance/queryability decision, not a correctness requirement. The
   correctness requirement is that every substrate have an explicit
   `reaction_role` (not fall through to `UNKNOWN` by omission). This
   replaces the earlier framing that treated "no matching node" as
   inherently a gap to close.

2. **Concrete follow-up work item (schema-data only, not this ingestion
   commit): tag the untagged substrates explicitly.** Checking
   `reaction_role_semantics.get_reaction_role_sme_records()` today shows
   curated records only for `benzene`, `NDMA`, `vinyl_chloride`, `HCA`
   (x2 enzymes), and one pending `trichloroethylene` candidate — **zero**
   of the 21 CYP-phenotyping probe substrates (`caffeine, midazolam,
   dextromethorphan, debrisoquine, S_warfarin, S_mephenytoin, tolbutamide,
   coumarin, bufuralol, bupropion, nifedipine, erythromycin, diclofenac,
   cyclosporine, omeprazole, p_nitrophenol, 7_ethoxyresorufin, phenacetin,
   theophylline, acetaminophen, resveratrol`) and none of the endogenous
   substrates (`bilirubin`, `testosterone`) have any `ReactionRoleAnnotation`
   record. They currently fall through to `ReactionRole.UNKNOWN` by
   default, which is neutral by accident rather than by explicit design.
   Proposed follow-up (separate commit, SME sign-off needed before
   marking any record `SMEReviewStatus.CURATED`): add explicit
   `ReactionRoleAnnotation` records with `role=PROBE_ONLY` for the 21
   pharmacology probe substrates, and `role=CLEARANCE` for `bilirubin`
   (matches its existing JSON note: "Endogenous UGT1A1 substrate ...
   Gilbert syndrome") and `testosterone` (matches its role as an
   androgen-clearance-pathway phenotype marker). This is a data/registry
   change to `reaction_role_semantics.py`, not a knowledge-graph change,
   and is independent of when interaction_parameters.json ingestion
   happens.

3. **Resolved: `Edge.carcinogen` stays carcinogen-only; probe/endogenous
   kinetics stay off the edge model.** Generalizing `Edge.carcinogen` to a
   broader `competing_substrate` field (so probe substrates could also
   live in `Edge.kinetics`) was considered and rejected for the first
   ingestion pass: probe/endogenous substrates never need graph
   traversal — they are leaf inputs to a calculation, not entities any
   query needs to reach relationships from. Their Km/Vmax/product data
   stays sourced from a lightweight name + `reaction_role` keyed lookup
   outside the graph (i.e., still JSON/registry-backed, not a graph
   query). `Edge.kinetics` and `Edge.carcinogen` are reserved exclusively
   for the 10 (enzyme, carcinogen) pairs that already have a real edge
   and a real `Carcinogen` node backing them, keeping the edge schema's
   meaning unambiguous.

## Addendum 3: substrate-node identity, visibility, and provider wiring

### Identity: no new edge type, extend metadata instead

`grounding.py` already has an established mechanism for "is this the same
real-world entity as an existing node" — `Node.match_status`
(`CANONICAL`/`ALIAS`/`CUSTOM`/`UNMATCHED`) plus `canonical_id`/
`canonical_label` fields on the node itself. (Note: `match_status`'s exact
semantics are narrower/more specific than a general equivalence flag —
this doc treats it only as "the existing term-matching mechanism
grounding.py already uses," without relying on a precise definition
beyond that.) A new `IS_EQUAL_TO` edge type was considered and rejected:
it would create a second, competing identity mechanism alongside the
existing node-level one, and every consumer (`kinetic_resolver`,
`interaction_engine`, `endpoint_toxic_flux`, `mechanism_attribution`)
would need new hop-resolution logic before any lookup.

Instead, substrates split into three identity buckets:

- **Bucket A (13 of 54): substrate already equals an existing
  `Carcinogen` node** (8 exact + 5 case/naming aliases from Addendum 1).
  No new node. Per-mechanism data (Km/Vmax/product/`reaction_role`) is
  enzyme- and tissue-specific, so it belongs on the **edge**
  (`Edge.kinetics`, the same pattern already used for `Edge.carcinogen`
  disambiguation), not on the node.
- **Bucket B: no carcinogen counterpart exists at all** (21 CYP
  phenotyping probe substrates + `bilirubin`). New nodes under a new
  `NodeType.SUBSTRATE`, so the type itself signals "queryable kinetics
  entity, not a toxicological classification."
- **Bucket C: real carcinogens/toxicants with no node today**
  (`naphthalene, styrene, chloroform, trichloroethylene, nicotine,
  chrysene, benzo_a_anthracene, dibenz_ah_anthracene, 1_nitropyrene,
  3_methylindole, 4_aminobiphenyl, 6_aminochrysene, AalphaC, IQ,
  sterigmatocystin, methoxsalen, NNAL, cotinine`). Same `NodeType.SUBSTRATE`
  as an interim placeholder, using `Node.generate_id()` the same way a
  `Carcinogen` node would, so promoting one to `NodeType.CARCINOGEN`
  later is a type/metadata edit on the same node id, not a re-key.

### Visibility: filter by NodeType at the render boundary, not by match_status

The existing `GraphVisibility` filter (`graph_filters.filter_knowledge_graph`,
used by the D3 viewer tab) filters on `match_status`, which is a different
concern from "should this render at all." Reusing it to hide
`NodeType.SUBSTRATE` nodes would be fragile — a validated `SUBSTRATE`
node would still leak into `VALIDATED_ONLY` view. Proposed instead: an
explicit `NodeType` exclusion at the two rendering entry points
(`ui_map_viewer.py`'s bundled D3 export and
`exporter.to_interactive_html_string`), symmetric with the existing
`filter_knowledge_graph` pattern but keyed on type. `GraphEngine` query
methods (`get_data`/`get_node`) stay unaffected — the full graph remains
queryable; only the two rendering paths filter.

### Provider wiring: parameter_provider.py is the only required new implementation

`InteractionParameterProvider` (in `parameter_provider.py`) is already an
abstract base class; `JSONInteractionParameterProvider` is the only
concrete implementation today. A new `GraphInteractionParameterProvider(engine)`
implementing the same ABC — using `engine.get_data(enzyme, substrate,
key="kinetics")` for parameters and enumeration over enzyme-node edges for
listing — is a peer implementation, not a rewrite:

- `kinetic_resolver.py`: no logic changes. `KineticParameterResolver` only
  calls provider methods; the one touch point is that its constructor's
  type hint is currently the *concrete* `JSONInteractionParameterProvider`
  class rather than the ABC, so it needs widening — a one-line signature
  change, called out for the audit trail even though it's trivial.
- `interaction_equations.py`: no changes. It only operates on already-
  resolved numeric values, never touches JSON or the graph.
- `interaction_schema.py`: no changes. Its dataclasses are already
  provider-agnostic; a graph-backed provider just populates the same
  shapes.
- **Latent inconsistency surfaced (not created) by this work:**
  `reaction_role` is dual-sourced today. Verified 0 of 74
  `competitive_inhibition` JSON entries carry a `reaction_role` key, so
  `JSONInteractionParameterProvider._build_competitive_interaction`
  always defaults `CompetitiveInteraction.reaction_role` to `UNKNOWN`.
  The actual curated roles (the benzene/NDMA/vinyl_chloride/HCA/TCE
  records from Addendum 2) live in a separate registry,
  `reaction_role_semantics.get_default_reaction_role_registry()`, queried
  directly by `interaction_engine.py`/`endpoint_toxic_flux.py`, bypassing
  the provider layer entirely. **Open question to resolve before
  implementation:** should the graph become the single source of truth
  for `reaction_role` (provider surfaces it from a graph edge/node
  attribute), or does `reaction_role_semantics.py` stay as a Python-code
  overlay the graph-backed provider still consults?
- **Small additive engine gap:** `GraphEngine` has no method to enumerate
  all edges of a given type from a node (`neighbors()` returns bare node
  ids only, no type filter). A new `edges_from(node_id, edge_type=None)`
  getter is needed, in the same spirit as the existing `get_node`/
  `get_edge`/`get_data` getters added earlier on this branch.

## Implementation commit plan (scoping only — none of this has been started)

Proposed discrete, independently auditable commits, in dependency order:

1. **Add `NodeType.SUBSTRATE`** to `models.py` + tests. Schema-only, no
   data population, no behavior change — mirrors the existing
   `NodeType.TISSUE` precedent (added but zero-instance until used).
2. **Add `GraphEngine.edges_from(node_id, edge_type=None)`** + tests.
   Mirrors the symmetric-getters precedent from commit `c489766`. No
   callers yet.
3. **Add type-based rendering exclusion** for `NodeType.SUBSTRATE` at the
   two render entry points + tests confirming the bundled reference
   graph's render output is unchanged today (it has zero `SUBSTRATE`
   nodes yet, so this should be a no-op until commit 5).
4. **Populate Bucket A**: extend the 10 existing (enzyme, carcinogen)
   edges with `Edge.kinetics` from `competitive_inhibition`, resolving
   the 5 case/naming aliases from Addendum 1 to the right edge endpoints;
   create the 11 missing (enzyme, carcinogen) edges for the remaining
   matched pairs, with warnings for anything still unresolved (mirroring
   the tissue-expression warning pattern from commit `d86ad99`).
5. **Create `NodeType.SUBSTRATE` nodes for Buckets B and C** (21 probe
   substrates + `bilirubin`, plus the real missing carcinogens as
   placeholders) and their Enzyme→Substrate edges with kinetics. The
   bulk data-population commit, kept separate from provider wiring.
6. **Add explicit `ReactionRoleAnnotation` records** (or the graph-
   sourced equivalent, depending on how the open `reaction_role`
   single-source-of-truth question above is resolved) for the
   newly-tagged substrates — SME-reviewed, per the follow-up scoped in
   Addendum 2.
7. **Add `GraphInteractionParameterProvider`** implementing
   `InteractionParameterProvider`, backed by `GraphEngine` queries over
   commits 4–6's data; widen `KineticParameterResolver`'s constructor
   type hint to the ABC. Additive and independently tested — no
   call-site swap yet.
8. **Cut over `interaction_engine.py`** to construct its
   `KineticParameterResolver` with `GraphInteractionParameterProvider`
   instead of `JSONInteractionParameterProvider`. Parity tests confirming
   identical output to the JSON-backed path are required before this
   commit, since it's the actual behavior-affecting cutover the whole
   plan has been building toward.

## Addendum 4: `SUBSTRATE_OF` edge type + deferred direct-edge cleanup (TODO: 18 unclear cases)

As part of fixing the overloaded `ACTIVATES`/`DETOXIFIES` edge type (direct
`Carcinogen/Metabolite -> Carcinogen/Metabolite` edges that duplicated the
information already carried by a parallel `Enzyme -> Metabolite` edge), a
full audit of all 68 such direct edges was run and categorized:

- **33 `REDUNDANT_PARALLEL_ENZYME_EDGE`** and **6 `ENZYME_NAMED_IN_PROSE_NOT_MODELED`**
  (39 total) were restructured: the direct edge was removed, a new
  `EdgeType.SUBSTRATE_OF` edge (`Carcinogen/Metabolite -> Enzyme`) was added
  in its place, and (for the 6 prose-only cases) the missing
  `Enzyme -> Metabolite` `ACTIVATES` edge was added so the chain reads
  `Carcinogen/Metabolite -[SUBSTRATE_OF]-> Enzyme -[ACTIVATES]-> Metabolite`
  end to end. This lets a query on an enzyme's neighbors show every
  carcinogen/metabolite competing for that enzyme (motivating use case:
  reasoning about competitive inhibition).
- **11 `LIKELY_LEGITIMATE_NON_ENZYMATIC`** were left as direct edges
  (spontaneous, non-enzymatic transformations — no enzyme to restructure
  through).
- **18 `UNCLEAR_NEEDS_SME_REVIEW` were deliberately deferred, not yet fixed.**
  These still have the old direct-edge shape and need a domain expert to
  confirm whether they are (a) genuinely spontaneous/non-enzymatic (→ leave
  as-is, reclassify as `LIKELY_LEGITIMATE_NON_ENZYMATIC`), (b) enzyme-mediated
  but with the catalyzing enzyme not yet identified in the literature review
  done for this pass (→ restructure via `SUBSTRATE_OF` like the 39 above once
  the enzyme is known), or (c) something else entirely (e.g. two distinct
  sequential steps mis-modeled as one edge). The 18 are:

  | Source | Target | Edge type |
  |---|---|---|
  | LeadInorganicCompounds | GSH | ACTIVATES |
  | HydroxyE2 | E2_quinone | ACTIVATES |
  | Benzene_oxide | HQ | ACTIVATES |
  | HQ | Benzoquinone | ACTIVATES |
  | Cd | Cd_MT | ACTIVATES |
  | Cd | ROS_metal | ACTIVATES |
  | CrVI | Cr_III | ACTIVATES |
  | CrVI | Cr_V | ACTIVATES |
  | CrVI | Asc_Cr_III | ACTIVATES |
  | NickelCompounds | ROS_metal | ACTIVATES |
  | MMA_III | ROS_metal | ACTIVATES |
  | Cr_V | ROS_metal | ACTIVATES |
  | Chloral_hydrate | TCA | ACTIVATES |
  | Radon | Radon_decay_products | ACTIVATES |
  | EstrogenProgestogenTherapy | E2 | ACTIVATES |
  | Naphthalene_1_2_oxide | Naphthoquinone | ACTIVATES |
  | MethylmercuryCompounds | MeHg_GSH | ACTIVATES |
  | GSH | MeHg_GSH | ACTIVATES |

  Full categorization detail (why each of the 68 landed where it did) is in
  `audit/activates_detoxifies_direct_edge_audit.csv` in the working tree
  this pass was done in (not yet committed to the repo as of this addendum —
  should be moved under `docs/` or `data/audits/` if kept long-term).

**Companion fixes required by the restructuring above** (also applied this
pass, since they'd otherwise silently regress existing functionality):

- `ExposoGraph/tissue_subgraphs.py`: two edge-type allow-lists (tissue-subgraph
  node inclusion logic) needed `SUBSTRATE_OF` added, or carcinogens connected
  to a tissue-relevant enzyme *only* via the now-removed direct edge would
  silently drop out of tissue subgraphs.
- `ExposoGraph/graph_analysis.py`: `metabolism_chain`'s `_METABOLISM_EDGE_TYPES`
  allow-list needed `SUBSTRATE_OF` added for the same reason — otherwise the
  first hop of the metabolism chain for all 39 restructured pathways would be
  invisible to that function.
- `ExposoGraph/engine.py`'s `paths_from_carcinogen`/`paths_to_carcinogen`
  needed no change — their default (`edge_types=None`) already walks every
  edge type, `SUBSTRATE_OF` included.
- `ExposoGraph/unified_api.py`'s kinetic-parameter attachment (keyed on
  `{EdgeType.ACTIVATES, EdgeType.DETOXIFIES}`) needed no change — `SUBSTRATE_OF`
  correctly carries no reaction role / kinetics of its own (see
  `ExposoGraph/reaction_role_rules.py`'s docstring).
- Edge-color maps in `exporter.py` / `_app_shared.py` (`ui_data.py`,
  `ui_preview.py` consumers) all use `.get(edge_type, <default color>)`, so
  `SUBSTRATE_OF` renders in the generic fallback color rather than crashing.
  Adding a dedicated color for it is a cosmetic follow-up, not yet done.

## Addendum 5: directional-schema redesign — `ACTIVATES`/`DETOXIFIES`/`TRANSPORTS`/`REPAIRS` split into role-specific directional edge types

Following on from Addendum 4, `ACTIVATES`, `DETOXIFIES`, `TRANSPORTS`, and
`REPAIRS` were themselves overloaded: each mixed together several distinct
reaction roles (Enzyme→product formation, receptor agonism, spontaneous
non-enzymatic transformation, and biologically ambiguous/undetermined
mechanism) under one label with an arbitrary direction convention. This pass
splits every edge of these four legacy types into a new, direction- and
role-explicit edge type, based on a fresh audit of all 327 nodes / 557 edges
in `graph-data.json` (edge count unchanged — this is a pure retyping pass,
not a data-scope change).

### New edge types (added to `EdgeType` in `models.py`)

| New type | Meaning | Direction |
|---|---|---|
| `PRODUCES` | Enzyme → product it forms | Enzyme → Metabolite/Substrate |
| `DETOXIFIED_BY` | Compound is cleared by an enzyme (self-referential — target node IS the substrate) | Compound → Enzyme |
| `AGONIZES` | Carcinogen/Metabolite activates a nuclear receptor | Carcinogen/Metabolite → Receptor |
| `TRANSPORTED_BY` | Compound moved by a transporter | Compound → Transporter |
| `REPAIRED_BY` | DNA adduct repaired by an enzyme | DNA_Adduct → Enzyme |
| `TRANSFORMS_SPONTANEOUSLY` | Non-enzymatic, spontaneous chemical conversion | Compound → Compound |
| `MECHANISM_UNCLEAR` | Compound→compound conversion whose mechanism is not yet resolved (superset of Addendum 4's deferred-18) | Compound → Compound |

Also added `NodeType.RECEPTOR = "Receptor"`; 8 nodes previously typed
`Enzyme` (`AHR, CAR, PXR, PPARA, ESR1, AR, RyR, HLA_DPB1`) were reclassified
to `Receptor`, since they are the sole targets of the new `AGONIZES` edges
and were never actually enzymes.

**The old `ACTIVATES`, `DETOXIFIES`, `TRANSPORTS`, `REPAIRS` enum members
were deliberately kept, not removed**, per explicit decision (see below) —
`graph-data.json` itself no longer contains any edge of these four types
after this migration, but the enum members remain load-bearing for code
outside this migration's scope.

### Reclassification breakdown

- **`ACTIVATES` (151 edges)** →
  - 104 Enzyme→X edges → `PRODUCES`
  - 18 X→Receptor edges (15 Carcinogen + 3 Metabolite, to the 8 receptor
    nodes) → `AGONIZES`
  - 29 compound↔compound edges → 11 → `TRANSFORMS_SPONTANEOUSLY`
    (non-enzymatic, mechanism understood), 18 → `MECHANISM_UNCLEAR`
    (verbatim the same 18 pairs deferred as `UNCLEAR_NEEDS_SME_REVIEW` in
    Addendum 4 — that backlog is still authoritative and unresolved, just
    retyped)
- **`DETOXIFIES` (87 edges, 100% Enzyme-sourced)**, split by whether the
  target is self-referential (same substance being cleared) or a genuinely
  distinct downstream product with its own, separately-modeled precursor
  node:
  - **72 edges → `DETOXIFIED_BY`** (`Compound -[DETOXIFIED_BY]-> Enzyme`):
    52 fully self-referential (36 Enzyme→Substrate, 16 Enzyme→Carcinogen)
    plus 20 of the 35 Enzyme→Metabolite edges confirmed self-referential
    after individual biochemical verification (not name-heuristics —
    e.g. `ALDH2 -> Acetaldehyde_int` looked like a distinct-product case at
    first pass but `Acetaldehyde_int` already has its own formation edges
    from `ADH1B`/`ADH1C`, so it is in fact self-referential).
  - **15 edges → `PRODUCES`** (kept as `Enzyme -[PRODUCES]-> Product`,
    direction unchanged from the reclassified-`ACTIVATES` convention above):
    genuine, distinct downstream products whose true precursor node exists
    in the graph but currently has **no edge to the producing enzyme** —
    a backlog gap, not fixed in this pass. The 15, with their missing
    precursor→enzyme entry edge:

    | Enzyme | Product (kept) | Missing precursor edge |
    |---|---|---|
    | CYP3A4 | HydroxyTestosterone | `Testosterone -> CYP3A4` |
    | CYP3A5 | HydroxyTestosterone | `Testosterone -> CYP3A5` |
    | GSTM1 | BPDE_GSH | `BPDE -> GSTM1` |
    | GSTP1 | BPDE_GSH | `BPDE -> GSTP1` |
    | GSTM1 | AFB1_GSH | `AFB1_epoxide -> GSTM1` |
    | GSTT1 | AFB1_GSH | `AFB1_epoxide -> GSTT1` |
    | NAT2 | PhIP_NAc | `PhIP -> NAT2` |
    | UGT1A1 | PhIP_gluc | `PhIP -> UGT1A1` |
    | COMT | E2_methyl | `HydroxyE2 -> COMT` (precursor identity tentative) |
    | UGT2B17 | Testosterone_gluc | `Testosterone -> UGT2B17` |
    | UGT2B17 | DHT_gluc | `DHT -> UGT2B17` |
    | UGT2B15 | DHT_gluc | `DHT -> UGT2B15` |
    | AKR1C2 | 3aAdiol | `DHT -> AKR1C2` |
    | ALDH2 | Formate | `Formaldehyde -> ALDH2` (node existence unchecked) |
    | ADH5 | Formate | `Formaldehyde -> ADH5` (node existence unchecked) |

    **This 15-row list supersedes an earlier, incorrect first-pass
    Category-B list** that had been built from name-matching heuristics
    alone (missed `BPDE_GSH`/`AFB1_GSH`, and would have mis-flipped them to
    `DETOXIFIED_BY`, which is biochemically wrong — a glutathione conjugate
    is a product, not the substrate being cleared). Fixed before the
    migration script was applied to `graph-data.json`; the numbers above are
    already the corrected, applied state.
- **`TRANSPORTS` (13 edges, 100% Enzyme-sourced)** → `TRANSPORTED_BY`.
- **`REPAIRS` (44 edges, 100% Enzyme→DNA_Adduct)** → `REPAIRED_BY`.
- `FORMS_ADDUCT` (76), `PATHWAY` (120), `SUBSTRATE_OF` (56), `INHIBITS` (7),
  `INDUCES` (3) — unchanged.

Post-migration edge-type counts in `graph-data.json` (557 total, unchanged):
`PATHWAY` 120, `PRODUCES` 119, `FORMS_ADDUCT` 76, `DETOXIFIED_BY` 72,
`SUBSTRATE_OF` 56, `REPAIRED_BY` 44, `AGONIZES` 18, `MECHANISM_UNCLEAR` 18,
`TRANSPORTED_BY` 13, `TRANSFORMS_SPONTANEOUSLY` 11, `INHIBITS` 7,
`INDUCES` 3.

The `canonical_predicate` field (present on 198/557 edges, always equal to
`type` when present) was updated in lockstep with `type` on every edge this
migration touched. The `carcinogen` field (informational — traces which
top-level carcinogen an edge relates to) was left untouched.

### Explicit decision: `figure_architecture.py`, `unified_api.py`, `seeder.py` are NOT updated

A blast-radius check found the legacy edge types are load-bearing beyond
`tissue_subgraphs.py`/`graph_analysis.py` (patched below):

- `figure_architecture.py` — dedicated figure-generation module with
  hardcoded per-edge-type colors/linestyles and hardcoded edge counts
  (`ACTIVATES: 57, DETOXIFIES: 31, TRANSPORTS: 9, REPAIRS: 19`) for what
  looks like a validated/published architecture diagram.
- `unified_api.py` — kinetic parameter attachment gated on
  `edge.type in {EdgeType.ACTIVATES, EdgeType.DETOXIFIES}`.
- `seeder.py` — default edge-type inference logic returns the legacy types
  for future data loads.

**Decision (explicit, user-confirmed): migrate `graph-data.json` only, leave
these three files untouched.** This is why the legacy `EdgeType` enum
members were kept rather than removed — deleting them would break these
three files' imports/references even though `graph-data.json` no longer
contains any edge instantiating them. `reaction_role_rules.py` is separate,
uncommitted, unrelated pending work and was not touched in this pass either.

### Companion fixes required (regression prevention, same pattern as Addendum 4)

- `ExposoGraph/tissue_subgraphs.py`: both edge-type allow-lists (tissue node
  inclusion, ~line 1147 and ~line 1247) extended to include `PRODUCES`,
  `DETOXIFIED_BY`, `TRANSPORTED_BY`, `REPAIRED_BY`, `AGONIZES`,
  `TRANSFORMS_SPONTANEOUSLY`, `MECHANISM_UNCLEAR` alongside the existing
  legacy/SUBSTRATE_OF entries — otherwise nodes connected to a tissue-
  relevant enzyme only via a newly-retyped edge would silently drop out of
  tissue subgraphs.
- `ExposoGraph/graph_analysis.py`: `_METABOLISM_EDGE_TYPES` extended with
  the same 7 new types, for the same reason (`metabolism_chain` would
  otherwise treat every retyped edge as invisible).
- `ExposoGraph/engine.py`'s `paths_from_carcinogen`/`paths_to_carcinogen`
  needed no change (default `edge_types=None` walks every type).

### Validation performed

- Migration script (`audit/migrate_directional_schema.py`) dry-run then
  `--apply`, editing `graph-data.json` via raw `json.load`/`json.dump`
  (no pydantic `model_dump()` reserialization, per standing convention).
  Confirmed zero duplicate `(source, target, type)` tuples introduced, and
  total edge count unchanged (557 before/after).
- Reloaded via `build_reference_engine()` (from
  `ExposoGraph/reference_data.py`) — graph loads and validates cleanly:
  327 nodes / 557 edges, node-type counts `Metabolite 80, Enzyme 69,
  Carcinogen 66, Substrate 49, DNA_Adduct 39, Pathway 16, Receptor 8`.
- `tools/graph_role_consistency_check.py` (standalone QA script, not
  pytest) run post-migration: 4 PASS, 0 WARN, 0 FAIL.
- `paths_from_carcinogen` sanity-checked for `BaP` (3,614 maximal paths),
  `PhIP` (1,201 paths), and `Cyclophosphamide` (3,537 paths) — all traverse
  the new edge types correctly (e.g. `Cyclophosphamide
  -[TRANSFORMS_SPONTANEOUSLY]-> Acrolein_CP -[FORMS_ADDUCT]-> Acr_dG
  -[REPAIRED_BY]-> XRCC1`).

### Still outstanding (backlog, not fixed in this pass)

1. The 18 `MECHANISM_UNCLEAR` edges — same list as Addendum 4's
   deferred-18, now retyped but still unresolved; needs SME review to
   determine enzymatic-vs-spontaneous mechanism.
2. The 15-row Category B precursor-gap table above — each needs a new
   precursor→enzyme entry edge added once curated/sourced. Note
   `Formaldehyde` node existence for the `ALDH2`/`ADH5` → `Formate` rows was
   not checked in this pass.
3. `figure_architecture.py` (hardcoded legacy edge-type counts/colors),
   `unified_api.py` (kinetic-parameter gating keyed on legacy types), and
   `seeder.py` (legacy-type inference for future loads) still reference the
   pre-migration schema and were explicitly left unmodified — any future
   work that touches these should be aware `graph-data.json` no longer has
   any edge of the legacy types they check for.
4. Dedicated colors/linestyles for the 7 new edge types in exporter/UI code
   (`exporter.py`, `_app_shared.py`, `ui_data.py`, `ui_preview.py`) are not
   yet added — same fallback-color situation noted for `SUBSTRATE_OF` in
   Addendum 4.

## Addendum 6: filled in missing AHR/CAR/PXR `INDUCES` edges

Addendum 5 introduced `AGONIZES` (ligand→receptor binding) as distinct from
`INDUCES` (receptor→enzyme transcriptional induction), but a review found
several `AGONIZES` edges whose own evidence text already named a downstream
induced enzyme that had no corresponding `INDUCES` edge. Six were added,
each carrying the `carcinogen` field of the specific `AGONIZES` edge whose
evidence text motivated it (all target enzyme/transporter nodes already
existed in the graph — no new nodes added):

| New edge | Sourced from (`carcinogen` field) | Citation |
|---|---|---|
| `AHR -[INDUCES]-> CYP1A2` | HCB | IARC Monograph 79, 2001 |
| `CAR -[INDUCES]-> CYP2B6` | PCB_non_dioxin | IARC Monograph 107 |
| `CAR -[INDUCES]-> CYP3A4` | PCB_non_dioxin | IARC Monograph 107 |
| `PXR -[INDUCES]-> CYP3A4` | PCB_non_dioxin | IARC Monograph 107 |
| `PXR -[INDUCES]-> ABCB1` | PCB_non_dioxin | IARC Monograph 107 |
| `PXR -[INDUCES]-> ABCC2` | PCB_non_dioxin | IARC Monograph 107 |

`PPARA` was reviewed but not extended: its `AGONIZES` edges (`PFOA`, `PFOS`,
`PCE`) only cite vague "peroxisome proliferation" language without naming a
specific enzyme, and the textbook PPARα target (`ACOX1`) does not yet exist
as a node in the graph — adding that `INDUCES` edge would require adding a
new node first, which was out of scope for this pass.

Total edge count: 557 → 563. Applied via
`audit/add_induces_edges.py` (raw `json.load`/`json.dump`, no pydantic
`model_dump`). Re-validated: `build_reference_engine()` loads cleanly
(327 nodes / 563 edges), `tools/graph_role_consistency_check.py` still
4 PASS / 0 WARN / 0 FAIL, and `paths_from_carcinogen` for `HCB` now
correctly traverses the new `AGONIZES -> INDUCES` chain (e.g.
`HCB -[AGONIZES]-> AHR -[INDUCES]-> CYP1A1 -[PRODUCES]-> BaP_epoxide -> ...`).
