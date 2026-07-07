# ExposoGraph 2.0 Manuscript Update Checklist

This checklist is for human manuscript authors. It does not state that the
manuscript has already been updated.

## Methods

- [ ] Distinguish Module 3 as the simple individual-carcinogen workflow and
  Module 5 as the advanced multi-carcinogen mechanism-resolved interaction
  workflow.
- [ ] State that Module 3 remains in v2.0 and is not deprecated.
- [ ] State that Module 5 does not replace Module 3 in v2.0.
- [ ] Describe Module 5 inhibition modes where implemented: competitive, pure
  non-competitive, uncompetitive, and mixed.
- [ ] Describe centralized Ki/Km/IC50 handling, including warning, provenance,
  uncertainty/status, and review-required behavior.
- [ ] Describe GSH activation scaling and fallback behavior.
- [ ] Present the adjusted-risk formula:
  `baseline_relative_risk = BASELINE_RISK_SCORES[carcinogen] * exposure_multiplier`;
  `final_mechanism_multiplier = induction_multiplier * inhibition_burden_multiplier * matrix_gsh_penalty`;
  `adjusted_relative_risk = baseline_relative_risk * final_mechanism_multiplier`.
- [ ] State that diagnostic biological output is not authoritative adjusted
  risk.
- [ ] Describe Shapley/eight-state mechanism attribution across induction,
  competition, and GSH states.

## Results

- [ ] Update Figure 3 text to match the regenerated outputs from WP10.
- [ ] Interpret Figure 3 as a Module 5 figure.
- [ ] Treat the pairwise synergy heatmap as descriptive if retained.
- [ ] Refer to Shapley/eight-state fields as the mechanism-attribution
  interpretation where implemented.

## Limitations

- [ ] State that HIF1A/VEGF gene-level treatment and cobalt-HIF-VEGF pathway
  modeling are deferred.
- [ ] State that expanded edge ontology design is deferred.
- [ ] State the GSH v2 model boundaries, including local bounded behavior,
  fallbacks, warnings, and absence of a full dynamic redox pathway model.
- [ ] State that population modeling, population-genomics expansion, broad
  adapter expansion, full knowledge-graph migration of runtime parameters, and
  new external data-source integrations are deferred to v3 or later unless
  explicitly supported by the accepted v2 package.

## Software Availability And Reproducibility

- [ ] Include the local Figure 3 regeneration command:
  `python tools/generate_figure3.py`.
- [ ] Include the local platform integration check:
  `python tools/local_platform_integration_check.py`.
- [ ] Identify the Module 3 example:
  `python examples/module3_simple_workflow.py`.
- [ ] Identify the Module 5 example:
  `python examples/module5_advanced_interaction_workflow.py`.

## Supplement

- [ ] Document model-card fields for Module 5 transparency.
- [ ] Document mechanism-resolved risk fields, warning fields, provenance
  fields, and review-required flags.
- [ ] Document Figure 3 CSV outputs and the distinction between descriptive
  heatmap values and mechanism-attribution fields.
