#!/usr/bin/env python3
"""
ExposoGraph Population-Scale Batch Runner
============================================
Parallelized execution of ExposoGraph patient_risk_query across
an entire All of Us cohort with checkpointing, progress tracking,
and memory-efficient batch processing.

Designed for:
- All of Us Researcher Workbench (GCP Jupyter, 4-96 CPUs)
- Local testing with synthetic cohorts

Authors: Kenneth J. Pienta (JHU), Julhash U. Kazi (Lund University)
Date: March 2026
"""

import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, TypeAlias, cast

JsonDict: TypeAlias = dict[str, Any]
ConfigDict: TypeAlias = dict[str, Any]
ParticipantRecord: TypeAlias = dict[str, Any]
ParticipantRecords: TypeAlias = list[ParticipantRecord]
SummaryRecord: TypeAlias = dict[str, Any]
ErrorRecord: TypeAlias = dict[str, Any]
WorkerInitArgs: TypeAlias = dict[str, Any]
WorkerArgs: TypeAlias = tuple[ParticipantRecord, str, WorkerInitArgs]
WorkerResult: TypeAlias = dict[str, Any]

# ══════════════════════════════════════════════════════════
#  CONFIGURATION
# ══════════════════════════════════════════════════════════

DEFAULT_CONFIG: ConfigDict = {
    "batch_size": 500,          # Participants per batch
    "max_workers": 4,           # Parallel worker processes
    "checkpoint_interval": 100,  # Save checkpoint every N participants
    "output_dir": "./simulation_output",
    "tissues": ["Liver"],        # Target tissues per participant
    "timeout_per_participant": 30,  # seconds
    "save_individual_results": False,  # Save per-participant JSON
    "save_summary_only": True,   # Only save aggregated results
    "log_level": "INFO",         # INFO, DEBUG, ERROR
}

_HIGH_RISK_LABELS = {"HIGH", "VERY HIGH", "ELEVATED"}


# ══════════════════════════════════════════════════════════
#  SINGLE-PARTICIPANT WORKER
# ══════════════════════════════════════════════════════════

def _run_single_participant(args: WorkerArgs) -> WorkerResult:
    """
    Worker function for a single participant.
    Runs in a separate process via ProcessPoolExecutor.

    Parameters
    ----------
    args : tuple
        (participant_dict, tissue, extensions_init_args)

    Returns
    -------
    dict : {person_id, status, result_summary, error}
    """
    participant, tissue, init_args = args
    del init_args
    pid = participant.get("person_id", "unknown")

    try:
        from ..unified_api import patient_risk_query

        genotypes = cast(dict[str, Any], participant.get("genotypes", {}))
        lifestyle = dict(cast(dict[str, Any], participant.get("lifestyle", {})))
        if "occupational_exposure" not in lifestyle:
            lifestyle["occupational_exposure"] = bool(lifestyle.get("occupational_risk", False))
        if "occupational_risk" not in lifestyle:
            lifestyle["occupational_risk"] = bool(lifestyle.get("occupational_exposure", False))
        scenario = participant.get("exposure_scenario", "general_population")

        # Run the full patient query
        result = patient_risk_query(
            genotypes=genotypes,
            tissue=tissue,
            lifestyle=lifestyle,
            exposure_scenario=scenario,
            chip_status=participant.get("chip_status"),
        )

        # Extract summary metrics (keep output compact for large cohorts)
        summary = _extract_summary_metrics(result, participant)

        return {
            "person_id": pid,
            "status": "success",
            "summary": summary,
            "error": None,
        }

    except Exception as e:
        return {
            "person_id": pid,
            "status": "error",
            "summary": None,
            "error": f"{type(e).__name__}: {str(e)}",
        }


def _get_attr_or_key(value: object, key: str) -> Any | None:
    """Return a field from either an object or mapping."""
    if isinstance(value, dict):
        return value.get(key)
    return getattr(value, key, None)


def _is_high_risk_label(value: object) -> bool:
    """Return whether a serialized risk label denotes high concern."""
    normalized = str(value or "").replace("RiskClassification.", "").strip().upper()
    return normalized in _HIGH_RISK_LABELS


def _extract_flux_summary_from_evidence(evidence: object) -> dict[str, JsonDict]:
    """Build compact flux-class summaries from unified-api evidence entries."""
    if not isinstance(evidence, list):
        return {}

    flux_summary: dict[str, JsonDict] = {}
    for item in evidence:
        class_name = _get_attr_or_key(item, "carcinogen_class")
        if not class_name:
            continue
        flux_summary[class_name] = {
            "net_ratio": _get_attr_or_key(item, "net_ratio"),
            "susceptibility_score_log2": _get_attr_or_key(item, "susceptibility_score_log2"),
            "reactive_intermediate_uM": _get_attr_or_key(item, "reactive_intermediate_uM"),
            "time_to_steady_state_days": _get_attr_or_key(
                item,
                "time_to_steady_state_days",
            ),
            "risk": _get_attr_or_key(item, "risk_classification"),
            "model_kind": _get_attr_or_key(item, "model_kind"),
            "parameter_source": _get_attr_or_key(item, "parameter_source"),
        }
    return flux_summary


def _extract_summary_metrics(result: object, participant: ParticipantRecord) -> SummaryRecord:
    """
    Extract compact summary metrics from a patient_risk_query result.
    Keeps memory manageable for 400K+ participants.
    """
    summary: SummaryRecord = {
        "person_id": participant.get("person_id"),
        "ancestry": participant.get("ancestry"),
        "exposure_scenario": participant.get("exposure_scenario"),
    }

    interactions = getattr(result, "interactions", None)
    if interactions is None and isinstance(result, dict):
        inter = result.get("interactions", {})
        if isinstance(inter, dict) and "error" not in inter:
            summary["interaction_factor"] = inter.get("interaction_factor")
            summary["total_independent_risk"] = inter.get("total_independent_risk")
            summary["total_interaction_risk"] = inter.get("total_interaction_risk")

            gsh = inter.get("gsh_status", {})
            if isinstance(gsh, dict):
                summary["gsh_fraction"] = gsh.get("fraction_normal")
                summary["gsh_tipping_point"] = gsh.get("tipping_point_reached", False)
                summary["gsh_mM"] = gsh.get("steady_state_gsh_mM")
            else:
                summary["gsh_fraction"] = None
                summary["gsh_tipping_point"] = None
                summary["gsh_mM"] = None
        else:
            summary["interaction_factor"] = None
            summary["total_independent_risk"] = None
            summary["total_interaction_risk"] = None
            summary["gsh_fraction"] = None
            summary["gsh_tipping_point"] = None
            summary["gsh_mM"] = None
    else:
        summary["interaction_factor"] = getattr(interactions, "interaction_factor", None)
        summary["total_independent_risk"] = getattr(interactions, "total_independent_risk", None)
        summary["total_interaction_risk"] = getattr(interactions, "total_interaction_risk", None)
        gsh = getattr(interactions, "gsh_status", None)
        summary["gsh_fraction"] = getattr(gsh, "fraction_normal", None)
        summary["gsh_tipping_point"] = getattr(gsh, "tipping_point_reached", False)
        summary["gsh_mM"] = getattr(gsh, "steady_state_gsh_mM", None)

    flux_summary = _extract_flux_summary_from_evidence(
        getattr(result, "flux_class_evidence", None)
    )
    if not flux_summary and isinstance(result, dict):
        flux_summary = _extract_flux_summary_from_evidence(result.get("flux_class_evidence"))

    flux_profile = getattr(result, "flux_profile", None)
    if not flux_summary and flux_profile is None and isinstance(result, dict):
        flux = result.get("flux_profile", {})
        per_class = flux.get("per_class_results", flux) if isinstance(flux, dict) else {}
        if isinstance(per_class, dict):
            for cls_name, cls_result in per_class.items():
                if isinstance(cls_result, dict) and "net_ratio" in cls_result:
                    ss = cls_result.get("steady_state_concentrations_uM", {})
                    ss_model = cls_result.get("steady_state_model", {})
                    flux_summary[cls_name] = {
                        "net_ratio": cls_result.get("net_ratio"),
                        "susceptibility_score_log2": cls_result.get("susceptibility_score_log2"),
                        "reactive_intermediate_uM": (
                            ss.get("reactive_intermediate_uM") if isinstance(ss, dict) else None
                        ),
                        "time_to_steady_state_days": (
                            ss_model.get("time_to_steady_state_days")
                            if isinstance(ss_model, dict)
                            else None
                        ),
                        "risk": cls_result.get("risk_classification"),
                        "model_kind": cls_result.get("model_kind"),
                        "parameter_source": cls_result.get("parameter_source"),
                    }
    elif not flux_summary:
        per_class = (
            getattr(flux_profile, "per_class_results", {}) if flux_profile is not None else {}
        )
        for cls_name, cls_result in per_class.items():
            risk = getattr(cls_result, "risk_classification", None)
            ss = getattr(cls_result, "steady_state_concentrations_uM", {})
            ss_model = getattr(cls_result, "steady_state_model", {})
            flux_summary[cls_name] = {
                "net_ratio": getattr(cls_result, "net_ratio", None),
                "susceptibility_score_log2": getattr(cls_result, "susceptibility_score_log2", None),
                "reactive_intermediate_uM": (
                    ss.get("reactive_intermediate_uM") if isinstance(ss, dict) else None
                ),
                "time_to_steady_state_days": (
                    ss_model.get("time_to_steady_state_days")
                    if isinstance(ss_model, dict)
                    else None
                ),
                "risk": getattr(risk, "value", risk),
                "model_kind": getattr(cls_result, "model_kind", None),
                "parameter_source": getattr(cls_result, "parameter_source", None),
            }
    summary["flux_classes"] = flux_summary

    # Count high-risk pathways
    high_risk = [
        k for k, v in flux_summary.items()
        if isinstance(v, dict) and _is_high_risk_label(v.get("risk"))
    ]
    measured_high_risk = [
        k for k in high_risk
        if flux_summary.get(k, {}).get("model_kind") == "measured_kinetics"
    ]
    proxy_high_risk = [
        k for k in high_risk
        if flux_summary.get(k, {}).get("model_kind")
        and flux_summary.get(k, {}).get("model_kind") != "measured_kinetics"
    ]
    summary["n_high_risk_pathways"] = len(high_risk)
    summary["high_risk_pathways"] = high_risk
    summary["n_measured_high_risk_pathways"] = len(measured_high_risk)
    summary["measured_high_risk_pathways"] = measured_high_risk
    summary["n_proxy_high_risk_pathways"] = len(proxy_high_risk)
    summary["proxy_high_risk_pathways"] = proxy_high_risk

    # Critical warnings
    warnings = getattr(result, "critical_warnings", None)
    if warnings is None and isinstance(result, dict):
        warnings = result.get("critical_warnings", [])
    if isinstance(warnings, list):
        summary["n_critical_warnings"] = len(warnings)
        summary["critical_severities"] = [
            warning.get("severity")
            if isinstance(warning, dict)
            else getattr(warning, "severity", None)
            for warning in warnings[:5]
        ]
    else:
        summary["n_critical_warnings"] = 0

    # Key genotypes for downstream analysis
    geno = cast(dict[str, Any], participant.get("genotypes", {}))
    summary["gstm1"] = geno.get("GSTM1")
    summary["gstt1"] = geno.get("GSTT1")
    summary["aldh2"] = geno.get("ALDH2")
    summary["nat2"] = geno.get("NAT2")
    summary["cyp1a1"] = geno.get("CYP1A1")
    summary["cyp2d6"] = geno.get("CYP2D6")
    summary["cyp2e1"] = geno.get("CYP2E1")

    # CHIP somatic-modifier summary
    chip_status = getattr(result, "chip_status", None)
    chip_effect = getattr(result, "chip_effect", None)
    if chip_status is not None:
        summary["chip_positive"] = bool(getattr(chip_status, "present", False))
        summary["chip_max_vaf"] = getattr(chip_status, "max_vaf", 0.0)
    else:
        summary["chip_positive"] = False
        summary["chip_max_vaf"] = 0.0
    if chip_effect is not None:
        summary["chip_aggregate_modifier"] = getattr(chip_effect, "aggregate_modifier", 1.0)
        summary["chip_driver_genes"] = list(getattr(chip_effect, "driver_genes", []))
    else:
        summary["chip_aggregate_modifier"] = 1.0
        summary["chip_driver_genes"] = []
    summary["chip_adjusted_risks"] = dict(getattr(result, "chip_adjusted_risks", {}) or {})

    return summary


# ══════════════════════════════════════════════════════════
#  BATCH RUNNER
# ══════════════════════════════════════════════════════════

class PopulationSimulation:
    """
    Population-scale ExposoGraph simulation engine.

    Processes participants in batches with:
    - Multiprocessing parallelism
    - Checkpoint/resume capability
    - Progress tracking
    - Memory-efficient incremental output

    Usage:
        sim = PopulationSimulation(participants, config)
        sim.run()
        results = sim.get_results()
    """

    def __init__(self, participants: ParticipantRecords, config: ConfigDict | None = None):
        """
        Parameters
        ----------
        participants : list
            List of participant dicts from allofus_adapter
        config : dict, optional
            Override DEFAULT_CONFIG settings
        """
        self.participants = participants
        self.config = {**DEFAULT_CONFIG, **(config or {})}
        self.results: list[SummaryRecord] = []
        self.errors: list[ErrorRecord] = []
        self.start_time: float | None = None
        self.end_time: float | None = None

        # Checkpoint state
        self.output_dir = Path(self.config["output_dir"])
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.results_file = self.output_dir / "results.jsonl"
        self.completed_ids: set[str] = set()

        # Resume from checkpoint if available
        self._load_checkpoint()

    def _load_checkpoint(self) -> None:
        """Load previous run state if checkpoint exists."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file) as f:
                loaded = json.load(f)
            state = cast(JsonDict, loaded) if isinstance(loaded, dict) else {}
            self.completed_ids = {str(pid) for pid in state.get("completed_ids", [])}
            print(f"  Resuming from checkpoint: {len(self.completed_ids)} "
                  f"already completed")

    def _save_checkpoint(self) -> None:
        """Save current progress to checkpoint file."""
        state: JsonDict = {
            "completed_ids": list(self.completed_ids),
            "n_completed": len(self.completed_ids),
            "n_total": len(self.participants),
            "n_errors": len(self.errors),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        with open(self.checkpoint_file, "w") as f:
            json.dump(state, f, indent=2)

    def _append_result(self, result: SummaryRecord) -> None:
        """Append a single result to the JSONL output file."""
        with open(self.results_file, "a") as f:
            f.write(json.dumps(result) + "\n")

    def run(self, tissue: str = "Liver") -> None:
        """
        Execute the population simulation.

        Parameters
        ----------
        tissue : str
            Target tissue for risk assessment
        """
        start_time = time.time()
        self.start_time = start_time

        # Filter out already-completed participants
        remaining = [p for p in self.participants
                     if p["person_id"] not in self.completed_ids]

        n_total = len(self.participants)
        n_remaining = len(remaining)
        n_workers = int(self.config["max_workers"])
        batch_size = int(self.config["batch_size"])
        checkpoint_interval = int(self.config["checkpoint_interval"])

        print("\n" + "=" * 70)
        print("  ExposoGraph Population Simulation")
        print("=" * 70)
        print(f"  Total participants: {n_total}")
        print(f"  Already completed:  {len(self.completed_ids)}")
        print(f"  Remaining:          {n_remaining}")
        print(f"  Workers:            {n_workers}")
        print(f"  Batch size:         {batch_size}")
        print(f"  Target tissue:      {tissue}")
        print(f"  Output:             {self.output_dir}")
        print("=" * 70)

        if n_remaining == 0:
            print("  All participants already processed.")
            self.end_time = time.time()
            return

        # Process in batches
        completed_in_run = 0
        errors_in_run = 0

        for batch_start in range(0, n_remaining, batch_size):
            batch = remaining[batch_start:batch_start + batch_size]
            batch_num = batch_start // batch_size + 1
            total_batches = (n_remaining + batch_size - 1) // batch_size

            print(f"\n  Batch {batch_num}/{total_batches} "
                  f"({len(batch)} participants)")

            # Prepare worker arguments
            worker_args: list[WorkerArgs] = [(p, tissue, {}) for p in batch]

            # Run with ProcessPoolExecutor
            if n_workers > 1:
                batch_results = self._run_parallel(worker_args, n_workers)
            else:
                batch_results = self._run_sequential(worker_args)

            # Process results
            for result in batch_results:
                pid = result["person_id"]
                self.completed_ids.add(pid)
                completed_in_run += 1

                if result["status"] == "success" and result["summary"]:
                    self.results.append(result["summary"])
                    self._append_result(result["summary"])
                else:
                    errors_in_run += 1
                    self.errors.append({
                        "person_id": pid,
                        "error": result.get("error", "Unknown"),
                    })

                # Checkpoint
                if completed_in_run % checkpoint_interval == 0:
                    self._save_checkpoint()
                    elapsed = time.time() - start_time
                    rate = completed_in_run / elapsed if elapsed > 0 else 0
                    eta = (n_remaining - completed_in_run) / rate if rate > 0 else 0
                    print(f"    Progress: {completed_in_run}/{n_remaining} "
                          f"({rate:.1f}/sec, ETA {eta/60:.1f} min)")

        # Final checkpoint
        self._save_checkpoint()
        end_time = time.time()
        self.end_time = end_time

        elapsed = end_time - start_time
        print(f"\n  Completed: {completed_in_run} participants in "
              f"{elapsed:.1f}s ({completed_in_run/elapsed:.1f}/sec)")
        print(f"  Errors: {errors_in_run}")
        print(f"  Results saved to: {self.results_file}")

    def _run_parallel(self, worker_args: list[WorkerArgs], n_workers: int) -> list[WorkerResult]:
        """Run participants in parallel using ProcessPoolExecutor."""
        results: list[WorkerResult] = []
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_run_single_participant, args): args[0]["person_id"]
                for args in worker_args
            }
            for future in as_completed(futures):
                pid = futures[future]
                try:
                    result = future.result(
                        timeout=float(self.config["timeout_per_participant"])
                    )
                    results.append(result)
                except Exception as e:
                    results.append({
                        "person_id": pid,
                        "status": "error",
                        "summary": None,
                        "error": f"Worker error: {str(e)}",
                    })
        return results

    def _run_sequential(self, worker_args: list[WorkerArgs]) -> list[WorkerResult]:
        """Run participants sequentially (for debugging or single-CPU)."""
        results: list[WorkerResult] = []
        for args in worker_args:
            result = _run_single_participant(args)
            results.append(result)
        return results

    def get_results(self) -> list[SummaryRecord]:
        """Return all collected results."""
        return self.results

    def get_errors(self) -> list[ErrorRecord]:
        """Return all error records."""
        return self.errors

    def load_results_from_disk(self) -> list[SummaryRecord]:
        """Load results from the JSONL output file."""
        results: list[SummaryRecord] = []
        if self.results_file.exists():
            with open(self.results_file) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        payload = json.loads(line)
                        if isinstance(payload, dict):
                            results.append(cast(SummaryRecord, payload))
        return results

    def get_run_stats(self) -> JsonDict:
        """Get summary statistics for the run."""
        elapsed = (
            (self.end_time - self.start_time)
            if self.end_time is not None and self.start_time is not None
            else 0.0
        )
        return {
            "n_total": len(self.participants),
            "n_completed": len(self.completed_ids),
            "n_errors": len(self.errors),
            "elapsed_seconds": elapsed,
            "rate_per_second": len(self.completed_ids) / elapsed if elapsed > 0 else 0,
            "output_dir": str(self.output_dir),
        }

    def save_run_report(self) -> None:
        """Save a JSON summary of the entire run."""
        report: JsonDict = {
            "stats": self.get_run_stats(),
            "config": self.config,
            "errors": self.errors[:100],  # cap at 100
        }
        report_file = self.output_dir / "run_report.json"
        with open(report_file, "w") as f:
            json.dump(report, f, indent=2)
        print(f"  Run report saved to: {report_file}")


# ══════════════════════════════════════════════════════════
#  CONVENIENCE FUNCTIONS
# ══════════════════════════════════════════════════════════

def run_synthetic_simulation(
    n_participants: int = 1000,
    n_workers: int = 1,
    tissue: str = "Liver",
    output_dir: str = "./simulation_output",
    seed: int = 42,
) -> PopulationSimulation:
    """
    Run a complete simulation on a synthetic cohort.
    Convenience function for testing the full pipeline.

    Parameters
    ----------
    n_participants : int
    n_workers : int
    tissue : str
    output_dir : str
    seed : int

    Returns
    -------
    PopulationSimulation : the completed simulation object
    """
    from .allofus_adapter import generate_synthetic_cohort

    print("=" * 70)
    print("  ExposoGraph Population Simulation — Synthetic Cohort Test")
    print("=" * 70)

    # Generate synthetic cohort
    participants = generate_synthetic_cohort(n=n_participants, seed=seed)

    # Configure and run
    config = {
        **DEFAULT_CONFIG,
        "max_workers": n_workers,
        "output_dir": output_dir,
        "checkpoint_interval": max(10, n_participants // 20),
    }

    sim = PopulationSimulation(participants, config)
    sim.run(tissue=tissue)
    sim.save_run_report()

    return sim


# ══════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="ExposoGraph Population-Scale Batch Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test with 100 synthetic participants (sequential)
  python -m ExposoGraph.population_simulation.batch_runner --synthetic 100 --workers 1

  # Full synthetic simulation with parallelism
  python -m ExposoGraph.population_simulation.batch_runner --synthetic 1000 --workers 4 --tissue Lung

  # Run on pre-built cohort JSON
  python -m ExposoGraph.population_simulation.batch_runner --input cohort.json --workers 8 --output ./results/

  # Resume a previous run
  python -m ExposoGraph.population_simulation.batch_runner --input cohort.json --output ./results/ --resume
        """
    )
    parser.add_argument("--synthetic", type=int, metavar="N",
                        help="Run on N synthetic participants")
    parser.add_argument("--input", type=str,
                        help="Path to participant cohort JSON file")
    parser.add_argument("--output", type=str, default="./simulation_output",
                        help="Output directory (default: ./simulation_output)")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of parallel workers (default: 1)")
    parser.add_argument("--tissue", type=str, default="Liver",
                        help="Target tissue (default: Liver)")
    parser.add_argument("--batch-size", type=int, default=500,
                        help="Batch size (default: 500)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for synthetic data")
    args = parser.parse_args()

    if args.synthetic:
        sim = run_synthetic_simulation(
            n_participants=args.synthetic,
            n_workers=args.workers,
            tissue=args.tissue,
            output_dir=args.output,
            seed=args.seed,
        )
        return

    if args.input:
        with open(args.input) as f:
            loaded = json.load(f)
        if not isinstance(loaded, list):
            raise ValueError("Participant cohort JSON must be a list of participant records")
        participants = cast(ParticipantRecords, loaded)
        print(f"  Loaded {len(participants)} participants from {args.input}")

        config = {
            **DEFAULT_CONFIG,
            "max_workers": args.workers,
            "output_dir": args.output,
            "batch_size": args.batch_size,
        }
        sim = PopulationSimulation(participants, config)
        sim.run(tissue=args.tissue)
        sim.save_run_report()
        return

    parser.print_help()


if __name__ == "__main__":
    main()
