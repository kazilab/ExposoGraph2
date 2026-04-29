"""Check and optionally normalize biomarker_mapping.json documents."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .loader import load_json_mapping
from .mapping_document import normalize_mapping_document, validate_mapping_document

DEFAULT_MAPPING_PATH = Path(__file__).resolve().parents[3] / "data" / "biomarker_mapping.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate biomarker mapping structure and enforce trace/update compatibility."
    )
    parser.add_argument(
        "--mapping",
        default=str(DEFAULT_MAPPING_PATH),
        help=f"Path to biomarker_mapping.json (default: {DEFAULT_MAPPING_PATH})",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Normalize mapping in-memory before validating (adds entry_id/trace/update metadata).",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write normalized mapping back to --mapping (requires --fix).",
    )
    args = parser.parse_args(argv)

    if args.write and not args.fix:
        parser.error("--write requires --fix")

    mapping_path = Path(args.mapping)
    document = load_json_mapping(mapping_path)
    if args.fix:
        document = normalize_mapping_document(document)

    issues = validate_mapping_document(document)
    if args.fix and args.write:
        mapping_path.write_text(json.dumps(document, indent=2, sort_keys=False) + "\n", encoding="utf-8")

    if issues:
        print(f"FAILED: {len(issues)} issue(s) found")
        for issue in issues:
            print(f"- {issue}")
        return 1

    print("OK: mapping document is valid, traceable, and update-compatible")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
