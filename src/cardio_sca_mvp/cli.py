from __future__ import annotations

import argparse
import json
from pathlib import Path

from .pipeline import run_pipeline, verify_release
from .sources import fetch_missing_sources, load_registry, preflight_sources


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(prog="cardio-sca-mvp")
    command.add_argument("--root", type=Path, default=Path.cwd())
    subcommands = command.add_subparsers(dest="command", required=True)
    subcommands.add_parser("fetch", help="download and verify missing frozen source assets")
    subcommands.add_parser("preflight", help="verify source sizes and SHA-256 checksums")
    run = subcommands.add_parser("run", help="build the MVP release")
    run.add_argument("--force-reference", action="store_true", help="rebuild pseudobulk reference caches")
    subcommands.add_parser("verify", help="verify release artifact hashes and controls")
    return command


def main() -> None:
    arguments = parser().parse_args()
    root = arguments.root.resolve()
    if arguments.command == "fetch":
        registry = load_registry(root / "config/sources.json")
        downloaded = fetch_missing_sources(registry, root)
        print(json.dumps({"downloaded": [str(path) for path in downloaded]}, indent=2))
    elif arguments.command == "preflight":
        registry = load_registry(root / "config/sources.json")
        result = preflight_sources(registry, root)
        print(result.to_json(orient="records", indent=2))
        if not result["status"].eq("admitted").all():
            raise SystemExit(1)
    elif arguments.command == "run":
        release = run_pipeline(root, force_reference=arguments.force_reference)
        print(json.dumps(release, indent=2))
    elif arguments.command == "verify":
        result = verify_release(root)
        print(json.dumps(result, indent=2))
        if result["status"] != "passed":
            raise SystemExit(1)


if __name__ == "__main__":
    main()
