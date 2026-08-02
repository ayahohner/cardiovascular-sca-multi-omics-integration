from __future__ import annotations

import argparse
import json
from pathlib import Path

from .pipeline import run_pipeline, verify_release
from .sources import (
    fetch_assets,
    load_analysis_config,
    preflight_assets,
)


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(prog="hspc-fate-mvp")
    command.add_argument("--root", type=Path, default=Path.cwd())
    subcommands = command.add_subparsers(dest="command", required=True)
    subcommands.add_parser("fetch", help="download required open HSPC archives")
    subcommands.add_parser("preflight", help="verify required HSPC archive sizes and checksums")
    run = subcommands.add_parser("run", help="build the donor-level HSPC fate MVP release")
    run.add_argument("--force-preprocess", action="store_true")
    subcommands.add_parser("verify", help="verify HSPC release hashes and scientific controls")
    return command


def main() -> None:
    arguments = parser().parse_args()
    root = arguments.root.resolve()
    config = load_analysis_config(root / "research/hspc-fate/analysis.json")
    if arguments.command == "fetch":
        downloaded = fetch_assets(config, root)
        print(json.dumps({"downloaded": [str(path) for path in downloaded]}, indent=2))
    elif arguments.command == "preflight":
        result = preflight_assets(config, root)
        print(result.to_json(orient="records", indent=2))
        if not result["status"].eq("admitted").all():
            raise SystemExit(1)
    elif arguments.command == "run":
        print(json.dumps(run_pipeline(root, force_preprocess=arguments.force_preprocess), indent=2))
    elif arguments.command == "verify":
        result = verify_release(root)
        print(json.dumps(result, indent=2))
        if result["status"] != "passed":
            raise SystemExit(1)


if __name__ == "__main__":
    main()
