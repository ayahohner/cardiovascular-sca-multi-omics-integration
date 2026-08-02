.PHONY: setup fetch preflight mvp verify hspc-fetch hspc-preflight hspc-mvp hspc-verify test test-fast clean-cache

setup:
	python3 -m venv .venv
	.venv/bin/python -m pip install --upgrade pip
	.venv/bin/python -m pip install -r requirements.lock
	.venv/bin/python -m pip install -e . --no-deps

fetch:
	.venv/bin/cardio-sca-mvp fetch

preflight:
	.venv/bin/cardio-sca-mvp preflight

mvp:
	.venv/bin/cardio-sca-mvp run

verify:
	.venv/bin/cardio-sca-mvp verify

hspc-fetch:
	.venv/bin/hspc-fate-mvp fetch

hspc-preflight:
	.venv/bin/hspc-fate-mvp preflight

hspc-mvp:
	.venv/bin/hspc-fate-mvp run

hspc-verify:
	.venv/bin/hspc-fate-mvp verify

test:
	.venv/bin/pytest --cov=cardio_sca_mvp --cov=hspc_fate_mvp --cov-report=term-missing

test-fast:
	.venv/bin/pytest -m "not data"

clean-cache:
	@echo "Remove data/derived manually only when a full reference rebuild is intended."
