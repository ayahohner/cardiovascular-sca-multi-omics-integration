.PHONY: setup fetch preflight mvp verify test test-fast clean-cache

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

test:
	.venv/bin/pytest --cov=cardio_sca_mvp --cov-report=term-missing

test-fast:
	.venv/bin/pytest -m "not data"

clean-cache:
	@echo "Remove data/derived manually only when a full reference rebuild is intended."
