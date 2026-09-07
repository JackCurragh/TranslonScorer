.PHONY: help clean clean-build clean-pyc clean-test lint format typecheck gate test test-cov build dist preflight bump-patch bump-minor bump-major

# Every target runs through one interpreter. Override for a venv or a specific
# version:  make test PYTHON=.venv/bin/python
#
# Tools are invoked as `$(PYTHON) -m <tool>`, never as bare executables: a
# `pip install --user` puts them in a scripts dir that is usually not on PATH,
# which is what made `make bump-*` fail with "command not found" (docs/RELEASE.md
# §1). The module form is equivalent and PATH-independent.
PYTHON ?= python3

help:
	@echo "PYTHON=$(PYTHON)"
	@echo ""
	@echo "gate         run the refactor regression gate (GAPDH golden + scalar≡vec)"
	@echo "clean        remove build, test and Python artifacts"
	@echo "lint         run ruff + black --check (the CI lint gate)"
	@echo "format       apply black + ruff --fix"
	@echo "typecheck    run mypy (report-only, as in CI)"
	@echo "test         run the test suite"
	@echo "test-cov     run tests with coverage"
	@echo "build/dist   build sdist + wheel"
	@echo "bump-patch   bump patch version, commit and tag (vX.Y.Z)"
	@echo "bump-minor   bump minor version, commit and tag"
	@echo "bump-major   bump major version, commit and tag"
	@echo "preflight    build + twine check + wheel-name check (what CI asserts)"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +

clean-pyc:
	find . -name '*.pyc' -delete
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .pytest_cache .mypy_cache coverage.xml htmlcov/

gate:
	$(PYTHON) -m pytest tests/test_golden.py tests/test_bam_provider.py -q

lint:
	$(PYTHON) -m ruff check TranslonScorer tests
	$(PYTHON) -m black --check TranslonScorer tests

format:
	$(PYTHON) -m black TranslonScorer tests
	$(PYTHON) -m ruff check --fix TranslonScorer tests

typecheck:
	$(PYTHON) -m mypy --config-file mypy.ini TranslonScorer

test:
	$(PYTHON) -m pytest -q -n auto

test-cov:
	$(PYTHON) -m pytest -n auto --cov=TranslonScorer --cov-report=term-missing --cov-report=xml

build dist: clean
	$(PYTHON) -m build
	ls -l dist

bump-patch:
	$(PYTHON) -m bumpversion patch

bump-minor:
	$(PYTHON) -m bumpversion minor

bump-major:
	$(PYTHON) -m bumpversion major

# Mirrors ci.yml's `build` job exactly. There is deliberately no TestPyPI
# target: it exercises neither OIDC nor the trusted-publisher binding, so it
# cannot rehearse the real path (docs/RELEASE.md §4).
preflight: dist
	$(PYTHON) -m twine check dist/*
	@for f in dist/*.whl dist/*.tar.gz; do \
		case "$$(basename $$f)" in \
			translonscorer-*) ;; \
			*) echo "ERROR: $$(basename $$f) is not normalised; PyPI needs the 'translonscorer-' prefix"; exit 1 ;; \
		esac; \
	done
	@echo "dist filenames normalised:"; ls -1 dist
