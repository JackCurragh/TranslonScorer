.PHONY: help clean clean-build clean-pyc clean-test lint gate test test-cov build dist release-test bump-patch bump-minor bump-major

help:
	@echo "gate         run the refactor regression gate (GAPDH golden + scalar≡vec)"
	@echo "clean        remove build, test and Python artifacts"
	@echo "lint         run mypy type checks"
	@echo "test         run the test suite"
	@echo "test-cov     run tests with coverage"
	@echo "build/dist   build sdist + wheel"
	@echo "bump-patch   bump patch version, commit and tag (vX.Y.Z)"
	@echo "bump-minor   bump minor version, commit and tag"
	@echo "bump-major   bump major version, commit and tag"
	@echo "release-test build then upload to TestPyPI"

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
	python3 -m pytest tests/test_golden.py tests/test_bam_provider.py -q

lint:
	ruff check TranslonScorer tests
	black --check TranslonScorer tests

format:
	black TranslonScorer tests
	ruff check --fix TranslonScorer tests

typecheck:
	mypy --config-file mypy.ini TranslonScorer

test:
	pytest -q -n auto

test-cov:
	pytest -n auto --cov=TranslonScorer --cov-report=term-missing --cov-report=xml

build dist: clean
	python -m build
	ls -l dist

bump-patch:
	bump2version patch

bump-minor:
	bump2version minor

bump-major:
	bump2version major

release-test: dist
	twine check dist/*
	twine upload --repository testpypi dist/*
