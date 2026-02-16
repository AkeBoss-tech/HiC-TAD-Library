# Makefile for HiC-TAD-Library

.PHONY: help install run clean setup test test-unit test-integration test-coverage test-verbose

help:
	@echo "Available commands:"
	@echo "  install          Install python dependencies"
	@echo "  run              Run the visualization script"
	@echo "  clean            Remove all generated PNG files"
	@echo "  setup            Install dependencies and run visualizations"
	@echo "  test             Run all tests"
	@echo "  test-unit        Run only unit tests"
	@echo "  test-integration Run only integration tests"
	@echo "  test-coverage    Run tests with coverage report"
	@echo "  test-verbose     Run tests with verbose output"

install-alphagenome:
	rm -rf ./externals/alphagenome
	pip install -e ./externals/alphagenome
	pip install python-dotenv

install: install-alphagenome
	pip install cooler cooltools bioframe matplotlib pandas numpy

run:
	python visualize_tads.py
	python visualize_boundaries.py
	python visualize_alphagenome.py

clean:
	rm -f media/*.png
	rm -rf htmlcov/
	rm -f .coverage
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true

test:
	pytest

test-unit:
	pytest -m unit

test-integration:
	pytest -m integration

test-coverage:
	pytest --cov=src --cov-report=html --cov-report=term

test-verbose:
	pytest -vv

setup: install run
