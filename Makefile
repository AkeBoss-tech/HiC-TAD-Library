# Makefile for HiC-TAD-Library

.PHONY: help install run clean setup run-haircell

help:
	@echo "Available commands:"
	@echo "  install         Install python dependencies"
	@echo "  run             Run the visualization scripts"
	@echo "  run-haircell    Run the hair cell integration workflow"
	@echo "  clean           Remove all generated PNG files"
	@echo "  setup           Install dependencies and run visualizations"

install-alphagenome:
	rm -rf ./externals/alphagenome
	pip install -e ./externals/alphagenome
	pip install python-dotenv

install: install-alphagenome
	pip install cooler cooltools bioframe matplotlib pandas numpy scipy pyliftover

run:
	python visualize_tads.py
	python visualize_boundaries.py
	python visualize_alphagenome.py

run-haircell:
	python haircell_workflow.py $(GENE)

clean:
	rm -f media/*.png
	rm -f media/haircell/*.png

setup: install run
