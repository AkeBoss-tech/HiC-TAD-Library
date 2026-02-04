# Makefile for HiC-TAD-Library

.PHONY: help install run clean setup

help:
	@echo "Available commands:"
	@echo "  install   Install python dependencies"
	@echo "  run       Run the visualization script"
	@echo "  clean     Remove all generated PNG files"
	@echo "  setup     Install dependencies and run visualizations"

install-alphagenome:
	rm -rf ./externals/alphagenome
	pip install -e ./externals/alphagenome
	pip install python-dotenv

install: install-alphagenome
	pip install cooler cooltools bioframe matplotlib pandas numpy

run:
	python visualize_tads.py
	python visualize_alphagenome.py

clean:
	rm -f media/*.png

setup: install run
