# Makefile for HiC-TAD-Library

.PHONY: help install run clean setup

help:
	@echo "Available commands:"
	@echo "  install   Install python dependencies"
	@echo "  run       Run the visualization script"
	@echo "  clean     Remove all generated PNG files"
	@echo "  setup     Install dependencies and run visualizations"

install:
	pip install cooler cooltools bioframe matplotlib pandas numpy

run:
	python visualize_tads.py

clean:
	rm -f media/*.png

setup: install run
