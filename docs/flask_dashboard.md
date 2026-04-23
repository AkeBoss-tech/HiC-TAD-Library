# Flask Dashboard Architecture

This repository now includes a first-pass Flask dashboard scaffold for a custom
HTML control room around the regulatory activity pipeline.

## Page Map

- `/`
  Main analysis page. Includes:
  - readiness checks for keys and reference assets
  - preset loci for mouse and human runs
  - reference and perturbation form inputs
  - quick links to recent saved runs

- `/result/<slug>`
  Unified result page for one saved regulatory activity manifest. Shows:
  - generated overview figures
  - top candidate elements
  - predicted affected genes
  - model and data provenance

- `/saved`
  Browser for existing `data/processed/regulatory_activity_*.json` files.

- `/api/result/<slug>`
  Raw JSON response for one saved run.

- `/download/<slug>`
  Download a manifest as an attachment.

## Application Structure

- `flask_dashboard/__init__.py`
  Flask app factory and registration.

- `flask_dashboard/routes.py`
  Route handlers and page orchestration.

- `flask_dashboard/services.py`
  Thin service layer that:
  - builds `LocusRequest`
  - runs the regulatory activity pipeline
  - loads saved manifests
  - exposes readiness metadata and presets

- `flask_dashboard/templates/`
  Jinja templates for the shell and pages.

- `flask_dashboard/static/dashboard.css`
  Shared visual design system for the Flask UI.

- `run_flask_dashboard.py`
  Local development entrypoint.

## Design Direction

The Flask dashboard is meant to replace the fragmented app-by-app Streamlit
experience with one task-oriented product flow:

1. select locus and mode
2. validate environment readiness
3. run reference or perturbation analysis
4. inspect structure, activity, and gene effects together
5. export or revisit saved manifests

## Next Implementation Steps

- add background job execution and polling for long AlphaGenome runs
- replace static figure-only rendering with interactive Plotly panels
- add linked selection between candidates, genes, and genomic plots
- add richer report export and provenance audit views
