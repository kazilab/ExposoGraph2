ExposoGraph
===========

**ExposoGraph** is a Python toolkit for building, curating, and exporting
carcinogen metabolism knowledge graphs using LLM-powered extraction and
manual entry.

Part of the **CarcinoGenomic Platform** — a 5-layer computational pipeline
for individualized carcinogen metabolism risk assessment from germline DNA.

Version: **0.0.5**
Developed by: **Data analysis team @ KaziLab**
Contact: **exposograph@kazilab.se**
Copyright: **KaziLab**

The current release separates two control layers:

- **Graph mode** for ingestion: ``exploratory`` or ``strict``
- **Graph visibility** for viewing and export: ``all``, ``validated_only``, or ``exploratory_only``

.. note::

   This project is under active development.

.. note::

   Parts of this documentation and code were created with assistance from
   ChatGPT Codex and Claude Code.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   quickstart
   deployment

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   schema
   gene-panels

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/config
   api/models
   api/engine
   api/graph_analysis
   api/exporter
   api/cytoscape_adapter
   api/viewer_dash
   api/graph_filters
   api/grounding
   api/llm_extractor
   api/db_clients
   api/seeder
   api/interaction_engine
   api/example_graphs
   api/reference_data
   api/storage

.. toctree::
   :maxdepth: 1
   :caption: Project

   production-readiness
   changelog
