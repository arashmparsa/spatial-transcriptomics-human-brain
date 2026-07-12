# Makefile for the spatialLIBD DLPFC Visium pipeline.
# Steps read/write intermediate AnnData under results/ so they are resumable.

PYTHON ?= python
CONFIG ?= config.yaml

.PHONY: setup download qc normalize cluster svg deconv all clean

setup:
	conda env create -f environment.yml || conda env update -f environment.yml

download:
	$(PYTHON) scripts/download_data.py --config $(CONFIG)

qc:
	$(PYTHON) run_pipeline.py --config $(CONFIG) --steps load qc

normalize:
	$(PYTHON) run_pipeline.py --config $(CONFIG) --steps normalize

cluster:
	$(PYTHON) run_pipeline.py --config $(CONFIG) --steps cluster

svg:
	$(PYTHON) run_pipeline.py --config $(CONFIG) --steps svg

deconv:
	$(PYTHON) run_pipeline.py --config $(CONFIG) --steps deconv

all:
	$(PYTHON) run_pipeline.py --config $(CONFIG) --steps load qc normalize cluster svg

clean:
	rm -rf results/*.h5ad results/figures/* results/*.csv
