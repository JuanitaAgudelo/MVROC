.PHONY: help install get_gdrive_data get_public_data

# cwd for these targets is pipeline/ (same as the stage drivers).

help:
	@echo "MVROC pipeline"
	@echo ""
	@echo "  make install          - Install Python deps from .bin/requirements.txt"
	@echo "  make get_gdrive_data  - Download pipeline data from Google Drive"
	@echo "  make get_public_data  - Optional: refresh obs + MPCORB from public URLs"
	@echo ""

install:
	python3 -m pip install -r .bin/requirements.txt

get_public_data:
	python3 .bin/get_public_data.py

get_gdrive_data:
	python3 .bin/get_gdrive_data.py
