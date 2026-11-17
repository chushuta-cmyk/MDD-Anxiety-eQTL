# Set PYTHON to an interpreter that has pytest, e.g. make check PYTHON=/usr/local/bin/python3
PYTHON  ?= python3
CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra -std=c99

.PHONY: all c-engine check clean

all: c-engine

c-engine: scripts/09a_stream_gse213982_pseudobulk.c scripts/09b_extract_cell_totals.c
	$(CC) $(CFLAGS) scripts/09a_stream_gse213982_pseudobulk.c -o scripts/stream_pseudobulk
	$(CC) $(CFLAGS) scripts/09b_extract_cell_totals.c -o scripts/extract_cell_totals
	@echo "C engines built."

check:
	$(CC) $(CFLAGS) -fsyntax-only scripts/09a_stream_gse213982_pseudobulk.c
	$(CC) $(CFLAGS) -fsyntax-only scripts/09b_extract_cell_totals.c
	$(PYTHON) -m pytest -q tests/

clean:
	rm -f scripts/stream_pseudobulk scripts/extract_cell_totals
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
