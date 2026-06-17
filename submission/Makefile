SHELL := /usr/bin/env bash

PAPER = paper
TEX   = pdflatex
BIB   = bibtex
LATEX_FLAGS = -interaction=nonstopmode -halt-on-error

REQUIRED_FILES := $(PAPER).tex bibliography.bib jnl.cls mathphys.bst compile.sh

.DEFAULT_GOAL := all

.PHONY: all help debug rebuild clean

all: $(PAPER).pdf ## Build paper.pdf

help: ## Show available targets
	@printf "Usage: make [target]\\n\\nTargets:\\n"
	@awk 'BEGIN {FS = ":.*## "}; /^[[:alnum:]_.-]+:.*## / {printf "  %-10s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

$(PAPER).pdf: $(REQUIRED_FILES)
	PAPER="$(PAPER)" LATEX="$(TEX)" BIBTEX="$(BIB)" LATEX_FLAGS="$(LATEX_FLAGS)" bash compile.sh

debug: ## Check tools, force a rebuild, and print log diagnostics
	@status=0; \
	echo "=== Build configuration ==="; \
	echo "PAPER=$(PAPER)"; \
	echo "TEX=$(TEX)"; \
	echo "BIB=$(BIB)"; \
	echo "LATEX_FLAGS=$(LATEX_FLAGS)"; \
	echo; \
	echo "=== Toolchain ==="; \
	for tool in "$(TEX)" "$(BIB)"; do \
		if command -v "$$tool" >/dev/null 2>&1; then \
			printf "%s: " "$$tool"; command -v "$$tool"; \
			"$$tool" --version | sed -n '1p'; \
		else \
			echo "MISSING: $$tool"; status=1; \
		fi; \
	done; \
	echo; \
	echo "=== Required files ==="; \
	for file in $(REQUIRED_FILES); do \
		if [ -f "$$file" ]; then \
			echo "ok: $$file"; \
		else \
			echo "missing: $$file"; status=1; \
		fi; \
	done; \
	if [ "$$status" -ne 0 ]; then exit "$$status"; fi; \
	echo; \
	echo "=== Rebuild ==="; \
	$(MAKE) --always-make all; status=$$?; \
	echo; \
	echo "=== LaTeX/BibTeX diagnostics ==="; \
	if [ -f "$(PAPER).log" ]; then \
		grep -nE '(^!|LaTeX Warning|Package .* Warning|Citation .*undefined|Reference .*undefined|There were undefined|Rerun to get)' "$(PAPER).log" || true; \
	else \
		echo "missing: $(PAPER).log"; \
	fi; \
	if [ -f "$(PAPER).blg" ]; then \
		grep -nE '(Warning|Error|---line|I found no)' "$(PAPER).blg" || true; \
	fi; \
	if [ "$$status" -ne 0 ] && [ -f "$(PAPER).log" ]; then \
		echo; \
		echo "=== Last 80 lines of $(PAPER).log ==="; \
		tail -n 80 "$(PAPER).log"; \
	fi; \
	exit "$$status"

rebuild: ## Remove generated files and rebuild
	$(MAKE) clean
	$(MAKE) all

clean: ## Remove generated LaTeX outputs
	rm -f $(PAPER).pdf $(PAPER).aux $(PAPER).bbl $(PAPER).blg \
	      $(PAPER).log $(PAPER).out $(PAPER).toc \
	      Meth.aux Meth.bbl Meth.blg \
	      Supp.aux Supp.bbl Supp.blg
