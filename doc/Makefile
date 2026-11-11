# Documentation-only build; independent of AMPS and SWMF build rules.
PYTHON ?= python3
LATEX_ENGINE ?= auto
.PHONY: all MANUAL manual PDF pdf check clean distclean help
all: PDF
MANUAL manual PDF pdf:
	$(PYTHON) -B build_manual.py --engine "$(LATEX_ENGINE)"
check:
	$(PYTHON) -B check_manual.py
clean:
	$(PYTHON) -B build_manual.py --clean
distclean:
	$(PYTHON) -B build_manual.py --clean --remove-pdf
help:
	@echo "make PDF                   Build doc/AMPS.pdf"
	@echo "make MANUAL                Synonym for make PDF"
	@echo "make PDF LATEX_ENGINE=...  Choose pdflatex/tectonic or an executable path"
	@echo "make check                 Check local inputs, assets and references"
	@echo "make clean                 Remove generated manual build files only"
	@echo "make distclean             Also remove generated AMPS.pdf"
