SHELL := /bin/sh
PYTHON ?= python3
RUNNER := $(PYTHON) run_tests.py

.PHONY: test smoke nightly extended list coverage strict-coverage

test: smoke

smoke:
	$(RUNNER) --tier smoke

nightly:
	$(RUNNER) --tier nightly

extended:
	$(RUNNER) --tier extended

list:
	$(RUNNER) --all --list

coverage:
	$(RUNNER) --coverage

# This target is deliberately red until every inventory entry is supported by
# qualified tests.  It is the final migration gate, not the routine nightly.
strict-coverage:
	$(RUNNER) --coverage --require-complete-coverage
