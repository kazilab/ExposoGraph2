PYTHON ?= python
PYTEST ?= $(PYTHON) -m pytest
SPHINX ?= $(PYTHON) -m sphinx
MYPY ?= $(PYTHON) -m mypy
RUFF ?= ruff

.PHONY: test test-cov lint typecheck docs docs-html check-biomarker-mapping ci ci-advisory
.PHONY: build-biomarker-mapping compare-biomarker-mapping

test:
	$(PYTEST) --no-cov

test-cov:
	$(PYTEST)

lint:
	$(RUFF) check .

typecheck:
	$(MYPY) ExposoGraph

docs:
	$(SPHINX) -b dummy docs docs/_build/dummy

docs-html:
	$(SPHINX) -b html docs docs/_build/html

check-biomarker-mapping:
	$(PYTHON) -m ExposoGraph._biomarker_scaffold.scripts.registries.check_mapping --mapping ExposoGraph/data/biomarker_mapping.json

build-biomarker-mapping:
	$(PYTHON) -m ExposoGraph._biomarker_scaffold.scripts.registries.build_mapping --source ExposoGraph/_biomarker_scaffold/data/registries/biomarkers_master.yaml --out ExposoGraph/data/biomarker_mapping.json --old ExposoGraph/data/biomarker_mapping_old.json

compare-biomarker-mapping:
	$(PYTHON) -m ExposoGraph._biomarker_scaffold.scripts.registries.build_mapping --compare-only --old ExposoGraph/data/biomarker_mapping_old.json --out ExposoGraph/data/biomarker_mapping.json

ci:
	@status=0; \
	$(MAKE) test || status=$$?; \
	$(MAKE) docs || status=$$?; \
	$(MAKE) check-biomarker-mapping || status=$$?; \
	$(MAKE) test-cov || status=$$?; \
	$(MAKE) lint || status=$$?; \
	$(MAKE) typecheck || status=$$?; \
	exit $$status

# Historical alias: every gate is now required, so this runs the same checks
# as ``make ci``. Kept so existing local workflows keep working.
ci-advisory: ci
