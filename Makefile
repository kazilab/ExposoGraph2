PYTHON ?= python
PYTEST ?= $(PYTHON) -m pytest
SPHINX ?= $(PYTHON) -m sphinx
MYPY ?= $(PYTHON) -m mypy
RUFF ?= ruff

.PHONY: test test-cov lint typecheck docs docs-html ci ci-advisory

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

ci:
	@status=0; \
	$(MAKE) test || status=$$?; \
	$(MAKE) docs || status=$$?; \
	$(MAKE) test-cov || status=$$?; \
	$(MAKE) lint || status=$$?; \
	$(MAKE) typecheck || status=$$?; \
	exit $$status

# Historical alias: every gate is now required, so this runs the same checks
# as ``make ci``. Kept so existing local workflows keep working.
ci-advisory: ci
