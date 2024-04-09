.PHONY: lint
lint:
	ruff --exit-zero check .
	snakefmt --check .

.PHONY: format
format:
	ruff --fix .
	ruff format .
	snakefmt .

.PHONY: pre-commit
pre-commit:
	pre-commit run --all-files

.PHONY: run-demo-pipeline
run-demo-pipeline:
	snakemake --snakefile Snakefile --software-deployment-method conda -n
