# Algebraic Genomics -- convenience targets (equivalent python commands work on any OS)
PY ?= python

.PHONY: quick full figures papers inventory clean

quick:        ## papers + fast checks + figures (about 30-40 min)
	$(PY) -u build_all.py --quick

full:         ## papers + every computation (about 4 h)
	$(PY) -u build_all.py --full

figures:      ## the five overview figures only
	$(PY) code/generate_figures.py

papers:       ## compile all papers, run no scripts
	$(PY) -u build_all.py --no-scripts

inventory:    ## list what would be run
	$(PY) build_all.py --inventory

clean:
	$(PY) -c "import glob,os;[os.remove(f) for e in ('aux','log','out') for f in glob.glob('**/*.'+e,recursive=True)]"
