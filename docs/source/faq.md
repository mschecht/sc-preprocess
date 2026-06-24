# FAQ

## My workflow is locked, what do I do?

`Snakemake` will lock a directory for multiple reasons to ensure data does not get overwritten. You can read about it [here](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-lock-the-working-directory).

To unlock the workflow directory pass `--unlock` directly to `Snakemake` through `--snakemake-args` like this: 

```
sc-preprocess run --config-file <your_config.yaml> --cores 1 --snakemake-args --unlock
````

## How do `.done` files track finished steps?

Each pipeline step writes an empty flag file to `00_LOGS/` when it completes successfully — for example, `1_L001_cellranger_gex.done`. Snakemake uses these files as targets: if the `.done` file already exists, the step is skipped on re-run.

To force a step to re-run, delete its `.done` file, or use `--forcerun` (see below).

## How do I re-run the workflow from a specific step?

Pass `--forcerun <rule_name>` to Snakemake via `--snakemake-args`. For example, to restart from `create_gex_anndata`:

```bash
sc-preprocess run --config-file pipeline_config.yaml \
                             --cores 1 \
                             --snakemake-args --forcerun create_gex_anndata
```

To find available rule names, visualize the pipeline DAG:

```bash
sc-preprocess run --config-file pipeline_config.yaml --cores 1 --dag | dot -Tpng > dag.png
```

## `render-docs` fails with "address already in use"

This happens when a previous `render-docs` process was not cleanly stopped and left a zombie sphinx-autobuild process holding the port.

Check what is already running:

```bash
lsof -i :8000
```

Kill all lingering sphinx-autobuild processes:

```bash
pkill -f sphinx-autobuild
```

Then re-run `sc-preprocess render-docs`.

## What happens when a Cell Ranger cluster-mode job is killed?

When `cluster-mode` is enabled with `jobmode: slurm` (or another scheduler), Cell Ranger's Martian runtime saves checkpoints after each completed stage. If the job is killed mid-run, Cell Ranger leaves a `_lock` file in the output directory. You can read more about Cell Ranger cluster mode [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode).

The pipeline automatically removes this `_lock` file when resubmitted, so the job resumes from the last completed stage rather than starting over.
