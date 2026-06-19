"""Main Snakemake workflow for single-cell preprocessing pipeline."""

import sys
from pathlib import Path

# Add scripts and utils directories to path
workflow_dir = Path(workflow.basedir)
sys.path.insert(0, str(workflow_dir / "scripts"))
sys.path.insert(0, str(workflow_dir.parent / "utils"))

# Import helper functions
from parse_config import (
    get_enabled_steps,
    get_step_method,
    get_output_dir,
    get_resources,
    get_hpc_config
)
from build_targets import build_all_targets
from custom_logger import custom_logger

# Parse config to determine enabled steps
ENABLED_STEPS = get_enabled_steps(config)
OUTPUT_DIR = get_output_dir(config)
RESOURCES = get_resources(config)
HPC = get_hpc_config(config)

# Log enabled steps
custom_logger.info("="*60)
custom_logger.info("Single-Cell Preprocessing Pipeline")
custom_logger.info("="*60)
custom_logger.info(f"Project: {config.get('project_name', 'Unnamed')}")
custom_logger.info(f"Output directory: {OUTPUT_DIR}")
custom_logger.info(f"Enabled steps: {', '.join(ENABLED_STEPS) if ENABLED_STEPS else 'None'}")
custom_logger.info("="*60)

# Validate at least one step is enabled
if not ENABLED_STEPS:
    raise ValueError(
        "No pipeline steps are enabled! "
        "Please enable at least one step in your configuration."
    )

# Conditionally include rule modules based on enabled steps
if any(step in ENABLED_STEPS for step in ["cellranger_gex", "cellranger_atac", "cellranger_arc", "cellranger_multi"]):
    include: "rules/cellranger.smk"
    # Phase 1: Include object creation rules after cellranger
    include: "rules/object_creation.smk"
    # Phase 5: Include batch aggregation rules
    include: "rules/batch_aggregation.smk"

if "demultiplexing" in ENABLED_STEPS:
    include: "rules/demultiplexing.smk"

if "doublet_detection" in ENABLED_STEPS:
    include: "rules/doublet_detection.smk"


# Master rule - collects all final outputs
rule all:
    input:
        build_all_targets(config, ENABLED_STEPS)
    default_target: True
    message:
        "Pipeline complete! All enabled steps have finished successfully."


# Utility rule to show pipeline configuration
rule show_config:
    run:
        import yaml
        custom_logger.info("="*60)
        custom_logger.info("Current Pipeline Configuration")
        custom_logger.info("="*60)
        custom_logger.info(yaml.dump(config, default_flow_style=False, sort_keys=False))
        custom_logger.info("="*60)


# Utility rule to clean outputs
rule clean:
    shell:
        """
        echo "Cleaning output directory: {OUTPUT_DIR}"
        rm -rf {OUTPUT_DIR}/*
        echo "Clean complete!"
        """
