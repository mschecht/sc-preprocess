"""Batch aggregation rules - merge per-capture objects into batch-level objects."""

import os
import sys
from pathlib import Path
from tempfile import gettempdir
import pandas as pd

# Import utilities
sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from sc_preprocess.config_validator import parse_output_directories

# Get output directories
OUTPUT_DIRS = parse_output_directories(config)
ANNDATA_DIR = OUTPUT_DIRS["anndata_dir"]
BATCH_DIR = os.path.join(config.get("output_dir", "output"), "04_BATCH_OBJECTS")
LOGS_DIR = OUTPUT_DIRS["logs_dir"]

# Helper function to get all captures for a batch
def get_captures_for_batch(batch_id, modality):
    """Get all captures for a given batch."""
    if modality == "gex" and config.get("cellranger_gex"):
        libraries_file = config["cellranger_gex"]["libraries"]
    elif modality == "atac" and config.get("cellranger_atac"):
        libraries_file = config["cellranger_atac"]["libraries"]
    elif modality == "arc" and config.get("cellranger_arc"):
        libraries_file = config["cellranger_arc"]["libraries"]
    elif modality == "multi" and config.get("cellranger_multi"):
        libraries_file = config["cellranger_multi"]["libraries"]
    else:
        return []

    df = pd.read_csv(libraries_file, sep="\t")
    captures = df[df['batch'] == int(batch_id)]['capture'].tolist()
    return captures


# ============================================================================
# GEX Batch Aggregation
# ============================================================================

if config.get("cellranger_gex"):

    # Get unique batches
    libraries_file = config["cellranger_gex"]["libraries"]
    gex_df = pd.read_csv(libraries_file, sep="\t")
    GEX_BATCHES = gex_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(GEX_BATCHES)} GEX batch(es)")

    def get_gex_percapture_objects(wildcards):
        """Get all per-capture .h5ad files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "gex")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5ad") for cap in captures]

    def get_gex_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "gex")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_gex_anndata.done") for cap in captures]

    rule aggregate_gex_batch:
        """Aggregate per-capture GEX objects into batch-level object."""
        input:
            h5ads = get_gex_percapture_objects,
            done_flags = get_gex_percapture_done
        output:
            h5ad = os.path.join(BATCH_DIR, "{batch}_gex.h5ad"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_gex_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "gex"
        threads: config["cellranger_gex"].get("batch_aggregation", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_gex"].get("batch_aggregation", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_gex"].get("batch_aggregation", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_gex_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"


# ============================================================================
# ATAC Batch Aggregation
# ============================================================================

if config.get("cellranger_atac"):

    # Get unique batches
    libraries_file = config["cellranger_atac"]["libraries"]
    atac_df = pd.read_csv(libraries_file, sep="\t")
    ATAC_BATCHES = atac_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(ATAC_BATCHES)} ATAC batch(es)")

    def get_atac_percapture_objects(wildcards):
        """Get all per-capture .h5ad files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "atac")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5ad") for cap in captures]

    def get_atac_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "atac")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_atac_anndata.done") for cap in captures]

    rule aggregate_atac_batch:
        """Aggregate per-capture ATAC objects into batch-level object."""
        input:
            h5ads = get_atac_percapture_objects,
            done_flags = get_atac_percapture_done
        output:
            h5ad = os.path.join(BATCH_DIR, "{batch}_atac.h5ad"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_atac_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "atac"
        threads: config["cellranger_atac"].get("batch_aggregation", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_atac"].get("batch_aggregation", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_atac"].get("batch_aggregation", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_atac_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"


# ============================================================================
# ARC Batch Aggregation
# ============================================================================

if config.get("cellranger_arc"):

    # Get unique batches
    libraries_file = config["cellranger_arc"]["libraries"]
    arc_df = pd.read_csv(libraries_file, sep="\t")
    ARC_BATCHES = arc_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(ARC_BATCHES)} ARC batch(es)")

    def get_arc_percapture_objects(wildcards):
        """Get all per-capture .h5mu files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "arc")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5mu") for cap in captures]

    def get_arc_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "arc")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_arc_mudata.done") for cap in captures]

    rule aggregate_arc_batch:
        """Aggregate per-capture ARC objects into batch-level object."""
        input:
            h5mus = get_arc_percapture_objects,
            done_flags = get_arc_percapture_done
        output:
            h5mu = os.path.join(BATCH_DIR, "{batch}_arc.h5mu"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_arc_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "arc"
        threads: config["cellranger_arc"].get("batch_aggregation", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_arc"].get("batch_aggregation", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_arc"].get("batch_aggregation", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_arc_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"


# ============================================================================
# Multi Batch Aggregation
# ============================================================================

if config.get("cellranger_multi"):

    libraries_file = config["cellranger_multi"]["libraries"]
    multi_df = pd.read_csv(libraries_file, sep="\t")
    MULTI_BATCHES = multi_df['batch'].unique().tolist()

    custom_logger.info(f"Batch aggregation: Found {len(MULTI_BATCHES)} multi batch(es)")

    def get_multi_percapture_objects(wildcards):
        """Get all per-capture .h5mu files for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "multi")
        return [os.path.join(ANNDATA_DIR, f"{wildcards.batch}_{cap}.h5mu") for cap in captures]

    def get_multi_percapture_done(wildcards):
        """Get all per-capture done flags for a batch."""
        captures = get_captures_for_batch(wildcards.batch, "multi")
        return [os.path.join(LOGS_DIR, f"{wildcards.batch}_{cap}_multi_mudata.done") for cap in captures]

    rule aggregate_multi_batch:
        """Aggregate per-capture multi objects into batch-level MuData object."""
        input:
            h5mus = get_multi_percapture_objects,
            done_flags = get_multi_percapture_done
        output:
            h5mu = os.path.join(BATCH_DIR, "{batch}_multi.h5mu"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_multi_batch_aggregation.done"))
        params:
            batch = "{batch}",
            modality = "multi"
        threads: config["cellranger_multi"].get("batch_aggregation", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_multi"].get("batch_aggregation", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_multi"].get("batch_aggregation", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_multi_batch_aggregation.log")
        script:
            "../scripts/aggregate_batch.py"


# ============================================================================
# Metadata Enrichment (Final Objects)
# ============================================================================

# Final output directory
FINAL_DIR = OUTPUT_DIRS.get("final_dir")

# Directories for metadata sources
DEMUX_DIR = OUTPUT_DIRS.get("demultiplexing_dir")
DOUBLET_DIR = OUTPUT_DIRS.get("doublet_detection_dir")


if config.get("cellranger_gex"):

    rule enrich_gex_metadata:
        """Enrich GEX batch object with analysis metadata (demux, doublet, annotation)."""
        input:
            batch_object = os.path.join(BATCH_DIR, "{batch}_gex.h5ad"),
            done_flag = os.path.join(LOGS_DIR, "{batch}_gex_batch_aggregation.done")
        output:
            enriched_object = os.path.join(FINAL_DIR, "{batch}_gex.h5ad"),
            obs_table = os.path.join(FINAL_DIR, "{batch}_gex_obs.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_gex_enrichment.done"))
        params:
            batch = "{batch}",
            modality = "gex",
            demux_dir = DEMUX_DIR,
            doublet_dir = DOUBLET_DIR
        threads: config["cellranger_gex"].get("enrichment", {}).get("threads", 4)
        resources:
            mem_mb = config["cellranger_gex"].get("enrichment", {}).get("mem_gb", 16) * 1024,
            runtime = config["cellranger_gex"].get("enrichment", {}).get("runtime_minutes", 60),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_gex_enrichment.log")
        script:
            "../scripts/merge_metadata.py"


if config.get("cellranger_atac"):

    rule enrich_atac_metadata:
        """Enrich ATAC batch object with analysis metadata (demux, doublet, annotation)."""
        input:
            batch_object = os.path.join(BATCH_DIR, "{batch}_atac.h5ad"),
            done_flag = os.path.join(LOGS_DIR, "{batch}_atac_batch_aggregation.done")
        output:
            enriched_object = os.path.join(FINAL_DIR, "{batch}_atac.h5ad"),
            obs_table = os.path.join(FINAL_DIR, "{batch}_atac_obs.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_atac_enrichment.done"))
        params:
            batch = "{batch}",
            modality = "atac",
            demux_dir = DEMUX_DIR,
            doublet_dir = DOUBLET_DIR
        threads: config["cellranger_atac"].get("enrichment", {}).get("threads", 4)
        resources:
            mem_mb = config["cellranger_atac"].get("enrichment", {}).get("mem_gb", 16) * 1024,
            runtime = config["cellranger_atac"].get("enrichment", {}).get("runtime_minutes", 60),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_atac_enrichment.log")
        script:
            "../scripts/merge_metadata.py"


if config.get("cellranger_arc"):

    rule enrich_arc_metadata:
        """Enrich ARC batch object with analysis metadata (demux, doublet, annotation)."""
        input:
            batch_object = os.path.join(BATCH_DIR, "{batch}_arc.h5mu"),
            done_flag = os.path.join(LOGS_DIR, "{batch}_arc_batch_aggregation.done")
        output:
            enriched_object = os.path.join(FINAL_DIR, "{batch}_arc.h5mu"),
            obs_table = os.path.join(FINAL_DIR, "{batch}_arc_obs.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_arc_enrichment.done"))
        params:
            batch = "{batch}",
            modality = "arc",
            demux_dir = DEMUX_DIR,
            doublet_dir = DOUBLET_DIR
        threads: config["cellranger_arc"].get("enrichment", {}).get("threads", 4)
        resources:
            mem_mb = config["cellranger_arc"].get("enrichment", {}).get("mem_gb", 16) * 1024,
            runtime = config["cellranger_arc"].get("enrichment", {}).get("runtime_minutes", 60),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_arc_enrichment.log")
        script:
            "../scripts/merge_metadata.py"


if config.get("cellranger_multi"):

    rule enrich_multi_metadata:
        """Enrich multi batch object with analysis metadata (demux, doublet, annotation)."""
        input:
            batch_object = os.path.join(BATCH_DIR, "{batch}_multi.h5mu"),
            done_flag = os.path.join(LOGS_DIR, "{batch}_multi_batch_aggregation.done")
        output:
            enriched_object = os.path.join(FINAL_DIR, "{batch}_multi.h5mu"),
            obs_table = os.path.join(FINAL_DIR, "{batch}_multi_obs.tsv.gz"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_multi_enrichment.done"))
        params:
            batch = "{batch}",
            modality = "multi",
            demux_dir = DEMUX_DIR,
            doublet_dir = DOUBLET_DIR
        threads: config["cellranger_multi"].get("enrichment", {}).get("threads", 4)
        resources:
            mem_mb = config["cellranger_multi"].get("enrichment", {}).get("mem_gb", 16) * 1024,
            runtime = config["cellranger_multi"].get("enrichment", {}).get("runtime_minutes", 60),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_multi_enrichment.log")
        script:
            "../scripts/merge_metadata.py"
