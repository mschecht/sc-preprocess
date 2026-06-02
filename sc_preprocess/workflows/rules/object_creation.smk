"""Object creation rules - convert Cell Ranger outputs to AnnData/MuData objects."""

import os
import sys
from pathlib import Path
from tempfile import gettempdir

sys.path.insert(0, str(Path(workflow.basedir).parent / "utils"))
from custom_logger import custom_logger
from sc_preprocess.config_validator import parse_output_directories

# Get centralized output directories
OUTPUT_DIRS = parse_output_directories(config)
ANNDATA_DIR = OUTPUT_DIRS["anndata_dir"]
LOGS_DIR = OUTPUT_DIRS["logs_dir"]


# ============================================================================
# GEX ANNDATA CREATION
# ============================================================================

if config.get("cellranger_gex"):
    GEX_COUNT_DIR = OUTPUT_DIRS["cellrangergex_count_dir"]

    rule create_gex_anndata:
        """Create AnnData object from GEX Cell Ranger count output (per-capture)."""
        input:
            h5 = os.path.join(GEX_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix.h5"),
            count_done = os.path.join(LOGS_DIR, "{batch}_{capture}_gex_count.done")
        output:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.done"))
        params:
            batch = "{batch}",
            capture = "{capture}"
        threads: config["cellranger_gex"].get("anndata", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_gex"].get("anndata", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_gex"].get("anndata", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_gex_anndata.log")
        script:
            "../scripts/create_gex_anndata.py"


# ============================================================================
# ATAC ANNDATA CREATION
# ============================================================================

if config.get("cellranger_atac"):
    ATAC_COUNT_DIR = OUTPUT_DIRS["cellrangeratac_count_dir"]

    rule create_atac_anndata:
        """Create AnnData object from ATAC Cell Ranger count output (per-capture)."""
        input:
            fragments = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "fragments.tsv.gz"),
            peak_matrix = os.path.join(ATAC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_peak_bc_matrix.h5"),
            count_done = os.path.join(LOGS_DIR, "{batch}_{capture}_atac_count.done")
        output:
            h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5ad"),
            snap_h5ad = os.path.join(ANNDATA_DIR, "{batch}_{capture}_snap.h5ad"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_atac_anndata.done"))
        params:
            batch = "{batch}",
            capture = "{capture}"
        threads: config["cellranger_atac"].get("anndata", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_atac"].get("anndata", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_atac"].get("anndata", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_atac_anndata.log")
        script:
            "../scripts/create_atac_anndata.py"


# ============================================================================
# ARC MUDATA CREATION
# ============================================================================

if config.get("cellranger_multi"):
    MULTI_COUNT_DIR = OUTPUT_DIRS["cellrangermulti_count_dir"]

    rule create_multi_mudata:
        """Create MuData object from Cell Ranger multi output (per-capture).

        Handles both non-multiplexed runs (outs/count/) and multiplexed Flex runs
        (outs/per_sample_outs/). VDJ-T and VDJ-B directories are read if present.
        """
        input:
            h5 = os.path.join(MULTI_COUNT_DIR, "{batch}_{capture}", "outs", "count", "filtered_feature_bc_matrix.h5"),
            count_done = os.path.join(LOGS_DIR, "{batch}_{capture}_multi_count.done")
        output:
            h5mu = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5mu"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_multi_mudata.done"))
        params:
            batch = "{batch}",
            capture = "{capture}",
            outs_dir = os.path.join(MULTI_COUNT_DIR, "{batch}_{capture}", "outs"),
            vdj_t_dir = os.path.join(MULTI_COUNT_DIR, "{batch}_{capture}", "outs", "vdj_t"),
            vdj_b_dir = os.path.join(MULTI_COUNT_DIR, "{batch}_{capture}", "outs", "vdj_b")
        threads: config["cellranger_multi"].get("anndata", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_multi"].get("anndata", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_multi"].get("anndata", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_multi_mudata.log")
        script:
            "../scripts/create_multi_mudata.py"


if config.get("cellranger_arc"):
    ARC_COUNT_DIR = OUTPUT_DIRS["cellrangerarc_count_dir"]

    rule create_arc_mudata:
        """Create MuData object from ARC Cell Ranger count output (per-capture)."""
        input:
            barcodes = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
            fragments = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "atac_fragments.tsv.gz"),
            count_done = os.path.join(LOGS_DIR, "{batch}_{capture}_arc_count.done")
        output:
            h5mu = os.path.join(ANNDATA_DIR, "{batch}_{capture}.h5mu"),
            done = touch(os.path.join(LOGS_DIR, "{batch}_{capture}_arc_mudata.done"))
        params:
            batch = "{batch}",
            capture = "{capture}",
            mtx_dir = os.path.join(ARC_COUNT_DIR, "{batch}_{capture}", "outs", "filtered_feature_bc_matrix")
        threads: config["cellranger_arc"].get("anndata", {}).get("threads", 16)
        resources:
            mem_mb = config["cellranger_arc"].get("anndata", {}).get("mem_gb", 32) * 1024,
            runtime = config["cellranger_arc"].get("anndata", {}).get("runtime_minutes", 120),
            tmpdir = RESOURCES.get("tmpdir") or gettempdir()
        log:
            os.path.join(LOGS_DIR, "{batch}_{capture}_arc_mudata.log")
        script:
            "../scripts/create_arc_mudata.py"
